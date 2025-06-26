package confidentmappingpass

import (
	"strconv"
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/sirupsen/logrus"
)

func ConfidentMappingWorker(confidentChan *ConfidentPassChan, wgConfidentMapping *sync.WaitGroup, annotationChan chan<- map[int]*mapperutils.TargetAnnotation, index *index.GenomeIndex) {
	defer wgConfidentMapping.Done()

	cMapsPerSeq := make(map[int][]*ConfidentTask, 0)

	for {
		confidentResult, ok := confidentChan.Receive()
		if !ok {
			break
		}

		// we know that fw and rv belong together, otherwise the rp wouldn't be sent to confidentChan
		targetMainId := confidentResult.ResultFw.SequenceIndex / 2
		cMapsPerSeq[targetMainId] = append(cMapsPerSeq[targetMainId], confidentResult)
	}

	// mainTargetId -> annotation
	// annotation then stores introns in plus orientation and minus orientation
	annotation := make(map[int]*mapperutils.TargetAnnotation)

	logrus.Info("Finished collecting all confident maps")
	for targetId, cMaps := range cMapsPerSeq {
		logrus.Infof("%s confident maps for target region %s", strconv.Itoa(len(cMaps)), strconv.Itoa(targetId))
		// get Introns per seqId
		// NOTE: Introns are 0 based, start inclusive and end exclusive
		annotation[targetId] = InferIntronsOfTarget(targetId, cMaps, index)
	}

	annotationChan <- annotation
	logrus.Info("Done with Annotation")
}

func getCoverageSlice(mappedReadPairs []*mapperutils.ReadPairMatchResults, geneLength int, index *index.GenomeIndex) []int {
	coverage := make([]int, geneLength)
	for _, mappedReadPair := range mappedReadPairs {
		// forward reads
		for _, region := range mappedReadPair.Fw[0].MatchedGenome.Regions {
			if index.IsSequenceForward(mappedReadPair.Fw[0].SequenceIndex) {
				for pos := region.Start; pos < region.End; pos++ {
					coverage[pos]++
				}
			} else {
				// reverse strand reads need complementary position
				// as done in FormatMappedReadPairToSAM
				for pos := region.Start; pos < region.End; pos++ {
					complementaryPos := geneLength - pos - 1
					coverage[complementaryPos]++
				}
			}
		}

		// reverse reads
		for _, region := range mappedReadPair.Rv[0].MatchedGenome.Regions {
			if index.IsSequenceForward(mappedReadPair.Rv[0].SequenceIndex) {
				for pos := region.Start; pos < region.End; pos++ {
					coverage[pos]++
				}
			} else {
				// get compl pos
				for pos := region.Start; pos < region.End; pos++ {
					complementaryPos := geneLength - pos - 1
					coverage[complementaryPos]++
				}
			}
		}
	}
	return coverage
}

func getHighCoverageRegons(coverageSlice []int) [][2]int {
	highCov := make([][2]int, 0)
	var deltaOne int
	var deltaTwo int
	foundStart := false
	var start int
	var covAtStart int
	for i := 0; i < len(coverageSlice)-8; i++ {
		deltaOne = coverageSlice[i+1] - coverageSlice[i]
		deltaTwo = coverageSlice[i+8] - coverageSlice[i]

		if !foundStart && deltaOne > 30 && deltaTwo > 200 {
			foundStart = true
			start = i
			covAtStart = coverageSlice[i]
		}

		if foundStart && coverageSlice[i] < covAtStart {
			foundStart = false
			highCov = append(highCov, [2]int{start, i})
		}

	}
	return highCov
}

func getMultimappingReads(mappedReadPairs []*mapperutils.ReadPairMatchResults, fwIntervals [][2]int, rvIntervals [][2]int) []*mapperutils.ReadPairMatchResults {
	multimappings := make([]*mapperutils.ReadPairMatchResults, 0)
main:
	for _, mappedReadPair := range mappedReadPairs {
		// check fw
		for _, interval := range fwIntervals {
			// append multimapp only if #mm > 0
			if spansInterval(mappedReadPair.Fw[0].MatchedGenome, interval) && len(mappedReadPair.Fw[0].MismatchesRead) != 0 {
				multimappings = append(multimappings, mappedReadPair)
				continue main
			}
		}

		// check rv
		for _, interval := range rvIntervals {
			// append multimapp only if #mm > 0
			if spansInterval(mappedReadPair.Rv[0].MatchedGenome, interval) && len(mappedReadPair.Rv[0].MismatchesRead) != 0 {
				multimappings = append(multimappings, mappedReadPair)
				continue main
			}
		}
	}
	return multimappings
}

func spansInterval(regions *regionvector.RegionVector, interval [2]int) bool {
	for _, region := range regions.Regions {
		if (region.Start <= interval[0] && region.End > interval[0]) || // region starts before interval and extends into it
			(region.Start < interval[1] && region.End >= interval[1]) || // region ends after interval and starts within it
			(region.Start >= interval[0] && region.End <= interval[1]) { // region is completely contained within interval
			return true
		}
	}
	return false
}

func getCoverageSlices(mappedReadPairs []*mapperutils.ReadPairMatchResults, geneLength int, index *index.GenomeIndex) ([]int, []int) {
	coverageFw := make([]int, geneLength)
	coverageRv := make([]int, geneLength)

	for _, mappedReadPair := range mappedReadPairs {
		// fw reads
		for _, region := range mappedReadPair.Fw[0].MatchedGenome.Regions {
			if index.IsSequenceForward(mappedReadPair.Fw[0].SequenceIndex) {
				for pos := region.Start; pos < region.End; pos++ {
					coverageFw[pos]++
				}
			} else {
				for pos := region.Start; pos < region.End; pos++ {
					coverageRv[pos]++
				}
			}
		}

		// rv reads
		for _, region := range mappedReadPair.Rv[0].MatchedGenome.Regions {
			if index.IsSequenceForward(mappedReadPair.Rv[0].SequenceIndex) {
				for pos := region.Start; pos < region.End; pos++ {
					coverageFw[pos]++
				}
			} else {
				for pos := region.Start; pos < region.End; pos++ {
					coverageRv[pos]++
				}
			}
		}
	}

	return coverageFw, coverageRv
}

func getTotalCoverage(coverageFw, coverageRv []int) []int {
	// this can be used to get the total coverage of fw and rv by mirrowing rv coords
	totalCoverage := make([]int, len(coverageFw))
	for i := 0; i < len(coverageFw); i++ {
		totalCoverage[i] = coverageFw[i] + coverageRv[len(coverageFw)-1-i]
	}
	return totalCoverage
}

func InferIntronsOfTarget(targetId int, confMaps []*ConfidentTask, index *index.GenomeIndex) *mapperutils.TargetAnnotation {
	// I. get all gaps of targetId in plus orientation
	plusOrientatedGaps, plusEvidence, minusEvidence := getGapsPlusOrientation(targetId, confMaps, index) // map[int][]*interval.Interval
	// II. count how often each gap exists (plus orientation gaps)
	countedGapsOfTarget := countGaps(plusOrientatedGaps)
	// III. cluster plus oriented gaps and get start/stop with highest evidence (eStart|eStop)
	plusOrientatedIntronsOfTarget := clusterGaps(countedGapsOfTarget)
	// IV. now we mirror these plus orientated introns to get minus coords
	minusOrientatedIntronsOfTarget := invertIntrons(targetId, plusOrientatedIntronsOfTarget, index)

	var preferredStrandedness int
	var confidence float32
	if plusEvidence > minusEvidence {
		preferredStrandedness = 0 // fw -> 0 rv -> 1 (plus configuration)
		confidence = float32(plusEvidence) / float32(plusEvidence+minusEvidence)
	} else {
		preferredStrandedness = 1 // fw -> 1 rv -> 0 (minus configuration)
		confidence = float32(minusEvidence) / float32(plusEvidence+minusEvidence)
	}

	if plusEvidence == minusEvidence {
		preferredStrandedness = -1
		confidence = 0.5
	}

	intronMap := make(map[int]*regionvector.RegionSet)
	intronMap[0] = regionvector.NewRegionSet(plusOrientatedIntronsOfTarget)
	intronMap[1] = regionvector.NewRegionSet(minusOrientatedIntronsOfTarget)

	targetAnnotation := mapperutils.TargetAnnotation{
		PreferedStrand: preferredStrandedness,
		Confidence:     confidence,
		Introns:        intronMap,
	}

	return &targetAnnotation
}

func invertIntrons(sId int, intronsOfTarget []*regionvector.Intron, index *index.GenomeIndex) []*regionvector.Intron {
	mirroredIntronsPerSeqId := make([]*regionvector.Intron, 0)
	geneLength := int(index.GetSequenceInfo(sId).EndGenomic - index.GetSequenceInfo(sId).StartGenomic)
	for _, intron := range intronsOfTarget {
		mirroredIntronsPerSeqId = append(mirroredIntronsPerSeqId, &regionvector.Intron{
			// Start:    geneLength - intron.End + 2*int(index.SequenceInfo[sId].StartGenomic),
			// End:      geneLength - intron.Start + 2*int(index.SequenceInfo[sId].StartGenomic),
			Start:          geneLength - intron.End,
			End:            geneLength - intron.Start,
			Evidence:       intron.Evidence,
			TrueSpliceSite: intron.TrueSpliceSite,
		})
	}

	return mirroredIntronsPerSeqId
}

// extracts all gaps in fw and rv mappings but mirrors the rv coords to match fw orientation
func getGapsPlusOrientation(sId int, cMapsPerSeq []*ConfidentTask, index *index.GenomeIndex) ([]*regionvector.Gap, int, int) {
	plusOrientatedGapsPerMainSeqId := make([]*regionvector.Gap, 0)
	plusStrandednessEvidence := 0
	minusStrandednessEvidence := 0

	// iterate over each confMap in target seq seqId
	for _, confMap := range cMapsPerSeq {
		confMap.ResultFw.MatchedGenome.MergeAlignmentBlocks()
		confMap.ResultRv.MatchedGenome.MergeAlignmentBlocks()

		// update strandedness counter
		if confMap.ResultFw.SequenceIndex%2 == 0 {
			plusStrandednessEvidence++
		} else {
			minusStrandednessEvidence++
		}

		// iterate over the fw matches (pairwise) and check if there is a gap
		for i := 0; i < len(confMap.ResultFw.MatchedGenome.Regions)-1; i++ {
			lStop := confMap.ResultFw.MatchedGenome.Regions[i].End
			rStart := confMap.ResultFw.MatchedGenome.Regions[i+1].Start
			knownSpliceSite := confMap.ResultFw.SpliceSitesInfo[i]
			if rStart > lStop {
				if confMap.ResultFw.SequenceIndex%2 == 0 {
					plusOrientatedGapsPerMainSeqId = append(plusOrientatedGapsPerMainSeqId, &regionvector.Gap{
						Start:           lStop,
						End:             rStart,
						KnownSpliceSite: knownSpliceSite,
					})
				} else {
					plusOrientatedGapsPerMainSeqId = append(plusOrientatedGapsPerMainSeqId, &regionvector.Gap{
						Start:           len(*index.Sequences[sId]) - rStart,
						End:             len(*index.Sequences[sId]) - lStop,
						KnownSpliceSite: knownSpliceSite,
					})
				}
			}
		}

		// the reverse gaps get mirrored to match fw coords
		for i := 0; i < len(confMap.ResultRv.MatchedGenome.Regions)-1; i++ {
			lStop := confMap.ResultRv.MatchedGenome.Regions[i].End
			rStart := confMap.ResultRv.MatchedGenome.Regions[i+1].Start
			knownSpliceSite := confMap.ResultRv.SpliceSitesInfo[i]
			if rStart > lStop {
				if confMap.ResultRv.SequenceIndex%2 == 0 {
					plusOrientatedGapsPerMainSeqId = append(plusOrientatedGapsPerMainSeqId, &regionvector.Gap{
						Start:           lStop,
						End:             rStart,
						KnownSpliceSite: knownSpliceSite,
					})
				} else {
					plusOrientatedGapsPerMainSeqId = append(plusOrientatedGapsPerMainSeqId, &regionvector.Gap{
						Start:           len(*index.Sequences[sId]) - rStart,
						End:             len(*index.Sequences[sId]) - lStop,
						KnownSpliceSite: knownSpliceSite,
					})
				}
			}
		}
	}

	return plusOrientatedGapsPerMainSeqId, plusStrandednessEvidence, minusStrandednessEvidence
}

func countGaps(gapsOfTarget []*regionvector.Gap) map[regionvector.Gap]int {
	gapMap := make(map[regionvector.Gap]int)
	// count unique gapsPerSeq per sid
	for _, g := range gapsOfTarget {
		// check if we have already seen the curr gap in sid
		_, exists := gapMap[*g]
		if !exists {
			// init new count with 1
			gapMap[*g] = 1
		} else {
			// increment
			gapMap[*g] = gapMap[*g] + 1
		}

		if g.KnownSpliceSite {
			gapMap[*g] = gapMap[*g] + 10 // reward splice sites
		}
	}
	return gapMap
}

type IntronCluster struct {
	lStart                       int // left most coord
	rStop                        int // right most coord
	maxEvidence                  int // how many gaps had eStart and eStop
	eStart                       int // start of gap with max evidence
	eStop                        int // stop of gap with max evidence
	maxEvidenceFollowsSpliceSite bool
}

func clusterGaps(targetGaps map[regionvector.Gap]int) []*regionvector.Intron {
	// cluster overlapping gaps
	clusters := make([]*IntronCluster, 0)
	for currGap, currGapEvidence := range targetGaps {
		if len(clusters) == 0 {
			// init slice with first cluster
			cluster := &IntronCluster{
				lStart:                       currGap.Start,
				rStop:                        currGap.End,
				maxEvidence:                  currGapEvidence,
				eStart:                       currGap.Start,
				eStop:                        currGap.End,
				maxEvidenceFollowsSpliceSite: currGap.KnownSpliceSite,
			}
			clusters = append(clusters, cluster)
			continue
		}
		// keep track if curr gap was added
		addedToExistingCluster := false

		for _, cluster := range clusters {
			// is the currGap overlapping with the cluster
			// I. -----##########------- cluster
			//    ------######---------- currGap → contained by cluster, add to cluster
			if currGap.Start >= cluster.lStart && currGap.End <= cluster.rStop {
				if cluster.maxEvidence < currGapEvidence {
					cluster.maxEvidence = currGapEvidence
					cluster.eStart = currGap.Start
					cluster.eStop = currGap.End
					cluster.maxEvidenceFollowsSpliceSite = currGap.KnownSpliceSite
				}
				addedToExistingCluster = true
				break
			}
			// II -----##########------- cluster
			//    ---#############------ currGap → contains cluster completely
			if currGap.Start <= cluster.lStart && currGap.End >= cluster.rStop {
				if cluster.maxEvidence < currGapEvidence {
					cluster.maxEvidence = currGapEvidence
					cluster.eStart = currGap.Start
					cluster.eStop = currGap.End
					cluster.lStart = currGap.Start
					cluster.rStop = currGap.End
					cluster.maxEvidenceFollowsSpliceSite = currGap.KnownSpliceSite
				} else {
					cluster.lStart = currGap.Start
					cluster.rStop = currGap.End
				}
				addedToExistingCluster = true
				break
			}
			// III -----##########------- cluster
			//     ---############------  currGap → partially contains cluster completely
			if currGap.Start <= cluster.lStart && currGap.End >= cluster.lStart && currGap.End <= cluster.rStop {
				if cluster.maxEvidence < currGapEvidence {
					cluster.maxEvidence = currGapEvidence
					cluster.eStart = currGap.Start
					cluster.eStop = currGap.End
					cluster.maxEvidenceFollowsSpliceSite = currGap.KnownSpliceSite
				}
				cluster.lStart = currGap.Start
				addedToExistingCluster = true
				break
			}
			// IV  -----##########------- cluster
			//     -----############----  currGap → partially contains cluster completely
			if currGap.Start >= cluster.lStart && currGap.Start <= cluster.rStop && currGap.End >= cluster.rStop {
				if cluster.maxEvidence < currGapEvidence {
					cluster.maxEvidence = currGapEvidence
					cluster.eStart = currGap.Start
					cluster.eStop = currGap.End
					cluster.maxEvidenceFollowsSpliceSite = currGap.KnownSpliceSite
				}
				cluster.rStop = currGap.End
				addedToExistingCluster = true
				break
			}
		}

		// if the gap wasn't added to any existing cluster, create a new one
		if !addedToExistingCluster {
			newCluster := &IntronCluster{
				lStart:                       currGap.Start,
				rStop:                        currGap.End,
				maxEvidence:                  currGapEvidence,
				eStart:                       currGap.Start,
				eStop:                        currGap.End,
				maxEvidenceFollowsSpliceSite: currGap.KnownSpliceSite,
			}
			clusters = append(clusters, newCluster)
		}
	}

	// now extract introns with highest evidence
	fwOrientatedGapsOfTarget := make([]*regionvector.Intron, 0)
	for _, intronCluster := range clusters {
		fwOrientatedGapsOfTarget = append(fwOrientatedGapsOfTarget, &regionvector.Intron{
			Start:          intronCluster.eStart,
			End:            intronCluster.eStop,
			Evidence:       intronCluster.maxEvidence,
			TrueSpliceSite: intronCluster.maxEvidenceFollowsSpliceSite,
		})
	}
	return fwOrientatedGapsOfTarget
}
