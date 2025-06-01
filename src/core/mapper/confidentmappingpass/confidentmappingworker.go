package confidentmappingpass

import (
	"fmt"
	"strconv"
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/sirupsen/logrus"
)

func ConfidentMappingWorker(confidentChan *ConfidentPassChan, wgConfidentMapping *sync.WaitGroup, annotationChan chan<- *mapperutils.TargetAnnotation, index *index.GenomeIndex) {
	defer wgConfidentMapping.Done()

	annotation := &mapperutils.TargetAnnotation{
		PreferedStrand: 100,
	}
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

	logrus.Info("Finished collecting all confident maps")
	for seqId, cMaps := range cMapsPerSeq {
		logrus.Infof("%s confident maps for taregt region %s", strconv.Itoa(len(cMaps)), strconv.Itoa(seqId))
	}

	// I. Get Introns
	InferIntrons(cMapsPerSeq, index)

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

func InferIntrons(confMapsPerMainSeqId map[int][]*ConfidentTask, index *index.GenomeIndex) {
	// I. get all gaps in fw direction and mirrow rv to fw dircetion
	fwOrientatedGapsPerMainSeqId := getGapsFwOrientation(confMapsPerMainSeqId, index) // map[int][]*interval.Interval
	// II. count how often each gap exists
	countedGapsPerSeq := countGaps(fwOrientatedGapsPerMainSeqId)
	// III. cluster gaps and get start/stop with highest evidence (eStart|eStop) NOTE: These introns match the seq orientation of the fw mappings
	fwOrientatedIntronsPerSeqId := clusterGaps(countedGapsPerSeq)
	// IV. now we mirror these fw orientated introns to get rv coords
	rvOrientatedIntronsPerSeqId := invertIntrons(fwOrientatedIntronsPerSeqId, index)
	fmt.Println(fwOrientatedIntronsPerSeqId)
	fmt.Println(rvOrientatedIntronsPerSeqId)
}

func invertIntrons(intronsPerSeq map[int][]*regionvector.Region, index *index.GenomeIndex) map[int][]*regionvector.Region {
	mirroredIntronsPerSeqId := make(map[int][]*regionvector.Region)
	for mainSid, introns := range intronsPerSeq {
		geneLength := int(index.GetSequenceInfo(mainSid).EndGenomic - index.GetSequenceInfo(mainSid).StartGenomic)
		for _, intron := range introns {
			mirroredIntronsPerSeqId[mainSid] = append(mirroredIntronsPerSeqId[mainSid], &regionvector.Region{
				// Start: geneLength - intron.End + 2*int(index.SequenceInfo[mainSid].StartGenomic),
				// End:   geneLength - intron.Start + 2*int(index.SequenceInfo[mainSid].StartGenomic), // since end exclusive
				Start: geneLength - intron.End + 2*int(index.SequenceInfo[mainSid].StartGenomic),
				End:   geneLength - intron.Start + 2*int(index.SequenceInfo[mainSid].StartGenomic), // since end exclusive
			})
		}
	}

	return mirroredIntronsPerSeqId
}

// extracts all gaps in fw and rv mappings but mirrors the rv coords to match fw orientation
func getGapsFwOrientation(cMapsPerSeq map[int][]*ConfidentTask, index *index.GenomeIndex) map[int][]*regionvector.Region { // TODO
	fwOrientatedGapsPerMainSeqId := make(map[int][]*regionvector.Region)

	// iterate over confMaps of target seqs
	for sId, confMaps := range cMapsPerSeq {
		gapsInFw := make([]*regionvector.Region, 0)
		gapsInRv := make([]*regionvector.Region, 0)
		// iterate over each confMap in target seq seqId
		for _, confMap := range confMaps {
			// iterate over the fw matches (pairwise) and check if there is a gap
			for i := 0; i < len(confMap.ResultFw.MatchedGenome.Regions)-1; i++ {
				lStop := confMap.ResultFw.MatchedGenome.Regions[i].End
				rStart := confMap.ResultFw.MatchedGenome.Regions[i+1].Start
				if rStart > lStop {
					gapsInFw = append(gapsInFw, &regionvector.Region{
						// Start: lStop + 1,
						// End:   rStart - 1,
						Start: lStop + int(index.SequenceInfo[sId].StartGenomic),
						End:   rStart + int(index.SequenceInfo[sId].StartGenomic),
					})
				}
			}

			// the reverse gaps get mirrored to match fw coords
			for i := 0; i < len(confMap.ResultRv.MatchedGenome.Regions)-1; i++ {
				lStop := confMap.ResultRv.MatchedGenome.Regions[i].End
				rStart := confMap.ResultRv.MatchedGenome.Regions[i+1].Start
				if rStart > lStop {
					gapsInRv = append(gapsInRv, &regionvector.Region{
						// Start: len(*index.Sequences[sId]) - rStart,
						// End:   len(*index.Sequences[sId]) - lStop,
						Start: len(*index.Sequences[sId]) - rStart + int(index.SequenceInfo[sId].StartGenomic),
						End:   len(*index.Sequences[sId]) - lStop + int(index.SequenceInfo[sId].StartGenomic),
					})
				}
			}
		}
		// the gaps stored for the main sId in fw orientation
		// this means that all introns infered from these gaps will match the gene orientation of the fw map
		// we later can mirror the infered introns and that way get the coords for the complement gene strand
		fwOrientatedGapsPerMainSeqId[sId] = append(gapsInFw, gapsInRv...)
	}

	return fwOrientatedGapsPerMainSeqId
}

func countGaps(gapsPerSeq map[int][]*regionvector.Region) map[int]map[regionvector.Region]int {
	gapMap := make(map[int]map[regionvector.Region]int)
	// count unique gapsPerSeq per sid
	for sid, gaps := range gapsPerSeq {
		for _, g := range gaps {
			// check if we already have seen the sid
			_, exists := gapMap[sid]
			if !exists {
				// if not, init new map
				gapMap[sid] = make(map[regionvector.Region]int)
			}
			// check if we have already seen the curr gap in sid
			gap := regionvector.Region{Start: g.Start, End: g.End}
			_, exists = gapMap[sid][gap]
			if !exists {
				// init new count with 1
				gapMap[sid][gap] = 1
			} else {
				// if the curr gap is already in gapMap[sid], increment counter
				gapMap[sid][gap] = gapMap[sid][gap] + 1
			}
		}
	}
	return gapMap
}

type IntronCluster struct {
	lStart      int // left most coord
	rStop       int // right most coord
	maxEvidence int // how many gaps had eStart and eStop
	eStart      int // start of gap with max evidence
	eStop       int // stop of gap with max evidence
}

func clusterGaps(gapMapPerSeq map[int]map[regionvector.Region]int) map[int][]*regionvector.Region {
	// cluster overlapping gaps
	clusters := make(map[int][]*IntronCluster, 0)
	for sid, gapMap := range gapMapPerSeq {
		if _, ok := clusters[sid]; !ok {
			clusters[sid] = []*IntronCluster{}
		}
		for currGap, currGapEvidence := range gapMap {
			if len(clusters[sid]) == 0 {
				// init slice with first cluster
				cluster := &IntronCluster{
					lStart:      currGap.Start,
					rStop:       currGap.End,
					maxEvidence: currGapEvidence,
					eStart:      currGap.Start,
					eStop:       currGap.End,
				}
				clusters[sid] = append(clusters[sid], cluster)
				continue
			}
			// keep track if curr gap was added
			addedToExistingCluster := false

			for _, cluster := range clusters[sid] {
				// is the currGap overlapping with the cluster
				// I. -----##########------- cluster
				//    ------######---------- currGap → contained by cluster, add to cluster
				if currGap.Start >= cluster.lStart && currGap.End <= cluster.rStop {
					if cluster.maxEvidence < currGapEvidence {
						cluster.maxEvidence = currGapEvidence
						cluster.eStart = currGap.Start
						cluster.eStop = currGap.End
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
					}
					cluster.rStop = currGap.End
					addedToExistingCluster = true
					break
				}
			}

			// if the gap wasn't added to any existing cluster, create a new one
			if !addedToExistingCluster {
				newCluster := &IntronCluster{
					lStart:      currGap.Start,
					rStop:       currGap.End,
					maxEvidence: currGapEvidence,
					eStart:      currGap.Start,
					eStop:       currGap.End,
				}
				clusters[sid] = append(clusters[sid], newCluster)
			}
		}
	}

	// now extract introns with highest evidence
	fwOrientatedGapsPerMainSeqId := make(map[int][]*regionvector.Region)
	for mainSid, intronClusters := range clusters {
		for _, intronCluster := range intronClusters {
			fwOrientatedGapsPerMainSeqId[mainSid] = append(fwOrientatedGapsPerMainSeqId[mainSid], &regionvector.Region{
				Start: intronCluster.eStart,
				End:   intronCluster.eStop,
			})
		}
	}
	return fwOrientatedGapsPerMainSeqId
}
