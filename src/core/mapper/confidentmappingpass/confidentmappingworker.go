package confidentmappingpass

import (
	"fmt"
	"os"
	"sort"
	"strconv"
	"sync"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/datastructure"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
)

func ConfidentMappingWorker(
	confidentChan *ConfidentPassChan,
	wgConfidentMapping *sync.WaitGroup,
	annotationChan chan<- map[int]*mapperutils.TargetAnnotation,
	index *index.GenomeIndex,
) {

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

		if len(annotation[targetId].Introns[0].Regions) == 0 {
			annotation[targetId] = nil // if no junctions in conf maps, no introns can be inferred and no graph can be created
		} else {
			// build graphs
			annotation[targetId].Introns[targetId].BuildTranscriptomeGraph(len(*index.Sequences[targetId]))
			annotation[targetId].Introns[targetId+1].BuildTranscriptomeGraph(len(*index.Sequences[targetId]))
			annotation[targetId].LogInfo()
		}
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

func repairIntrons(
	inferredIntrons []*regionvector.Intron,
	targetId int,
	genomeIndex *index.GenomeIndex,
) []*regionvector.Intron {

	repairedIntrons := make([]*regionvector.Intron, 0)

	// if IntronClusterRepairWindow == 10
	// we go from -5 to 5 (starting at i.Start and i.End)

	windowDelta := config.Mapper.GetConfidentIntronClusterRepairWindow() / 2
	refSeq := genomeIndex.Sequences[targetId]

	for _, intron := range inferredIntrons {
		if intron.SpliceSiteScore > 0 {
			repairedIntrons = append(repairedIntrons, intron)
			continue
		}

		var lookOnPlusStrand bool
		if genomeIndex.GetSequenceInfo(targetId).IsForwardStrand {
			if genomeIndex.IsSequenceForward(targetId) {
				lookOnPlusStrand = true
			} else {
				lookOnPlusStrand = false
			}
		} else {
			if genomeIndex.IsSequenceForward(targetId) {
				lookOnPlusStrand = false
			} else {
				lookOnPlusStrand = true
			}
		}

		repaired := false
		donorSiteStart := intron.Start
		acceptorSiteStart := intron.End
		for i := -windowDelta; i < windowDelta; i++ {
			donorSiteSeq := (*refSeq)[donorSiteStart+i : donorSiteStart+i+2]
			for j := -windowDelta; j < windowDelta; j++ {
				acceptorSiteSeq := (*refSeq)[acceptorSiteStart+j-2 : acceptorSiteStart+j]

				// if fw -> true else false
				// this can be done more efficiently if i split method in two
				score, isKnownSpliceSite := utils.ScoreSpliceSites(donorSiteSeq[0], donorSiteSeq[1],
					acceptorSiteSeq[0], acceptorSiteSeq[1], lookOnPlusStrand)
				if isKnownSpliceSite && intron.Start+i < intron.End+j-1 {
					repairedIntrons = append(repairedIntrons, &regionvector.Intron{
						Start:           intron.Start + i,
						End:             intron.End + j,
						Evidence:        intron.Evidence,
						Rank:            intron.Rank,
						SpliceSiteScore: 2 - score, // ScoreSpliceSites returns 0 for canonical splice site -> 2 = canonical, 1 = non-can, 0 = no splice site
					})
					repaired = true
					break
				}
				if repaired {
					break
				}
			}

		}
		// we checked all positions but could not repair intron: append unrepaired intron
		if !repaired {
			repairedIntrons = append(repairedIntrons, intron)
		}

	}
	return repairedIntrons
}

func InferIntronsOfTarget(
	targetId int,
	confMaps []*ConfidentTask,
	index *index.GenomeIndex,
) *mapperutils.TargetAnnotation {

	// I. get all gaps of targetId in plus orientation
	plusOrientatedGaps, plusEvidence, minusEvidence := getGapsPlusOrientation(
		targetId,
		confMaps,
		index,
	) // map[int][]*interval.Interval

	// II. count how often each gap exists (plus orientation gaps)
	countedGapsOfTarget := countGaps(plusOrientatedGaps)

	// III. cluster plus oriented gaps and get start/stop with highest evidence (eStart|eStop)
	plusOrientatedIntronsOfTarget := clusterGaps(countedGapsOfTarget)

	// Check if some introns which don't follow splice sites would follow them by modifying start/end coords
	if config.Mapper.Mapping.IsReadOriginRna {
		repairedIntrons := repairIntrons(plusOrientatedIntronsOfTarget, targetId, index)
		if repairedIntrons != nil {
			plusOrientatedIntronsOfTarget = repairedIntrons
		}
	}

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

	intronTreeMap := make(map[int]*datastructure.INTree)
	convertedPlusRegions := make([]datastructure.Bounds, len(plusOrientatedIntronsOfTarget))
	for i, intron := range plusOrientatedIntronsOfTarget {
		convertedPlusRegions[i] = &datastructure.IntronRegion{
			Lower: intron.Start,
			Upper: intron.End,
		}
	}
	convertedMinusRegions := make([]datastructure.Bounds, len(minusOrientatedIntronsOfTarget))
	for i, intron := range minusOrientatedIntronsOfTarget {
		convertedMinusRegions[i] = &datastructure.IntronRegion{
			Lower: intron.Start,
			Upper: intron.End,
		}
	}

	intronTreeMap[0] = datastructure.NewINTree(convertedPlusRegions)
	intronTreeMap[1] = datastructure.NewINTree(convertedMinusRegions)

	targetAnnotation := mapperutils.TargetAnnotation{
		PreferedStrand: preferredStrandedness,
		Confidence:     confidence,
		Introns:        intronMap,
		IntronTrees:    intronTreeMap,
	}

	return &targetAnnotation
}

func invertIntrons(
	sId int,
	intronsOfTarget []*regionvector.Intron,
	index *index.GenomeIndex,
) []*regionvector.Intron {

	mirroredIntronsPerSeqId := make([]*regionvector.Intron, 0)

	geneLength := int(
		index.GetSequenceInfo(sId).EndGenomic -
			index.GetSequenceInfo(sId).StartGenomic,
	)

	for _, intron := range intronsOfTarget {
		mirroredIntronsPerSeqId = append(
			mirroredIntronsPerSeqId,
			&regionvector.Intron{
				// Start:    geneLength - intron.End + 2*int(index.SequenceInfo[sId].StartGenomic),
				// End:      geneLength - intron.Start + 2*int(index.SequenceInfo[sId].StartGenomic),
				Start:           geneLength - intron.End,
				End:             geneLength - intron.Start,
				Evidence:        intron.Evidence,
				SpliceSiteScore: intron.SpliceSiteScore,
			},
		)
	}

	return mirroredIntronsPerSeqId
}

// extracts all gaps in fw and rv mappings but mirrors the rv coords to match fw orientation
func getGapsPlusOrientation(
	sId int,
	cMapsPerSeq []*ConfidentTask,
	index *index.GenomeIndex,
) (
	[]*regionvector.Gap,
	int,
	int,
) {

	plusOrientatedGapsPerMainSeqId := make([]*regionvector.Gap, 0)
	plusStrandednessEvidence := 0
	minusStrandednessEvidence := 0

	// iterate over each confMap in target seq seqId
	for _, confMap := range cMapsPerSeq {

		if len(confMap.ResultFw.MatchedRead.Regions) != len(confMap.ResultFw.MatchedGenome.Regions) {
			fmt.Println("error in get gaps plus orientation result fw")
			os.Exit(1)
		}
		if len(confMap.ResultRv.MatchedRead.Regions) != len(confMap.ResultRv.MatchedGenome.Regions) {
			fmt.Println("error in get gaps plus orientation result rv")
			os.Exit(1)
		}

		confMap.ResultFw.NormalizeRegions()
		confMap.ResultRv.NormalizeRegions()

		// update strandedness counter
		if confMap.ResultFw.SequenceIndex%2 == 0 {
			plusStrandednessEvidence++
		} else {
			minusStrandednessEvidence++
		}

		skippedGaps := 0
		// iterate over the fw matches (pairwise) and check if there is a gap
		for i := 0; i < len(confMap.ResultFw.MatchedGenome.Regions)-1; i++ {
			lStop := confMap.ResultFw.MatchedGenome.Regions[i].End
			rStart := confMap.ResultFw.MatchedGenome.Regions[i+1].Start
			if lStop == rStart {
				skippedGaps++
				continue // no gap in genome, only in read
			}

			var knownSpliceSite int = 0

			if config.Mapper.Mapping.IsReadOriginRna &&
				len(confMap.ResultFw.SpliceSitesInfo) != 0 {
				knownSpliceSite = confMap.ResultFw.SpliceSitesInfo[i-skippedGaps]
			}

			if rStart > lStop {
				if confMap.ResultFw.SequenceIndex%2 == 0 {
					plusOrientatedGapsPerMainSeqId = append(plusOrientatedGapsPerMainSeqId, &regionvector.Gap{
						Start:           lStop,
						End:             rStart,
						SpliceSiteScore: knownSpliceSite,
					})
				} else {
					plusOrientatedGapsPerMainSeqId = append(plusOrientatedGapsPerMainSeqId, &regionvector.Gap{
						Start:           len(*index.Sequences[sId]) - rStart,
						End:             len(*index.Sequences[sId]) - lStop,
						SpliceSiteScore: knownSpliceSite,
					})
				}
			}
		}

		// the reverse gaps get mirrored to match fw coords
		skippedGaps = 0
		for i := 0; i < len(confMap.ResultRv.MatchedGenome.Regions)-1; i++ {
			lStop := confMap.ResultRv.MatchedGenome.Regions[i].End
			rStart := confMap.ResultRv.MatchedGenome.Regions[i+1].Start

			if lStop == rStart {
				skippedGaps++
				continue // no gap in genome, only in read
			}

			var knownSpliceSite int = 0
			if config.Mapper.Mapping.IsReadOriginRna &&
				len(confMap.ResultRv.SpliceSitesInfo) != 0 {
				knownSpliceSite = confMap.ResultRv.SpliceSitesInfo[i-skippedGaps]
			}

			if rStart > lStop {
				if confMap.ResultRv.SequenceIndex%2 == 0 {
					plusOrientatedGapsPerMainSeqId = append(
						plusOrientatedGapsPerMainSeqId,
						&regionvector.Gap{
							Start:           lStop,
							End:             rStart,
							SpliceSiteScore: knownSpliceSite,
						},
					)
				} else {
					plusOrientatedGapsPerMainSeqId = append(
						plusOrientatedGapsPerMainSeqId,
						&regionvector.Gap{
							Start:           len(*index.Sequences[sId]) - rStart,
							End:             len(*index.Sequences[sId]) - lStop,
							SpliceSiteScore: knownSpliceSite,
						},
					)
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

		gapMap[*g] = gapMap[*g] + g.SpliceSiteScore // reward splice sites
	}
	return gapMap
}

type IntronCluster struct {
	lStart                     int // left most coord
	rStop                      int // right most coord
	maxEvidence                int // how many gaps had eStart and eStop
	eStart                     int // start of gap with max evidence
	eStop                      int // stop of gap with max evidence
	maxEvidenceSpliceSiteScore int
}

func clusterGaps(targetGaps map[regionvector.Gap]int) []*regionvector.Intron {
	// init slice
	sortedGaps := make([]regionvector.Gap, len(targetGaps))
	ind := 0
	for gap := range targetGaps {
		sortedGaps[ind] = gap
		ind++
	}
	// sort by length s.t. we cluster bottom up
	sort.Slice(sortedGaps, func(i, j int) bool {
		if sortedGaps[i].Start == sortedGaps[j].Start || sortedGaps[i].End == sortedGaps[j].End {
			return sortedGaps[i].Length() < sortedGaps[j].Length()
		}

		if sortedGaps[i].Length() == sortedGaps[j].Length() {
			return sortedGaps[i].Start < sortedGaps[j].Start
		}
		return sortedGaps[i].Length() < sortedGaps[j].Length()
	})

	// cluster overlapping gaps
	clusters := make([]*IntronCluster, 0)

	// iterate in order
	for _, currGap := range sortedGaps {
		currGapEvidence := targetGaps[currGap]
		if len(clusters) == 0 {
			// init slice with first cluster
			cluster := &IntronCluster{
				lStart:                     currGap.Start,
				rStop:                      currGap.End,
				maxEvidence:                currGapEvidence,
				eStart:                     currGap.Start,
				eStop:                      currGap.End,
				maxEvidenceSpliceSiteScore: currGap.SpliceSiteScore,
			}
			clusters = append(clusters, cluster)
			continue
		}
		// keep track if curr gap was added
		addedToExistingCluster := false

		for _, cluster := range clusters {
			if currGap.End > cluster.lStart && currGap.Start < cluster.rStop {

				// check if dist from gap start to cluster start is small enough to absorb
				if cluster.lStart > currGap.Start &&
					cluster.lStart-currGap.Start > config.Mapper.GetConfidentIntronClusterDelta() {
					continue
				}
				// check if dist from gap end to cluster end is small enough to absorb
				if cluster.rStop < currGap.End &&
					currGap.End-cluster.rStop > config.Mapper.GetConfidentIntronClusterDelta() {
					continue
				}
				if cluster.maxEvidence < currGapEvidence {
					cluster.maxEvidence = currGapEvidence
					cluster.eStart = currGap.Start
					cluster.eStop = currGap.End
					cluster.maxEvidenceSpliceSiteScore = currGap.SpliceSiteScore
				}
				// grow cluster
				cluster.rStop = currGap.End
				cluster.lStart = currGap.Start
				addedToExistingCluster = true
				break
			}
		}

		// if the gap wasn't added to any existing cluster, create a new one
		if !addedToExistingCluster {
			newCluster := &IntronCluster{
				lStart:                     currGap.Start,
				rStop:                      currGap.End,
				maxEvidence:                currGapEvidence,
				eStart:                     currGap.Start,
				eStop:                      currGap.End,
				maxEvidenceSpliceSiteScore: currGap.SpliceSiteScore,
			}
			clusters = append(clusters, newCluster)
		}
	}

	// now extract introns with highest evidence
	fwOrientatedGapsOfTarget := make([]*regionvector.Intron, 0)
	for _, intronCluster := range clusters {
		fwOrientatedGapsOfTarget = append(fwOrientatedGapsOfTarget, &regionvector.Intron{
			Start:           intronCluster.eStart,
			End:             intronCluster.eStop,
			Evidence:        intronCluster.maxEvidence,
			SpliceSiteScore: intronCluster.maxEvidenceSpliceSiteScore,
		})
	}
	return fwOrientatedGapsOfTarget
}
