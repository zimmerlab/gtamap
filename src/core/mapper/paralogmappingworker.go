package mapper

import (
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/sirupsen/logrus"
)

func ParalogMappingWorker(needParalogMapChan <-chan *mapperutils.ReadPairMatchResults, wg *sync.WaitGroup, index *index.GenomeIndex, finalMatchesChan chan<- *mapperutils.ReadPairMatchResults) {
	logrus.Info("Started paralog mapping")
	defer wg.Done()

	// per default we want to map all reads to paraog index (except for confident reads)
	for readPairResult := range needParalogMapChan {

		// init paralog map in readPairResult
		readPairResult.ParalogMappings = make(map[string]map[int]*mapperutils.ValidReadPairCombination)
		// projects 0,1 -> 0; 2,3 -> 1; 4,5 -> 2 to map to main target region
		parentSeqIndices := make(map[int]struct{})

		if len(index.Sequences) > 1 {
			// we have several target regions in out index -> find out which one mappedReadPair is mapping to and then
			// only remap on paralog index of that target region
			parentSeqIndices = getParentIndices(readPairResult.Fw, readPairResult.Rv)
		} else {
			// if we only have one seq in our index, we can only remap to its corresponding paralog index
			parentSeqIndices[0] = struct{}{}
		}

		logrus.Infof("Mapping readPair (fwId=%s) to paralog regions of parent indices: %s", readPairResult.ReadPair.ReadR1.Header, parentSeqIndices)

		for parentSeqIndex := range parentSeqIndices {
			parentSeqId := index.GetSequenceInfo(parentSeqIndex).GeneId
			paralogIndex, exists := index.ParalogRegions[parentSeqId]
			if !exists {
				logrus.Debugf("No paralog index provided for target region %s", parentSeqId)
				continue
			}

			// Map only to paralog regions of main region (parentSeqIndex)
			paralogMappingsFw, isMappableFw := MapRead(readPairResult.ReadPair.ReadR1, paralogIndex, true)
			paralogMappingsRv, isMappableRv := MapRead(readPairResult.ReadPair.ReadR2, paralogIndex, true)

			if !isMappableFw || !isMappableRv || len(paralogMappingsFw) == 0 || len(paralogMappingsRv) == 0 {
				logrus.WithFields(logrus.Fields{
					"isMappableFw":  isMappableFw,
					"isMappableRv":  isMappableRv,
					"num resultsFw": len(paralogMappingsFw),
					"num resultsRv": len(paralogMappingsRv),
				}).Debugf("readpair not mappable to paralog index of target region %s", parentSeqId)
				continue
			}
			bestCandidate := getParalogCandidates(paralogMappingsFw, paralogMappingsRv)
			if bestCandidate != nil {
				readPairResult.ParalogMappings[parentSeqId] = bestCandidate
			}
		}
		finalMatchesChan <- readPairResult
	}

	logrus.Debug("Finished paralog mapping")
}

// Returns a slice of parentSeqIndices for a given readPairResult. If we map only to one main target region,
// there is only one paralog index (and we do not need to do anything), but if the index holds more than one seq, we do not want to map to all possible
// indices, only to the paralogs of the parentSeqIndex
func getParentIndices(fwMatches []*mapperutils.ReadMatchResult, rvMatches []*mapperutils.ReadMatchResult) map[int]struct{} {
	mapPerParalogSeqIndex := make(map[int]struct{})

	// accumulate mapping results per paralog main id
	for _, mapping := range fwMatches {
		mapPerParalogSeqIndex[mapping.SequenceIndex/2] = struct{}{}
	}

	for _, mapping := range rvMatches {
		mapPerParalogSeqIndex[mapping.SequenceIndex/2] = struct{}{}
	}
	return mapPerParalogSeqIndex
}

// thre results from the paralog greedy map need to be binned and sorted. We only want to look at results that make sense
// (fw needs corresponding rv map on rv seq etc).
func getParalogCandidates(alternativeMappingsFw []*mapperutils.ReadMatchResult, alternativeMappingsRv []*mapperutils.ReadMatchResult) map[int]*mapperutils.ValidReadPairCombination {
	fwMapPerParalogSeqIndex := make(map[int][]*mapperutils.ReadMatchResult)
	rvMapPerParalogSeqIndex := make(map[int][]*mapperutils.ReadMatchResult)
	mappedParalogIds := make(map[int]struct{})

	// accumulate mapping results per paralog main id
	for _, mapping := range alternativeMappingsFw {
		paralogRegionIndex := mapping.SequenceIndex / 2 // maps back to the main sequence index
		fwMapPerParalogSeqIndex[paralogRegionIndex] = append(fwMapPerParalogSeqIndex[paralogRegionIndex], mapping)
		mappedParalogIds[paralogRegionIndex] = struct{}{}
	}

	for _, mapping := range alternativeMappingsRv {
		paralogRegionIndex := mapping.SequenceIndex / 2
		rvMapPerParalogSeqIndex[paralogRegionIndex] = append(rvMapPerParalogSeqIndex[paralogRegionIndex], mapping)
		mappedParalogIds[paralogRegionIndex] = struct{}{}
	}

	alternativeMappinsInParalogs := make(map[int]*mapperutils.ValidReadPairCombination)
	for paralogId := range mappedParalogIds {
		fwMapsOfParalogId := fwMapPerParalogSeqIndex[paralogId]
		rvMapsOfParalogId := rvMapPerParalogSeqIndex[paralogId]
		alternativeMappinsInParalogs[paralogId] = getPossibleMappingCombinations(fwMapsOfParalogId, rvMapsOfParalogId)
	}
	return alternativeMappinsInParalogs
}

// receives fw and rv matches of one seqID. Returns possible combinations of fw/rv.
// E.g. fw -> 25 and rv -> 25 doesnt work since they need to map to separate strands etc
// Currently returns combination with least amount of mm
func getPossibleMappingCombinations(fwMatches []*mapperutils.ReadMatchResult, rvMatches []*mapperutils.ReadMatchResult) *mapperutils.ValidReadPairCombination {
	var bestCombination *mapperutils.ValidReadPairCombination
	minMisMatches := 100000 // init

	for i := 0; i < len(fwMatches); i++ {
		fwMatch := fwMatches[i]
		for j := 0; j < len(rvMatches); j++ {
			rvMatch := rvMatches[j]
			if fwMatch.SequenceIndex-1 == rvMatch.SequenceIndex || fwMatch.SequenceIndex == rvMatch.SequenceIndex-1 {
				currMM := len(fwMatch.MismatchesRead) + len(rvMatch.MismatchesRead)
				if currMM < minMisMatches {
					minMisMatches = currMM
					if bestCombination == nil {
						bestCombination = &mapperutils.ValidReadPairCombination{
							Fw:            fwMatch,
							Rv:            rvMatch,
							NumMismatches: currMM,
						}
					} else {
						minMisMatches = currMM
						bestCombination.Fw = fwMatch
						bestCombination.Rv = rvMatch
						bestCombination.NumMismatches = currMM
					}
				}
			}
		}
	}
	return bestCombination
}
