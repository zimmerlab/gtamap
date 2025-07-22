package mapper

import (
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/core/mapper/thirdpass"
	"github.com/sirupsen/logrus"
)

func ParalogMappingWorker(needParalogMapChan <-chan *mapperutils.ReadPairMatchResults, wg *sync.WaitGroup, index *index.GenomeIndex,
	thirdPassChan *thirdpass.ThirdPassChannel,
) {
	logrus.Info("Started paralog mapping")
	defer wg.Done()

	// per default we want to map all reads to paraog index (except for confident reads)
	for readPairResult := range needParalogMapChan {

		// init
		paralogMappings := make(map[string]map[int]*mapperutils.ValidReadPairCombination)

		// projects 0,1 -> 0; 2,3 -> 1; 4,5 -> 2 to map to main target region
		parentSeqIndices := make(map[int]struct{})

		// Keeping logic for now (currently we agreed on only having one main target region)
		if len(index.Sequences) > 1 {
			// if we have several target regions in our index -> find out which ones the RP is mapping to and then
			// only remap on paralog index of these target regions
			parentSeqIndices = getParentIndices(readPairResult.Fw, readPairResult.Rv)
		} else {
			// if we only have one seq in our index, we can only remap to its corresponding paralog index
			parentSeqIndices[0] = struct{}{} // if only one seq in main index -> seqIndex == 0
		}

		logrus.Debugf("Mapping readPair (fwId=%s) to paralog regions of parent indices: %s", readPairResult.ReadPair.ReadR1.Header, parentSeqIndices)

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
			bestCandidates := getParalogCandidates(paralogMappingsFw, paralogMappingsRv)
			if len(bestCandidates) != 0 {
				// postprocess best candidates and add to ParalogMappings
				for _, mapping := range bestCandidates {
					postprocessReadMatch(index.ParalogRegions[parentSeqId], readPairResult.ReadPair.ReadR1, mapping.Fw)
					postprocessReadMatch(index.ParalogRegions[parentSeqId], readPairResult.ReadPair.ReadR2, mapping.Rv)
				}
				paralogMappings[parentSeqId] = bestCandidates
			}
		}
		thidPassTask := thirdpass.ThirdPassTask{
			ReadPairId: readPairResult.ReadPair.ReadR1.Header,
			//ParalogInfo: paralogMappings,
		}
		thirdPassChan.Send(&thidPassTask)
	}

	logrus.Info("Done with paralog mapping")
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

// the results from the paralog greedy map need to be binned and sorted. We only want to look at results that make sense
// (fw needs corresponding rv map on rv seq etc).
func getParalogCandidates(alternativeMappingsFw []*mapperutils.ReadMatchResult, alternativeMappingsRv []*mapperutils.ReadMatchResult) map[int]*mapperutils.ValidReadPairCombination {
	fwMapPerParalogSeqIndex, rvMapPerParalogSeqIndex, mappedParalogIds := mapperutils.AssignReadMatchResults(alternativeMappingsFw, alternativeMappingsRv)

	bestParalogMappings := make(map[int]*mapperutils.ValidReadPairCombination)
	for paralogId := range mappedParalogIds {
		fwMapsOfParalogId := fwMapPerParalogSeqIndex[paralogId]
		rvMapsOfParalogId := rvMapPerParalogSeqIndex[paralogId]
		bestCombination := mapperutils.GetBestPossibleMappingCombination(fwMapsOfParalogId, rvMapsOfParalogId)
		if bestCombination != nil {
			bestParalogMappings[paralogId] = bestCombination
		}
	}
	return bestParalogMappings
}
