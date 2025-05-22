package mapper

import (
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/sirupsen/logrus"
)

func ParalogMappingWorker(readPairChan <-chan *mapperutils.ReadPairMatchResults, wg *sync.WaitGroup, index *index.GenomeIndex) {
	logrus.Info("Started paralog mapping")

	// per default we want to map all reads to paraog index (except for confident reads)
	for readPairResult := range readPairChan {
		for mappedReadPair := range readPairChan {
			// projects 0,1 -> 0; 2,3 -> 1; 4,5 -> 2 to map to main taret seq id
			parentSeqIndices := make(map[int]struct{})

			if len(index.Sequences) > 1 {
				// we have several target regions in out index -> find out which one mappedReadPair is mapping to and then
				// only remap on paralog index of that target region
				parentSeqIndices = getParentIndex(mappedReadPair.Fw, mappedReadPair.Rv)
			} else {
				// if we only have one seq in our index, we can only remap to its corresponding paralog index
				parentSeqIndices[0] = struct{}{}
			}

			logrus.Debugf("Mapping readPair (fwId=%s) to paralog regions of parent indices: %s", readPairResult.ReadPair.ReadR1.Header, parentSeqIndices)
			continue
			for parentSeqIndex := range parentSeqIndices {
				parentSeqId := index.SequenceHeaders[parentSeqIndex]

				// Map only to paralog regions of main region (parentSeqIndex)
				alternativeMappingsFw, isMappableFw := MapRead(readPairResult.ReadPair.ReadR1, index.ParalogRegions[parentSeqId])
				alternativeMappingsRv, isMappableRv := MapRead(readPairResult.ReadPair.ReadR2, index.ParalogRegions[parentSeqId])

				if !isMappableFw || !isMappableRv || len(alternativeMappingsFw) == 0 || len(alternativeMappingsRv) == 0 {
					logrus.WithFields(logrus.Fields{
						"isMappableFw":  isMappableFw,
						"isMappableRv":  isMappableRv,
						"num resultsFw": len(alternativeMappingsFw),
						"num resultsRv": len(alternativeMappingsRv),
					}).Debug("readpair not mappable")
					return
				}
			}
		}
	}
	wg.Done()

	logrus.Debug("Finished paralog mapping")
}

// Returns a slice of parentSeqIndices for a given readPairResult. If we map only to one main target region,
// there is only one paralog index (and we do not need to do anything), but if the index holds more than one seq, we do not want to map to all possible
// indices, only to the paralogs of the parentSeqIndex
func getParentIndex(fwMatches []*mapperutils.ReadMatchResult, rvMatches []*mapperutils.ReadMatchResult) map[int]struct{} {
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
