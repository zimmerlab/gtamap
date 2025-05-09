package mapper

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core/algorithms"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mappedreadpair"
	"github.com/sirupsen/logrus"
	"sync"
)

func PostMappingWorker(mappedReadPairChan <-chan *mappedreadpair.ReadPairMatchResult, wg *sync.WaitGroup, outputChan chan<- string) {

	logrus.Debug("Started accumulating mapped readpairs")

	// in intervalsPerSeq, the key 0 references GenomeIndex.Sequences[0]
	// that's why we do intervalsPerSeq[mappedReadPair.Fw.SequenceIndex/2] and intervalsPerSeq[mappedReadPair.Rv.SequenceIndex/2]
	// it projects 0,1 -> 0; 2,3 -> 1; 4,5 -> 2
	results := make(map[int][]*mappedreadpair.ReadPairMatchResult)
	intervalsPerSeq := make(map[int][][2]int)
	for mappedReadPair := range mappedReadPairChan {
		results[mappedReadPair.Fw.SequenceIndex/2] = append(results[mappedReadPair.Fw.SequenceIndex/2], mappedReadPair)
		for _, region := range mappedReadPair.Fw.MatchedGenome.Regions {
			intervalsPerSeq[mappedReadPair.Fw.SequenceIndex/2] = append(intervalsPerSeq[mappedReadPair.Fw.SequenceIndex/2], [2]int{region.Start, region.End})
		}

		// we need to convert rv coords to fw coords for correct coverage
		for _, region := range mappedReadPair.Rv.MatchedGenome.Regions {
			convStart := len(*mappedReadPair.Index.Sequences[mappedReadPair.Fw.SequenceIndex]) - region.End
			convEnd := len(*mappedReadPair.Index.Sequences[mappedReadPair.Rv.SequenceIndex]) - region.Start
			intervalsPerSeq[mappedReadPair.Rv.SequenceIndex/2] = append(intervalsPerSeq[mappedReadPair.Rv.SequenceIndex/2], [2]int{convStart, convEnd})
		}
	}

	for seqId, intervals := range intervalsPerSeq {
		fmt.Println(seqId)
		res := algorithms.FindHighCoverageRegions(intervals)
		fmt.Println(res)
	}

	logrus.Debug("Identifying critical intervalsPerSeq")
	fmt.Println(len(results))

	for _, mappedReadPairs := range results {
		for _, mappedReadPair := range mappedReadPairs {
			fwRes := mappedReadPair.Fw
			rvRes := mappedReadPair.Fw
			output, isMapping := FormatMappedReadPairToSAM(*fwRes, *rvRes, mappedReadPair.ReadPair, mappedReadPair.Index)
			if isMapping {
				outputChan <- output
			}
		}
	}

	wg.Done()
}
