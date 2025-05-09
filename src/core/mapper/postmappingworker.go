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

	results := make(map[int][]*mappedreadpair.ReadPairMatchResult)
	intervalsPerSeq := make(map[int][][2]int)
	for mappedReadPair := range mappedReadPairChan {
		results[mappedReadPair.Fw.SequenceIndex/2] = append(results[mappedReadPair.Fw.SequenceIndex/2], mappedReadPair)
		for _, region := range mappedReadPair.Fw.MatchedGenome.Regions {
			intervalsPerSeq[mappedReadPair.Fw.SequenceIndex/2] = append(intervalsPerSeq[mappedReadPair.Fw.SequenceIndex/2], [2]int{region.Start, region.End})
		}
		for _, region := range mappedReadPair.Rv.MatchedGenome.Regions {
			convStart := 844 - region.End
			convEnd := 844 - region.Start
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
