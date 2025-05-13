package mapper

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core/algorithms"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/graph"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mappedreadpair"
	"github.com/sirupsen/logrus"
	"sync"
)

func PostMappingWorker(mappedReadPairChan <-chan *mappedreadpair.ReadPairMatchResult, wg *sync.WaitGroup, outputChan chan<- string) {

	logrus.Info("Started accumulating mapped readpairs")

	// in intervalsPerSeq, the key 0 references GenomeIndex.Sequences[0]
	// that's why we do intervalsPerSeq[mappedReadPair.Fw.SequenceIndex/2] and intervalsPerSeq[mappedReadPair.Rv.SequenceIndex/2]
	// it projects 0,1 -> 0; 2,3 -> 1; 4,5 -> 2
	resultsPerSeqIndex := make(map[int][]*mappedreadpair.ReadPairMatchResult)
	intervalsPerSeq := make(map[int][][2]int)
	for mappedReadPair := range mappedReadPairChan {
		resultsPerSeqIndex[mappedReadPair.Fw.SequenceIndex/2] = append(resultsPerSeqIndex[mappedReadPair.Fw.SequenceIndex/2], mappedReadPair)
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
	logrus.Info("Finished accumulating mapped readpairs")

	for seqId, intervals := range intervalsPerSeq {
		fmt.Println(seqId)
		res := algorithms.FindHighCoverageRegions(intervals)
		fmt.Println(res)
	}

	logrus.Debug("Identifying critical intervalsPerSeq")
	fmt.Println(len(resultsPerSeqIndex))

	for _, mappedReadPairs := range resultsPerSeqIndex {
		og := graph.OverlapGraph{
			Nodes:     map[int]*graph.Node{},
			Index:     mappedReadPairs[0].Index,
			NodeCount: 0,
		}
		for _, mappedReadPair := range mappedReadPairs {
			og.InsertNode(mappedReadPair)
		}
		fmt.Println(og.ToString())
		fmt.Println()
	}

	// WRITE TO SAM
	for _, mappedReadPairs := range resultsPerSeqIndex {
		for _, mappedReadPair := range mappedReadPairs {
			fwRes := mappedReadPair.Fw
			rvRes := mappedReadPair.Rv
			output, isMapping := FormatMappedReadPairToSAM(*fwRes, *rvRes, mappedReadPair.ReadPair, mappedReadPair.Index)
			if isMapping {
				outputChan <- output
			}
		}
	}

	wg.Done()
}
