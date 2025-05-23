package mapper

import (
	"strings"
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
)

func PostMappingWorker(finalMatchesChan <-chan *mapperutils.ReadPairMatchResults, wg *sync.WaitGroup, outputChan chan<- string, index *index.GenomeIndex) {
	defer wg.Done()

	// TODO Here comes the Ambiguity Logic

	// WRITE TO SAM
	for mappedReadPair := range finalMatchesChan {
		var builder strings.Builder

		for i := 0; i < len(mappedReadPair.Fw); i++ {
			for j := 0; j < len(mappedReadPair.Rv); j++ {
				s, isOk := readPairResultToSamString(index, mappedReadPair.ReadPair, mappedReadPair.Fw[i], mappedReadPair.Rv[j])
				if !isOk {
					continue
				}
				builder.WriteString(s)
			}
		}
		outputChan <- builder.String()
	}
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
