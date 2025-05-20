package mapper

import (
	"fmt"
	"strconv"
	"strings"
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/sirupsen/logrus"
)

func PostMappingWorker(mappedReadPairChan <-chan *ReadPairMatchResults, wg *sync.WaitGroup, outputChan chan<- string) {
	logrus.Info("Started accumulating mapped readpairs")
	defer wg.Done()

	// I. Extract and sort mappring readpairs into map: mainSeqIndex -> [mappings]
	resultsPerSeqIndex := make(map[int][]*ReadPairMatchResults)
	for mappedReadPair := range mappedReadPairChan {
		fmt.Println(mappedReadPair)

		// TODO: Handle multimappings
		if len(mappedReadPair.Fw) > 1 || len(mappedReadPair.Rv) > 1 {

			if len(mappedReadPair.Fw) != len(mappedReadPair.Rv) {
			}

			logrus.Fatalf("Multimapping readpair found. Currently not handled: Read ID: %s", mappedReadPair.ReadPair.ReadR1.Header)
		}

		if mappedReadPair.Fw[0].SequenceIndex/2 != mappedReadPair.Rv[0].SequenceIndex/2 {
			logrus.Infof("Mates of Readpair don't agree on seq index %s", mappedReadPair.ReadPair.ReadR1.Header)
		}

		// For now I simply take the first Fw and Rv
		// projects 0,1 -> 0; 2,3 -> 1; 4,5 -> 2 to map to main taret seq id
		parentSeqIndex := mappedReadPair.Fw[0].SequenceIndex / 2

		resultsPerSeqIndex[parentSeqIndex] = append(resultsPerSeqIndex[parentSeqIndex], mappedReadPair)
	}
	logrus.Info("Finished accumulating mapped readpairs")

	// II. extract coverage
	// Here we assume that there are already only one fw and rv per mappedREadPair
	// I will have to create a new struct which only holds one fw/rv but other than that is the exact same
	// as ReadPairMatchResults
	multiMapsPerSeqIndex := make(map[string][]*ReadPairMatchResults)
	for seqIndex, mappedReadPairs := range resultsPerSeqIndex {
		coverageSliceFw, coverageSliceRv := getCoverageSlices(mappedReadPairs, len((*mappedReadPairs[0].Index.Sequences[seqIndex])))
		highCovRegionsfw := getHighCoverageRegons(coverageSliceFw)
		highCovRegionsrv := getHighCoverageRegons(coverageSliceRv)
		multimaps := getMultimappingReads(mappedReadPairs, highCovRegionsfw, highCovRegionsrv)

		targetRegionId := mappedReadPairs[0].Index.GetSequenceInfo(seqIndex).GeneId
		multiMapsPerSeqIndex[targetRegionId] = multimaps
	}

	paralogReads := make(map[string]struct{})

	for seqIndex, ambiguousReadPairs := range multiMapsPerSeqIndex {
		for _, ambiguousReadPair := range ambiguousReadPairs {
			// currMM := len(ambiguousReadPair.Rv[0].MismatchesRead) + len(ambiguousReadPair.Fw[0].MismatchesRead)

			if !(ambiguousReadPair.ReadPair.ReadR1.Header == "A00604:202:HLYW3DSXY:3:2670:22363:10880 2:N:0:AACTCGGA+TCTGGACA") {
				continue
			}

			// var lastFw []mapperutils.ReadMatchResult
			// var lastRv []mapperutils.ReadMatchResult

			for i := 0; i < 10; i++ {
				// alternativeMappingsFw, isMappableFw := MapRead(ambiguousReadPair.ReadPair.ReadR1, ambiguousReadPair.Index.ParalogRegions[seqIndex])
				alternativeMappingsRv, _ := MapRead(ambiguousReadPair.ReadPair.ReadR2, ambiguousReadPair.Index.ParalogRegions[seqIndex])

				fmt.Println(alternativeMappingsRv)
				// if i == 0 {
				// 	lastFw = alternativeMappingsFw
				// 	lastRv = alternativeMappingsRv
				// 	continue
				// }
				//
				// if len(lastFw) != len(alternativeMappingsFw) {
				// 	fmt.Println("Fw Map differs")
				// 	fmt.Println(lastFw)
				// 	fmt.Println(alternativeMappingsFw)
				// }
				// if len(lastRv) != len(alternativeMappingsRv) {
				// 	fmt.Println("Rv Map differs")
				// 	fmt.Println(lastRv)
				// 	fmt.Println(alternativeMappingsRv)
				// }
				// lastFw = alternativeMappingsFw
				// lastRv = alternativeMappingsRv

				// if !isMappableFw || !isMappableRv {
				// 	continue
				// }

			}
			fmt.Println("-------------------------------------")

			// here we make use of the fact that we are using paired end reads
			// this means we only consider a valid mapping as: fw and rv agree on main seq index and are on
			// diff strands
			// if ambiguousReadPair.ReadPair.ReadR2.Header == "A00604:202:HLYW3DSXY:3:1275:17381:36855 1:N:0:AACTCGGA+TCTGGACA" {
			// 	fmt.Println()
			// }
			// bestParalogId, _, minMM := getPossibleMapping(alternativeMappingsFw, alternativeMappingsRv, currMM)
			//
			// if bestParalogId == -1 {
			// 	logrus.Debugf("No alternative mapping found for %s with %s mismatches", ambiguousReadPair.ReadPair.ReadR1.Header, strconv.Itoa(currMM))
			// } else {
			// 	logrus.Debugf("Alternative mapping found for %s in paralog region %s with %s less mismatches", ambiguousReadPair.ReadPair.ReadR1.Header, bestParalogId, strconv.Itoa(currMM-minMM))
			// 	paralogReads[ambiguousReadPair.ReadPair.ReadR2.Header] = struct{}{}
			// 	i++
			//
			// 	// fmt.Println(possibleMapping)
			// 	// fmt.Println(possibleMapping[0].SequenceIndex)
			// 	// fmt.Println(possibleMapping[1].SequenceIndex)
			// }
		}
		fmt.Println(len(ambiguousReadPairs))
	}
	// fmt.Println(i)

	// WRITE TO SAM
	for _, mappedReadPairs := range resultsPerSeqIndex {
		for _, mappedReadPair := range mappedReadPairs {
			resultFw := mappedReadPair.Fw
			resultRv := mappedReadPair.Rv

			if len(resultFw) > 1 || len(resultRv) > 1 {
				// TODO: handle multimapping reads

				//logrus.WithFields(logrus.Fields{
				//	"numResultsFw": len(resultFw),
				//	"numResultsRv": len(resultRv),
				//	"readFw":       readPair.ReadR1.Header,
				//	"readRv":       readPair.ReadR2.Header,
				//}).Warn("multimapping reads not handled yet")

				var builder strings.Builder

				// TODO: handle this appropriately
				// combine every first of pair read with every second of pair read

				for i := 0; i < len(resultFw); i++ {
					for j := 0; j < len(resultRv); j++ {
						s, isOk := readPairResultToSamString(mappedReadPair.Index, mappedReadPair.ReadPair, resultFw[i], resultRv[j])
						if !isOk {
							continue
						}
						builder.WriteString(s)
					}
				}

				//for _, resFw := range resultFw {
				//	s, isOk := readPairResultToSamString(genomeIndex, readPair, resFw, nil)
				//	if !isOk {
				//		continue
				//	}
				//	builder.WriteString(s)
				//}
				//
				//for _, resRv := range resultRv {
				//	s, isOk := readPairResultToSamString(genomeIndex, readPair, nil, resRv)
				//	if !isOk {
				//		continue
				//	}
				//	builder.WriteString(s)
				//}

				outputChan <- builder.String()
			}
			_, exists := paralogReads[mappedReadPair.ReadPair.ReadR2.Header]
			if exists {
				out, _ := readPairResultToSamString(mappedReadPair.Index, mappedReadPair.ReadPair, resultFw[0], resultRv[0])
				outputChan <- out
			}
		}
	}
}

// holds all relevant mapping information of potentially several hits per readpair
// bundled into one type
type ReadPairMatchResults struct {
	ReadPair *fastq.ReadPair
	Fw       []*mapperutils.ReadMatchResult
	Rv       []*mapperutils.ReadMatchResult
	Index    *index.GenomeIndex
}

func (i ReadPairMatchResults) String() string {
	var builder strings.Builder
	builder.Write([]byte("ReadPairR1 Header: "))
	builder.Write([]byte(i.ReadPair.ReadR1.Header))
	builder.Write([]byte("\n"))
	builder.Write([]byte("  <== FW MAPPINGS ==>"))
	builder.Write([]byte("\n"))
	for _, mapping := range i.Fw {
		builder.Write([]byte("\t SeqIndex: "))
		seqIndex := strconv.Itoa(mapping.SequenceIndex)
		builder.WriteString(seqIndex)
		builder.WriteString("\n")
		builder.WriteString("\t GENOME -> ")
		builder.WriteString(mapping.MatchedGenome.String())
		builder.WriteString("\n")
		builder.WriteString("\t READ   -> ")
		builder.WriteString(mapping.MatchedRead.String())
		builder.WriteString("\n")
		builder.WriteString("\t MISMAT -> ")
		ints := mapping.MismatchesRead
		strs := make([]string, len(ints))
		for i, v := range ints {
			strs[i] = strconv.Itoa(v)
		}
		builder.WriteString(strings.Join(strs, ","))
		builder.WriteString("\n")
	}
	builder.Write([]byte("  <== RV MAPPINGS ==>"))
	builder.Write([]byte("\n"))
	for _, mapping := range i.Rv {
		builder.Write([]byte("\t SeqIndex: "))
		seqIndex := strconv.Itoa(mapping.SequenceIndex)
		builder.WriteString(seqIndex)
		builder.WriteString("\n")
		builder.WriteString("\t GENOME -> ")
		builder.WriteString(mapping.MatchedGenome.String())
		builder.WriteString("\n")
		builder.WriteString("\t READ   -> ")
		builder.WriteString(mapping.MatchedRead.String())
		builder.WriteString("\n")
		builder.WriteString("\t MISMAT -> ")
		ints := mapping.MismatchesRead
		strs := make([]string, len(ints))
		for i, v := range ints {
			strs[i] = strconv.Itoa(v)
		}
		builder.WriteString(strings.Join(strs, ","))
		builder.WriteString("\n")
	}
	return builder.String()
}

func getCoverageSlice(mappedReadPairs []*ReadPairMatchResults, geneLength int) []int {
	coverage := make([]int, geneLength)
	for _, mappedReadPair := range mappedReadPairs {
		// forward reads
		for _, region := range mappedReadPair.Fw[0].MatchedGenome.Regions {
			if mappedReadPair.Index.IsSequenceForward(mappedReadPair.Fw[0].SequenceIndex) {
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
			if mappedReadPair.Index.IsSequenceForward(mappedReadPair.Rv[0].SequenceIndex) {
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

func getMultimappingReads(mappedReadPairs []*ReadPairMatchResults, fwIntervals [][2]int, rvIntervals [][2]int) []*ReadPairMatchResults {
	multimappings := make([]*ReadPairMatchResults, 0)
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

func getCoverageSlices(mappedReadPairs []*ReadPairMatchResults, geneLength int) ([]int, []int) {
	coverageFw := make([]int, geneLength)
	coverageRv := make([]int, geneLength)

	for _, mappedReadPair := range mappedReadPairs {
		// fw reads
		for _, region := range mappedReadPair.Fw[0].MatchedGenome.Regions {
			if mappedReadPair.Index.IsSequenceForward(mappedReadPair.Fw[0].SequenceIndex) {
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
			if mappedReadPair.Index.IsSequenceForward(mappedReadPair.Rv[0].SequenceIndex) {
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

func getPossibleMapping(alternativeMappingsFw []mapperutils.ReadMatchResult, alternativeMappingsRv []mapperutils.ReadMatchResult, mismatchesInTarget int) (int, []mapperutils.ReadMatchResult, int) {
	mapPerParalogSeqIndex := make(map[int][]mapperutils.ReadMatchResult)

	// accumulate mapping results per paralog main id
	for _, mapping := range alternativeMappingsFw {
		mapPerParalogSeqIndex[mapping.SequenceIndex/2] = append(mapPerParalogSeqIndex[mapping.SequenceIndex/2], mapping)
	}

	for _, mapping := range alternativeMappingsRv {
		mapPerParalogSeqIndex[mapping.SequenceIndex/2] = append(mapPerParalogSeqIndex[mapping.SequenceIndex/2], mapping)
	}

	bestId := -1
	minMis := 1000000
	for paralogSeqIndex, result := range mapPerParalogSeqIndex {
		// do fw and rv map?
		if len(result) != 2 {
			logrus.Debug("Mate unmapped or multimapped  in paralog mapping attempt.")
			continue
		}

		if len(result) > 2 {
			logrus.Fatalf("Read mapps to multiple locations on same region and strand: paralog id %s, result %s", paralogSeqIndex, result)
		}

		// only consider cases where fwId == rvId-1 || fwId-1 == rvId && fwId / 2 == rvId / 2
		if !(result[0].SequenceIndex == result[1].SequenceIndex-1 || result[1].SequenceIndex == result[0].SequenceIndex-1) {
			logrus.Debugf("Fw and Rv Reads map to same strand on same region: paralog id %s, result %s", paralogSeqIndex, result)
			continue
		}

		mismatches := len(result[0].MismatchesRead) + len(result[1].MismatchesRead)
		if bestId == -1 {
			bestId = paralogSeqIndex
			minMis = mismatches
			continue
		}

		if mismatches < minMis {
			bestId = paralogSeqIndex
			minMis = mismatches
		}

	}
	if minMis < mismatchesInTarget && bestId != 1 {
		return bestId, mapPerParalogSeqIndex[bestId], minMis
	}
	return -1, nil, -1
}
