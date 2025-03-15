package mapper

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/matchutils"
	"github.com/KleinSamuel/gtamap/src/core/timer"
	"github.com/KleinSamuel/gtamap/src/datawriter"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
	"github.com/sirupsen/logrus"
	"os"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

type MappingTask struct {
	ID       int
	ReadPair *fastq.ReadPair
}

func MapperWorker(workerId int, taskQueue <-chan MappingTask, taskQueueWriter chan<- string,
	wg *sync.WaitGroup, genomeIndex *index.GenomeIndex, timerChannel chan<- *timer.Timer) {

	logrus.WithFields(logrus.Fields{
		"workerId": workerId,
	}).Debug("Started MapperWorker")

	for task := range taskQueue {

		logrus.WithFields(logrus.Fields{
			"workerId": workerId,
			"task":     task.ID,
		}).Debug("Processing task")

		result, isMappable := MapReadPair(task.ReadPair, genomeIndex, timerChannel)

		if isMappable {
			taskQueueWriter <- result
		}
	}

	// signal that the worker has finished
	wg.Done()

	logrus.WithFields(logrus.Fields{
		"workerId": workerId,
	}).Debug("Finished MapperWorker")
}

func OutputWorker(taskQueue <-chan string, wg *sync.WaitGroup, writer *datawriter.Writer) {

	logrus.Debug("Started writeOutputWorker")

	for task := range taskQueue {
		writer.Write(task)
	}

	wg.Done()

	logrus.Debug("Finished writeOutputWorker")
}

func TimerWorker(timerChannel <-chan *timer.Timer, wg *sync.WaitGroup) {
	defer wg.Done()

	totalMapR1 := time.Duration(0)
	totalMapR2 := time.Duration(0)
	totalDetermineReadLocation := time.Duration(0)
	totalExactMatch := time.Duration(0)

	for t := range timerChannel {
		totalMapR1 += t.MapR1
		totalMapR2 += t.MapR2
		totalDetermineReadLocation += t.DetermineReadLocation
		totalExactMatch += t.ExactMatch
	}

	logrus.Info("Total mapping time: ", totalMapR1+totalMapR2)
	logrus.Info("Total mapping time R1: ", totalMapR1)
	logrus.Info("Total mapping time R2: ", totalMapR2)
	logrus.Info("Total determine read location time: ", totalDetermineReadLocation)
	logrus.Info("Total exact match time: ", totalExactMatch)
}

func MapAll(genomeIndex *index.GenomeIndex, reader *fastq.Reader, writer *datawriter.Writer) {

	runtime.GOMAXPROCS(runtime.NumCPU())

	timerStartTotal := time.Now()

	samHeader := sam.Header{
		Version:                 sam.Version,
		ReferenceSequenceName:   "test",
		ReferenceSequenceLength: int(len(*genomeIndex.Sequence)),
		GenomeAnnotationVersion: "",
		GenomeAssemblyVersion:   "",
		OrganismTaxId:           "",
		ToolVersion:             config.ToolVersion(),
	}

	writer.Write(samHeader.String())

	numWorkers := 1

	taskQueueMapping := make(chan MappingTask)
	taskQueueWriter := make(chan string)
	timerChannel := make(chan *timer.Timer)

	// wait group that keeps track of the mapping goroutines that are still running
	var waitgroupMapping sync.WaitGroup

	// start the mapping worker goroutine pool
	for i := 0; i < numWorkers; i++ {
		waitgroupMapping.Add(1)
		go MapperWorker(i, taskQueueMapping, taskQueueWriter, &waitgroupMapping, genomeIndex, timerChannel)
	}

	// wait group that keeps track of the writer goroutine that is still running
	var waitgroupWriter sync.WaitGroup
	// start the writer goroutine
	waitgroupWriter.Add(1)
	go OutputWorker(taskQueueWriter, &waitgroupWriter, writer)

	var waitGroupTimer sync.WaitGroup
	waitGroupTimer.Add(1)
	go TimerWorker(timerChannel, &waitGroupTimer)

	// the number of read pairs that have been processed
	taskCounter := 0

	// add each read pair as a mapping task to the task queue
	for readPair := reader.NextRead(); readPair != nil; readPair = reader.NextRead() {

		// TODO: remove after testing (only process specific read pair)
		//name := strings.Split(readPair.ReadR1.Header, " ")[0]
		//if name != "@3-0000/7-fw" {
		//	continue
		//}
		//if name != "@90" {
		//	continue
		//}

		mappingTask := MappingTask{
			ID:       taskCounter,
			ReadPair: readPair,
		}
		taskQueueMapping <- mappingTask

		taskCounter++

		// TODO: remove after testing
		if taskCounter == 400 {
			break
		}
	}

	close(taskQueueMapping)

	waitgroupMapping.Wait()

	close(timerChannel)
	waitGroupTimer.Wait()

	close(taskQueueWriter)

	waitgroupWriter.Wait()
	writer.Close()

	fmt.Println("mapping finished")
	fmt.Println("num tasks: ", taskCounter)

	totalDuration := time.Since(timerStartTotal)

	logrus.WithFields(logrus.Fields{
		"duration": totalDuration,
		"io":       reader.Duration,
	}).Info("Finished mapping")

	fmt.Println("total mapping time: ", totalDuration)
}

func MapReadPair(readPair *fastq.ReadPair, genomeIndex *index.GenomeIndex, timerChannel chan<- *timer.Timer) (string, bool) {

	keepFw := Filter(readPair.ReadR1.Sequence, genomeIndex)
	keepRw := Filter(readPair.ReadR2.Sequence, genomeIndex)

	if !keepFw || !keepRw {
		return "", false
	}

	resultFw, isMappableFw := MapRead(readPair.ReadR1, genomeIndex)
	resultRw, isMappableRw := MapRead(readPair.ReadR2, genomeIndex)

	if !isMappableFw || !isMappableRw {
		return "", false
	}

	if len(resultFw) == 0 || len(resultRw) == 0 {
		return "", false
	}

	if len(resultFw) > 1 || len(resultRw) > 1 {
		// TODO: handle multimapping reads
		return readPair.ReadR1.Header + "_fw\tMULTIMAPPING\t" + strconv.Itoa(len(resultFw)) + "\n" +
			readPair.ReadR1.Header + "_rw\tMULTIMAPPING\t" + strconv.Itoa(len(resultRw)), true
	}

	resFw := resultFw[0]
	resRw := resultRw[0]

	var builder strings.Builder

	// QNAME
	builder.WriteString(readPair.ReadR1.Header)
	builder.WriteString("\t")
	// FLAG
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// RNAME
	builder.WriteString(strconv.Itoa(resFw.SequenceIndex))
	builder.WriteString("\t")
	// POS
	startRelativeFw := resFw.MatchedGenome.GetFirstRegion().Start
	startGenomeFw := 45884425 + startRelativeFw

	builder.WriteString(strconv.Itoa(startGenomeFw))
	builder.WriteString("\t")
	// MAPQ
	builder.WriteString(strconv.Itoa(255))
	builder.WriteString("\t")
	// CIGAR
	builder.WriteString(resFw.GetCigar())
	builder.WriteString("\t")
	// PNEXT
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// RNEXT
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// TLEN
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// SEQ
	builder.WriteString(string(*readPair.ReadR1.Sequence))
	builder.WriteString("\t")
	// QUAL
	builder.WriteString(string(*readPair.ReadR1.Quality))
	builder.WriteString("\n")

	// R2 READ

	// QNAME
	builder.WriteString(readPair.ReadR2.Header)
	builder.WriteString("\t")
	// FLAG
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// RNAME
	builder.WriteString(strconv.Itoa(resRw.SequenceIndex))
	builder.WriteString("\t")
	// POS
	startRelativeRw := resRw.MatchedGenome.GetLastRegion().End
	startGenomeRw := 45903174 - startRelativeRw + 1
	builder.WriteString(strconv.Itoa(startGenomeRw))
	builder.WriteString("\t")
	// MAPQ
	builder.WriteString(strconv.Itoa(255))
	builder.WriteString("\t")
	// CIGAR
	builder.WriteString(resRw.GetCigar())
	builder.WriteString("\t")
	// PNEXT
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// RNEXT
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// TLEN
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// SEQ
	builder.WriteString(string(*readPair.ReadR2.Sequence))
	builder.WriteString("\t")
	// QUAL
	builder.WriteString(string(*readPair.ReadR2.Quality))
	builder.WriteString("\n")

	//resString := readPair.ReadR1.Header + "_fw\t"
	//for _, resFw := range resultFw {
	//	startRelative := resFw.MatchedGenome.GetFirstRegion().Start
	//	startGenome := 45884425 + startRelative
	//	resString += strconv.Itoa(resFw.SequenceIndex) + ":" + strconv.Itoa(startRelative) + "|" + strconv.Itoa(startGenome) + ":" + strconv.Itoa(resFw.MatchedRead.Length()) + " "
	//}
	//resString = resString + "\n"
	//
	//resString += readPair.ReadR2.Header + "_rw\t"
	//for _, resRw := range resultRw {
	//	startRelative := resRw.MatchedGenome.GetLastRegion().End
	//	startGenome := 45903174 - startRelative + 1
	//	resString += strconv.Itoa(resRw.SequenceIndex) + ":" + strconv.Itoa(startRelative) + "|" + strconv.Itoa(startGenome) + ":" + strconv.Itoa(resRw.MatchedRead.Length()) + " "
	//}
	//resString = resString + "\n"

	return builder.String(), true
}

// Filter the read sequence by checking if it contains at least 5 matching kmers
func Filter(readSequence *[]byte, genomeIndex *index.GenomeIndex) bool {

	numMatching := 0

	numMatchingFw := 0
	numMatchingRw := 0

	for i := 0; i <= len(*readSequence)-(int(config.KmerLength())); i += int(config.KmerLength()) {

		kmer := (*readSequence)[i : i+int(config.KmerLength())]

		//matches := genomeIndex.KeywordTree.FindKeyword(&kmer, i)
		matches := genomeIndex.GetKeywordFromMap(*(*[10]byte)(kmer))

		okFw := false
		okRv := false
		for _, m := range matches {
			if m.SequenceIndex == 0 {
				okFw = true
			} else {
				okRv = true
			}
		}

		if okFw {
			numMatchingFw++
		}
		if okRv {
			numMatchingRw++
		}

		if matches != nil {
			numMatching++
		}
	}

	//fmt.Println("numMatchingFw: ", numMatchingFw)
	//fmt.Println("numMatchingRw: ", numMatchingRw)

	return numMatching >= 7
}

func sortedIndicesDesc(list []int) []int {
	indices := make([]int, len(list))
	for i := range list {
		indices[i] = i
	}
	sort.Slice(indices, func(i, j int) bool {
		return list[indices[i]] > list[indices[j]]
	})
	return indices
}

func MapRead(read *fastq.Read, genomeIndex *index.GenomeIndex) ([]matchutils.ReadMatchResult, bool) {

	logrus.WithFields(logrus.Fields{
		"read": read.Header,
	}).Info("Mapping read")

	globalMatches := matchutils.GlobalMatchResult{
		MatchesPerSequence: make([]*matchutils.SequenceMatchResult, genomeIndex.KeywordTree.NumSequences),
	}

	// create all non-overlapping k-mers for the read pair
	for i := 0; i <= len(*read.Sequence)-(int(config.KmerLength())); i += int(config.KmerLength()) {

		kmer := (*read.Sequence)[i : i+int(config.KmerLength())]

		matches := genomeIndex.KeywordTree.FindKeyword(&kmer, i)

		if matches == nil {
			continue
		}

		for _, match := range matches {

			// add new sequence to global matchutils
			if globalMatches.MatchesPerSequence[match.SequenceIndex] == nil {
				globalMatches.MatchesPerSequence[match.SequenceIndex] = &matchutils.SequenceMatchResult{
					MatchesPerDiagonal: make(map[int][]*matchutils.Match),
				}
			}

			// add new diagonal to sequence match
			if globalMatches.MatchesPerSequence[match.SequenceIndex].MatchesPerDiagonal[match.StartGenome] == nil {
				globalMatches.MatchesPerSequence[match.SequenceIndex].MatchesPerDiagonal[match.StartGenome] =
					make([]*matchutils.Match, 0)
			}

			// add match to diagonal
			globalMatches.MatchesPerSequence[match.SequenceIndex].MatchesPerDiagonal[match.StartGenome] =
				append(globalMatches.MatchesPerSequence[match.SequenceIndex].MatchesPerDiagonal[match.StartGenome], match)
		}
	}

	// list has size of number of sequences in index
	// each element represents the maximum number of kmers that matched exactly on the same diagonal
	maxDiagonalHitsPerSequence := make([]int, genomeIndex.KeywordTree.NumSequences)

	for seqIndex, sequenceMatches := range globalMatches.MatchesPerSequence {

		maxHits := 0

		if sequenceMatches != nil {
			for _, matches := range sequenceMatches.MatchesPerDiagonal {
				if len(matches) > maxHits {
					maxHits = len(matches)
				}
			}
		}

		maxDiagonalHitsPerSequence[seqIndex] = maxHits
	}

	seqIndexSorted := sortedIndicesDesc(maxDiagonalHitsPerSequence)

	// minimum number of hits any diagonal must have to be considered
	MIN_HITS := 5

	for i, seqIndex := range seqIndexSorted {
		if maxDiagonalHitsPerSequence[seqIndex] < MIN_HITS {
			seqIndexSorted[i] = -1
		}
	}

	logrus.WithFields(logrus.Fields{
		"maxDiagonalHitsPerSequence": maxDiagonalHitsPerSequence,
		"seqIndexSorted":             seqIndexSorted,
	}).Info("Potential sequences")

	//for seqIndex, sequenceMatches := range globalMatches.MatchesPerSequence {
	//
	//	if sequenceMatches == nil {
	//		continue
	//	}
	//
	//	dSizes := make([]int, len(sequenceMatches.MatchesPerDiagonal))
	//	dIndices := make([]int, len(sequenceMatches.MatchesPerDiagonal))
	//	counter := 0
	//
	//	for diagonal, matches := range sequenceMatches.MatchesPerDiagonal {
	//		fmt.Println(seqIndex, diagonal, len(matches))
	//
	//		dSizes[counter] = len(matches)
	//		dIndices[counter] = diagonal
	//		counter++
	//	}
	//
	//	fmt.Println(dSizes)
	//	fmt.Println(dIndices)
	//}

	results := make([]matchutils.ReadMatchResult, 0)

sequenceLoop:
	for _, seqIndex := range seqIndexSorted {
		//for seqIndex, sequenceMatches := range globalMatches.MatchesPerSequence {

		if seqIndex == -1 {
			// -1 means this sequence was disregarded because of some constraints
			continue
		}

		sequenceMatches := globalMatches.MatchesPerSequence[seqIndex]

		if sequenceMatches == nil {
			continue
		}

		logrus.Info("computing best match for sequence with index: ", seqIndex)

		// 1. find the best diagonal (most matches)
		// 2. find gaps in diagonal (unmatched regions via regionvector)
		// 3. fill gaps in digonal (each position for a match or mismatch)
		// 4. iterate repeat 1-4 until there are no fitting diagonals left
		// now there are multiple options left:
		// - the read is fully matched already -> done
		// - there are multiple diagonals with unmatched positions in between
		//   (probably due to a junction -> extend both diagonals as far as possible)
		// - there are unmatched positions in the front or back of the read
		//   (extend but if not possible then search for another diagonal with shorter kmer)

		result := matchutils.ReadMatchResult{
			SequenceIndex:  seqIndex,
			MatchedRead:    regionvector.NewRegionVector(),
			MatchedGenome:  regionvector.NewRegionVector(),
			MismatchesRead: make([]int, 0),
		}

		diagonalHandler := matchutils.NewDiagonalHandlerWithData(sequenceMatches.MatchesPerDiagonal)

		for result.MatchedRead.Length() < len(*read.Sequence) {

			logrus.Info("searching for next best diagonal")

			bestDiagonal, bestDiagonalLength := diagonalHandler.GetBestDiagonal()

			if bestDiagonal == -1 {
				logrus.Info("no suitable diagonal found")
				break
			}

			diagonalRead := regionvector.NewRegionVector()
			diagonalGenome := regionvector.NewRegionVector()

			for _, match := range sequenceMatches.MatchesPerDiagonal[bestDiagonal] {
				// the region can not be part of another diagonal that is already used
				if match.Used {
					continue
				}
				diagonalRead.AddRegion(match.FromRead, match.ToRead)
				diagonalGenome.AddRegion(match.FromGenome, match.ToGenome)
			}

			logrus.WithFields(logrus.Fields{
				"posGenome": bestDiagonal,
				"length":    bestDiagonalLength,
				"read":      diagonalRead,
				"genome":    diagonalGenome,
				"gaps":      len(diagonalRead.Regions) - 1,
			}).Info("found best diagonal")

			// find gaps in diagonal and fill them (matches and mismatches)
			mismatches := make([]int, 0)

			foundGap := len(diagonalRead.Regions) > 1

			// more than one element means there is a gap
			// the regionvectors of read and genome should have the same length as they are coupled
			for len(diagonalRead.Regions) > 1 {

				gapRead := diagonalRead.GetFirstGap()
				gapGenome := diagonalGenome.GetFirstGap()

				logrus.WithFields(logrus.Fields{
					"read":   gapRead,
					"genome": gapGenome,
				}).Info("found gap")

				diagonalRead.AddRegion(gapRead.Start, gapRead.End)
				for i := gapRead.Start; i < gapRead.End; i++ {

					readByte := (*read.Sequence)[i]
					gIndex := diagonalGenome.GetFirstRegion().Start + i
					genomeByte := (*genomeIndex.Sequence)[gIndex]

					if readByte != genomeByte {
						mismatches = append(mismatches, i)

						logrus.Info("mismatch at read position: ", i)

						// TODO: make max allowed mismatches configurable
						if len(mismatches) > 5 {
							logrus.Info("too many mismatches")
							continue sequenceLoop
						}
					}
				}
				diagonalGenome.AddRegion(gapGenome.Start, gapGenome.End)
			}

			if foundGap {
				logrus.WithFields(logrus.Fields{
					"read":       diagonalRead,
					"genome":     diagonalGenome,
					"mismatches": mismatches,
				}).Info("filled gap")
			}

			// add digonal to result (diagonal rv should have length 1 after gap filling)
			result.MatchedRead.AddRegion(diagonalRead.GetFirstRegion().Start, diagonalRead.GetFirstRegion().End)
			result.MatchedGenome.AddRegion(diagonalGenome.GetFirstRegion().Start, diagonalGenome.GetFirstRegion().End)
			result.MismatchesRead = append(result.MismatchesRead, mismatches...)

			// remove diagonal from handler
			diagonalHandler.ConsumeDiagonal(bestDiagonal)
		}

		logrus.WithFields(logrus.Fields{
			"read":   result.MatchedRead,
			"genome": result.MatchedGenome,
		}).Info("done processing diagonals")

		if result.MatchedRead.Length() == len(*read.Sequence) {
			logrus.Info("MATCHED FULLY")
		} else {

			if len(result.MatchedRead.Regions) > 1 {

				logrus.Info("UNMATCHED POSITIONS WITHIN READ")

				// resolve the borders within the read where at least one junction was found
				for len(result.MatchedRead.Regions) > 1 {

					// position in read where the gap starts
					lRead := result.MatchedRead.GetFirstGap().Start
					// position in genome where the gap starts
					lGenome := result.MatchedGenome.GetFirstGap().Start
					// position in read where the gap ends
					rRead := result.MatchedRead.GetFirstGap().End
					// position in genome where the gap ends
					rGenome := result.MatchedGenome.GetFirstGap().End

					//readSequence := (*read.Sequence)[lRead:rRead]
					//genomeSequenceL := (*genomeIndex.Sequence)[lGenome : lGenome+15]
					//genomeSequenceR := (*genomeIndex.Sequence)[rGenome-15 : rGenome]
					//fmt.Println("READ:\t\t", string(readSequence))
					//fmt.Println("GENOME L:\t", string(genomeSequenceL))
					//fmt.Println("GENOME R:\t", string(genomeSequenceR))
					//fmt.Println(lRead, lGenome, rRead, rGenome)

					// the length of the extension based on the number if unmapped bases in the read
					extensionLength := rRead - lRead

					// cululative mismatch count for the left and right extensions
					// for lErrors the index i represents the number of mismatches for the first i positions of the extension
					// for rErrors the index i represents the number of mismatches for the last i positions of the extension
					lErrors := make([]int, extensionLength+1)
					rErrors := make([]int, extensionLength+1)

					lErrors[0] = 0
					rErrors[0] = 0

					for i := 1; i <= extensionLength; i++ {
						lErrors[i] = lErrors[i-1]
						if (*read.Sequence)[lRead+i-1] != (*genomeIndex.Sequence)[lGenome+i-1] {
							lErrors[i]++
						}

						rErrors[i] = rErrors[i-1]
						if (*read.Sequence)[rRead-i] != (*genomeIndex.Sequence)[rGenome-i] {
							rErrors[i]++
						}
					}

					// the minimum number of mismatches
					minErrors := lErrors[extensionLength] + rErrors[extensionLength]
					// the position of the split with the minimum number of mismatches
					minSplit := -1

					// TODO: keep track of the actual mismatch positions
					// TODO: if no suitable split is found then:
					// - maybe there is another exon in between if enough bases missing from read
					// - maybe keep the readpair for second pass
					for i := 0; i <= extensionLength; i++ {

						lPos := i
						rPos := extensionLength - i

						numMismatches := lErrors[lPos] + rErrors[rPos]

						donorSiteStart := result.MatchedGenome.GetFirstGap().Start + i
						donorSiteSeq := (*genomeIndex.Sequence)[donorSiteStart : donorSiteStart+2]

						splitRev := extensionLength - i
						acceptorSiteStart := result.MatchedGenome.GetFirstGap().End - splitRev
						acceptorSiteSeq := (*genomeIndex.Sequence)[acceptorSiteStart-2 : acceptorSiteStart]

						if donorSiteSeq[0] == byte('G') && donorSiteSeq[1] == byte('T') {
							// canonical splice donor site dinucleotide GT
							numMismatches += 0
						} else if donorSiteSeq[0] == byte('G') && donorSiteSeq[1] == byte('C') {
							// non-canonical splice donor site dinucleotide GC
							numMismatches += 1
						} else if donorSiteSeq[0] == byte('A') && donorSiteSeq[1] == byte('T') {
							// non-canonical splice donor site dinucleotide GC
							numMismatches += 1
						} else {
							// all other splice donor site dinucleotides
							numMismatches += 2
						}

						if acceptorSiteSeq[0] == byte('A') && acceptorSiteSeq[1] == byte('G') {
							// canonical splice acceptor site dinucleotide AG
							numMismatches += 0
						} else if acceptorSiteSeq[0] == byte('A') && acceptorSiteSeq[1] == byte('C') {
							numMismatches += 1
						} else {
							numMismatches += 2
						}

						if numMismatches < minErrors {
							minErrors = numMismatches
							minSplit = i
						}
					}

					// apply the best split
					if minSplit == -1 {
						logrus.Error("minSplit is zero")
						os.Exit(1)
					}

					result.MatchedRead.AddRegion(lRead, lRead+minSplit)
					result.MatchedGenome.AddRegion(lGenome, lGenome+minSplit)

					result.MatchedRead.AddRegion(rRead-(extensionLength-minSplit), rRead)
					result.MatchedGenome.AddRegion(rGenome-(extensionLength-minSplit), rGenome)
				}
			}

			if result.MatchedRead.GetFirstRegion().Start > 0 {
				logrus.Info("UNMATCHED POSITIONS IN FRONT OF READ")

				//os.Exit(1)

				//// read sequence
				//unmappedRead := (*read.Sequence)[0:result.MatchedRead.GetFirstRegion().Start]
				//// genome sequence
				//endGenome := result.MatchedGenome.GetFirstRegion().Start
				//startGenome := result.MatchedGenome.GetFirstRegion().Start - result.MatchedRead.GetFirstRegion().Start
				//unmappedGenome := (*genomeIndex.Sequence)[startGenome:endGenome]
				//
				//mismatches := make([]int, 0)
				//
				//for i := 0; i < len(unmappedRead); i++ {
				//	if unmappedRead[i] != unmappedGenome[i] {
				//		mismatches = append(mismatches, i)
				//	}
				//}
				//
				//result.MatchedRead.AddRegion(0, result.MatchedRead.GetFirstRegion().Start)
				//result.MatchedGenome.AddRegion(startGenome, endGenome)
				//result.MismatchesRead = append(result.MismatchesRead, mismatches...)
			}

			if result.MatchedRead.GetLastRegion().End < len(*read.Sequence) {
				logrus.Info("UNMATCHED POSITIONS IN BACK OF READ")

				//// read sequence
				//startRead := result.MatchedRead.GetLastRegion().End
				//endRead := len(*read.Sequence)
				//unmappedRead := (*read.Sequence)[startRead:endRead]
				//// genome sequence
				//startGenome := result.MatchedGenome.GetLastRegion().End
				//endGenome := result.MatchedGenome.GetLastRegion().End + (endRead - startRead)
				//unmappedGenome := (*genomeIndex.Sequence)[startGenome:endGenome]
				//
				//mismatches := make([]int, 0)
				//
				//for i := 0; i < len(unmappedRead); i++ {
				//	if unmappedRead[i] != unmappedGenome[i] {
				//		mismatches = append(mismatches, i)
				//	}
				//}
				//
				//result.MatchedRead.AddRegion(startRead, endRead)
				//result.MatchedGenome.AddRegion(startGenome, endGenome)
				//result.MismatchesRead = append(result.MismatchesRead, mismatches...)
			}

		}

		results = append(results, result)

		//break
	}

	if len(results) == 0 {
		return results, false
	}
	return results, true
}
