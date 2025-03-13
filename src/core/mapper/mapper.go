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
		name := strings.Split(readPair.ReadR1.Header, " ")[0]
		//if name != "@3-0000/7-fw" {
		//	continue
		//}
		if name != "@7" {
			continue
		}

		mappingTask := MappingTask{
			ID:       taskCounter,
			ReadPair: readPair,
		}
		taskQueueMapping <- mappingTask

		taskCounter++

		// TODO: remove after testing
		if taskCounter == 1 {
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

	resString := readPair.ReadR1.Header + "_fw\t"
	for _, resFw := range resultFw {
		resString += strconv.Itoa(resFw.SequenceIndex) + ":" + strconv.Itoa(resFw.MatchedGenome.GetFirstRegion().Start) + " "
	}
	resString = resString + "\n"

	resString += readPair.ReadR2.Header + "_rw\t"
	for _, resRw := range resultRw {
		resString += strconv.Itoa(resRw.SequenceIndex) + ":" + strconv.Itoa(resRw.MatchedGenome.GetFirstRegion().Start) + " "
	}
	resString = resString + "\n"

	return resString, true
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

	fmt.Println("numMatchingFw: ", numMatchingFw)
	fmt.Println("numMatchingRw: ", numMatchingRw)

	return numMatching >= 7
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

		match := genomeIndex.KeywordTree.FindKeyword(&kmer, i)

		if match == nil {
			continue
		}

		for _, m := range match {

			// add new sequence to global matchutils
			if globalMatches.MatchesPerSequence[m.SequenceIndex] == nil {
				globalMatches.MatchesPerSequence[m.SequenceIndex] = &matchutils.SequenceMatchResult{
					MatchesPerDiagonal: make(map[int][]*matchutils.Match),
				}
			}

			// add new diagonal to sequence match
			if globalMatches.MatchesPerSequence[m.SequenceIndex].MatchesPerDiagonal[m.StartGenome] == nil {
				globalMatches.MatchesPerSequence[m.SequenceIndex].MatchesPerDiagonal[m.StartGenome] =
					make([]*matchutils.Match, 0)
			}

			// add match to diagonal
			globalMatches.MatchesPerSequence[m.SequenceIndex].MatchesPerDiagonal[m.StartGenome] =
				append(globalMatches.MatchesPerSequence[m.SequenceIndex].MatchesPerDiagonal[m.StartGenome], m)
		}
	}

	results := make([]matchutils.ReadMatchResult, 0)

	for seqIndex, sequenceMatches := range globalMatches.MatchesPerSequence {

		if sequenceMatches == nil {
			continue
		}

		for diagonal, matches := range sequenceMatches.MatchesPerDiagonal {
			fmt.Println(seqIndex, diagonal, len(matches))
		}
	}

sequenceLoop:
	for seqIndex, sequenceMatches := range globalMatches.MatchesPerSequence {

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
			fmt.Println("MATCHED FULLY")
		} else {

			if len(result.MatchedRead.Regions) > 1 {
				fmt.Println("UNMATCHED POSITIONS WITHIN READ")

				for len(result.MatchedRead.Regions) > 1 {

					fmt.Println("resolving borders because of junction")

					// position in read where the gap starts
					lRead := result.MatchedRead.GetFirstGap().Start
					// position in genome where the gap starts
					lGenome := result.MatchedGenome.GetFirstGap().Start
					// position in read where the gap ends
					rRead := result.MatchedRead.GetFirstGap().End
					// position in genome where the gap ends
					rGenome := result.MatchedGenome.GetFirstGap().End

					readSequence := (*read.Sequence)[lRead:rRead]
					genomeSequenceL := (*genomeIndex.Sequence)[lGenome : lGenome+15]
					genomeSequenceR := (*genomeIndex.Sequence)[rGenome-15 : rGenome]

					fmt.Println("READ:\t\t", string(readSequence))
					fmt.Println("GENOME L:\t", string(genomeSequenceL))
					fmt.Println("GENOME R:\t", string(genomeSequenceR))

					fmt.Println(lRead, lGenome, rRead, rGenome)

					extensionLength := rRead - lRead

					fmt.Println("extension length: ", extensionLength)

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

					fmt.Println("lErrors: ", lErrors)
					fmt.Println("rErrors: ", rErrors)

					splits := make([]int, 0, extensionLength+1)

					for i := 0; i <= extensionLength; i++ {

						lPos := i
						rPos := extensionLength - i

						numMismatches := lErrors[lPos] + rErrors[rPos]

						if numMismatches == 0 {
							splits = append(splits, i)
						}

						fmt.Println("lPos: ", lPos)
						fmt.Println("rPos: ", rPos)
						fmt.Println("numMismatches: ", numMismatches)
						fmt.Println("--")
					}

					fmt.Println("splits: ", splits)

					for _, split := range splits {
						fmt.Println("split: ", split)

						donorSiteStart := result.MatchedGenome.GetFirstGap().Start + split
						donorSiteSeq := (*genomeIndex.Sequence)[donorSiteStart : donorSiteStart+2]

						splitRev := extensionLength - split

						acceptorSiteStart := result.MatchedGenome.GetFirstGap().End - splitRev
						acceptorSiteSeq := (*genomeIndex.Sequence)[acceptorSiteStart-2 : acceptorSiteStart]

						fmt.Println("donor", string(donorSiteSeq))
						fmt.Println("acceptor", string(acceptorSiteSeq))
					}

					os.Exit(1)
				}

				//for len(result.MatchedRead.Regions) > 1 {
				//
				//	startRead := result.MatchedRead.GetFirstGap().Start
				//	endRead := result.MatchedRead.GetFirstGap().End
				//	unmappedRead := (*read.Sequence)[startRead:endRead]
				//
				//	startGenome := result.MatchedGenome.GetFirstGap().Start
				//	endGenome := result.MatchedGenome.GetFirstGap().End
				//	unmappedGenome := (*genomeIndex.Sequence)[startGenome:endGenome]
				//
				//	mismatches := make([]int, 0)
				//
				//	if len(unmappedRead) != len(unmappedGenome) {
				//		break
				//	}
				//
				//	fmt.Println("unmappedRead: ", unmappedRead)
				//	fmt.Println("unmappedGenome: ", unmappedGenome)
				//
				//	for i := 0; i < len(unmappedRead); i++ {
				//		if unmappedRead[i] != unmappedGenome[i] {
				//			mismatches = append(mismatches, i)
				//		}
				//	}
				//
				//	result.MatchedRead.AddRegion(startRead, endRead)
				//	result.MatchedGenome.AddRegion(startGenome, endGenome)
				//	result.MismatchesRead = append(result.MismatchesRead, mismatches...)
				//}
			}

			if result.MatchedRead.GetFirstRegion().Start > 0 {
				fmt.Println("UNMATCHED POSITIONS IN FRONT OF READ")

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
				fmt.Println("UNMATCHED POSITIONS IN BACK OF READ")

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

	/*
		for _, result := range results {

			if result.MatchedRead.Length() != len(*read.Sequence) {
				continue
			}

			qname := read.Header
			flag := 0
			rname := "tbd"
			pos := 45884425 + result.MatchedGenome.GetFirstRegion().Start
			mapq := 0
			cigar := strconv.Itoa(len(*read.Sequence)) + "M"
			rnext := "*"
			pnext := 0
			tlen := 0
			seq := string(*read.Sequence)
			qual := string(*read.Quality)

			//mm := fmt.Sprint(result.MismatchesRead)

			res := fmt.Sprintf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n", qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual)

			//fmt.Println(result.MatchedGenome.GetFirstRegion().Start)
			//fmt.Println(res)
			return res, true
		}

		return "unmappable\n", false
	*/
}
