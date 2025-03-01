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
	"runtime"
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

		result := MapReadPair(task.ReadPair, genomeIndex, timerChannel)

		taskQueueWriter <- result
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
		if name != "@3-0000/6-fw" {
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

func MapReadPair(readPair *fastq.ReadPair, genomeIndex *index.GenomeIndex, timerChannel chan<- *timer.Timer) string {

	keepFw := Filter(readPair.ReadR1.Sequence, genomeIndex)
	keepRw := Filter(readPair.ReadR2.Sequence, genomeIndex)

	if !keepFw || !keepRw {
		return "DISCARD\n"
	}

	FindAnchors(readPair, genomeIndex)

	return "result\n"
}

// Filter the read sequence by checking if it contains at least 5 matching kmers
// TODO: replace this by map
func Filter(readSequence *[]byte, genomeIndex *index.GenomeIndex) bool {

	numMatching := 0

	for i := 0; i <= len(*readSequence)-(int(config.KmerLength())); i += int(config.KmerLength()) {

		kmer := (*readSequence)[i : i+int(config.KmerLength())]

		matches := genomeIndex.KeywordTree.FindKeyword(&kmer, i)

		if matches != nil {
			numMatching++
		}
	}

	return numMatching >= 5
}

func FindAnchors(pair *fastq.ReadPair, genomeIndex *index.GenomeIndex) {

	fmt.Println("FindAnchors")

	//diagonals := make(map[int]int)

	globalMatches := matchutils.GlobalMatchResult{
		MatchesPerSequence: make([]*matchutils.SequenceMatchResult, genomeIndex.KeywordTree.NumSequences),
	}

	// create all non-overlapping k-mers for the read pair
	for i := 0; i <= len(*pair.ReadR1.Sequence)-(int(config.KmerLength())); i += int(config.KmerLength()) {

		kmer := (*pair.ReadR1.Sequence)[i : i+int(config.KmerLength())]

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

		for result.MatchedRead.Length() < len(*pair.ReadR1.Sequence) {

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

					readByte := (*pair.ReadR1.Sequence)[i]
					gIndex := diagonalGenome.GetFirstRegion().Start + i
					genomeByte := (*genomeIndex.Sequence)[gIndex]

					if readByte != genomeByte {
						mismatches = append(mismatches, i)

						logrus.Info("mismatch at read position: ", i)
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

		if result.MatchedRead.Length() == len(*pair.ReadR1.Sequence) {
			fmt.Println("MATCHED FULLY")
			break
		}

		if len(result.MatchedRead.Regions) > 1 {
			fmt.Println("UNMATCHED POSITIONS WITHIN READ")
		}

		if result.MatchedRead.GetFirstRegion().Start > 0 {
			fmt.Println("UNMATCHED POSITIONS IN FRONT OF READ")
		}

		if result.MatchedRead.GetLastRegion().End < len(*pair.ReadR1.Sequence) {
			fmt.Println("UNMATCHED POSITIONS IN BACK OF READ")
		}

		break
	}
}
