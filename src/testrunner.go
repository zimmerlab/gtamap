package main

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapping"
	"github.com/KleinSamuel/gtamap/src/dataloader"
	"github.com/KleinSamuel/gtamap/src/datawriter"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/KleinSamuel/gtamap/src/formats/gtf"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
	"github.com/sirupsen/logrus"
	"os"
	"strconv"
	"strings"
	"sync"
	"time"
)

var KmerLength int = 8

func GetKmerLength() int {
	return KmerLength
}

func inspectAnnotation() {
	pathGtfZeroed := "../resources/ENSG00000173585.zeroed.gtf"
	pathFastaZeroed := "../resources/ENSG00000173585.fasta"
	pathFastaIndex := "../resources/ENSG00000173585.fasta.fai"

	var annotation *gtf.Annotation = dataloader.GenerateInputForIndexFromFile(pathGtfZeroed, pathFastaZeroed, pathFastaIndex)

	for _, trans := range annotation.Genes[0].Transcripts {
		fmt.Println(trans.SequenceDna)
	}
}

func buildAndSerializeIndex() {
	pathGtfZeroed := "../resources/ENSG00000173585.zeroed.gtf"
	pathFastaZeroed := "../resources/ENSG00000173585.zeroed.fasta"
	pathOutput := "../resources/ENSG00000173585.dev.gtai"

	//pathGtfZeroed := "/home/users/klein/Projects/GTA_Map/resources/ENSG00000173585.zeroed.gtf"
	//pathFastaZeroed := "/home/users/klein/Projects/GTA_Map/resources/ENSG00000173585.zeroed.fasta"
	//pathOutput := "/home/users/klein/Projects/GTA_Map/resources/ENSG00000173585.dev.gtai"

	gtfFile, errGtf := os.Open(pathGtfZeroed)
	if errGtf != nil {
		logrus.Fatal("Error reading gtf file", errGtf)
	}

	fastaFile, errFasta := os.Open(pathFastaZeroed)
	if errFasta != nil {
		logrus.Fatal("Error reading fasta file", errFasta)
	}

	outputFile, errIndex := os.Create(pathOutput)
	if errIndex != nil {
		logrus.Fatal("Error reading output file (.fai)", errIndex)
	}

	index.BuildAndSerializeIndex(gtfFile, fastaFile, outputFile)
}

func deserializeIndex() *index.GtaIndex {
	return index.DezerializeFromPath("../resources/ENSG00000173585.dev.gtai")
}

//func testArgparse() {
//
//	parser := argparse.NewParser("gtamap", "Gene-based Trancript-Aware readMAPping")
//
//	var cmdIndex *argparse.Command = parser.NewCommand("index", "Build the GTAMap index (.gtai).")
//
//	var gtfFile *os.File = cmdIndex.File("", "gtf", os.O_RDONLY, 0600, &argparse.Options{
//		Required: true,
//		Help:     "Gene transfer format (GTF) file.",
//	})
//	var fastaFile *os.File = cmdIndex.File("", "fasta", os.O_RDONLY, 0600, &argparse.Options{
//		Required: true,
//		Help:     "Nucleotide sequences (FASTA) file. The corresponding index (.fai) file must be present.",
//	})
//	var outputFile *os.File = cmdIndex.File("o", "output", os.O_WRONLY|os.O_CREATE, 0600, &argparse.Options{
//		Required: true,
//		Help:     "Output file (.gtai).",
//	})
//	var logLevelIndex *string = cmdIndex.Selector("", "loglevel", []string{"ERROR", "INFO", "DEBUG"}, &argparse.Options{
//		Required: false,
//		Help:     "Log output level.",
//		Default:  "ERROR",
//	})
//
//	var cmdMap *argparse.Command = parser.NewCommand("map", "Map reads to the GTAMap index.")
//	var indexFile *os.File = cmdMap.File("", "index", os.O_RDONLY, 0600, &argparse.Options{
//		Required: true,
//		Help:     "GTAMap index (.gtai) file",
//	})
//	var fastqFwFile *os.File = cmdMap.File("", "r1", os.O_RDONLY, 0600, &argparse.Options{
//		Required: true,
//		Help:     "FASTQ file containing the forward reads.",
//	})
//	var fastqRwFile *os.File = cmdMap.File("", "r2", os.O_RDONLY, 0600, &argparse.Options{
//		Required: false,
//		Help:     "FASTQ file containing the reverse reads.",
//	})
//	var logLevelMap *string = cmdMap.Selector("", "loglevel", []string{"ERROR", "INFO", "DEBUG"}, &argparse.Options{
//		Required: false,
//		Help:     "log output level",
//		Default:  "ERROR",
//	})
//
//	err := parser.Parse(os.Args)
//	if err != nil {
//		fmt.Print(parser.Usage(err))
//		return
//	}
//
//	if cmdIndex.Happened() {
//
//		fmt.Println("Indexing..")
//
//		fmt.Println("GTF file: ", *gtfFile)
//		fmt.Println("FASTA file: ", *fastaFile)
//		fmt.Println("Output file: ", *outputFile)
//		fmt.Println("loglevel: ", *logLevelIndex)
//
//	} else if cmdMap.Happened() {
//
//		fmt.Println("Mapping..")
//
//		fmt.Println("loglevel: ", *logLevelMap)
//		fmt.Println("index: ", *indexFile)
//		fmt.Println("fastq fw: ", *fastqFwFile)
//		fmt.Println("fastq rw: ", *fastqRwFile)
//
//	} else {
//		fmt.Println(parser.Usage("no valid command supplied"))
//		return
//	}
//}

func testFastqReader() {

	pathReadsFw := "../resources/reads_first_1.1.fq"
	//pathReadsRv := "../resources/reads_first_1.2.fq"

	reader := fastq.InitFromPaths(pathReadsFw, "")

	read := reader.NextRead()

	fmt.Println(read.ReadR1)
	fmt.Println(read.ReadR2)
}

func testMapping() {

	timerStartTotal := time.Now()

	timerStart := time.Now()

	gtaIndex := deserializeIndex()

	fmt.Println("read tree")
	fmt.Println("duration: ", time.Since(timerStart))

	samHeader := sam.Header{
		Version:                 sam.Version,
		ReferenceSequenceName:   gtaIndex.Gene.Chromosome,
		ReferenceSequenceLength: int(gtaIndex.Gene.EndGenomic - gtaIndex.Gene.StartGenomic + 1),
		GenomeAnnotationVersion: "",
		GenomeAssemblyVersion:   "",
		OrganismTaxId:           "",
		ToolVersion:             config.ToolVersion(),
		Transcripts:             make([]*sam.TranscriptInfo, len(gtaIndex.Transcripts)),
	}
	for i, transcript := range gtaIndex.Transcripts {
		samHeader.Transcripts[i] = &sam.TranscriptInfo{
			Id:                  "T" + strconv.Itoa(i),
			TranscriptEnsemblId: transcript.TranscriptIdEnsembl,
			TranscriptLength:    transcript.SequenceLength,
		}
	}

	timerStart = time.Now()

	pathReadsR1 := "../resources/reads/manual/reads_ccr9.1.fq"
	pathReadsR2 := "../resources/reads/manual/reads_ccr9.2.fq"

	//pathReadsR1 := "/usr/local/storage2/sam_fasta_test/fw.fastq"
	//pathReadsR2 := "/usr/local/storage2/sam_fasta_test/rw.fastq"

	reader := fastq.InitFromPaths(pathReadsR1, pathReadsR2)

	writer := datawriter.InitFromPath("../out/reads_ccr9.sam")
	writer.Write(samHeader.String())

	numWorkers := 1

	taskQueueMapping := make(chan ReadPairMappingTask)
	taskQueueWriter := make(chan string)

	// wait group that keeps track of the mapping goroutines that are still running
	var waitgroupMapping sync.WaitGroup

	// start the mapping worker goroutine pool
	for i := 0; i < numWorkers; i++ {
		waitgroupMapping.Add(1)
		go mapReadPairWorker(i, taskQueueMapping, taskQueueWriter, &waitgroupMapping, gtaIndex)
	}

	// wait group that keeps track of the writer goroutine that is still running
	var waitgroupWriter sync.WaitGroup
	// start the writer goroutine
	waitgroupWriter.Add(1)
	go writeOutputWorker(taskQueueWriter, &waitgroupWriter, writer)

	// the number of read pairs that have been processed
	taskCounter := 0

	for readPair := reader.NextRead(); readPair != nil; readPair = reader.NextRead() {

		// TODO: remove after testing (only process specific read pair)
		name := strings.Split(readPair.ReadR1.Header, " ")[0]
		if name != "@3-0002/1" {
			continue
		}

		mappingTask := ReadPairMappingTask{
			ID:       taskCounter,
			ReadPair: readPair,
		}
		taskQueueMapping <- mappingTask

		taskCounter++
	}

	close(taskQueueMapping)

	waitgroupMapping.Wait()

	fmt.Println("mapping finished")
	fmt.Println("num tasks: ", taskCounter)

	close(taskQueueWriter)

	waitgroupWriter.Wait()
	writer.Close()

	fmt.Println("writer finished")

	totalDuration := time.Since(timerStartTotal)

	logrus.WithFields(logrus.Fields{
		"duration": totalDuration,
	}).Info("Finished mapping")
}

type ReadPairMappingTask struct {
	ID       int
	ReadPair *fastq.ReadPair
}

func writeOutputWorker(taskQueue <-chan string, wg *sync.WaitGroup, writer *datawriter.Writer) {

	logrus.Info("Started writeOutputWorker")

	for task := range taskQueue {
		writer.Write(task)
	}

	wg.Done()

	logrus.Info("Finished writeOutputWorker")
}

func mapReadPairWorker(workerId int,
	taskQueue <-chan ReadPairMappingTask,
	taskQueueWriter chan<- string,
	wg *sync.WaitGroup,
	gtaIndex *index.GtaIndex) {

	logrus.WithFields(logrus.Fields{
		"workerId": workerId,
	}).Info("Started mapReadPairWorker")

	for task := range taskQueue {

		logrus.WithFields(logrus.Fields{
			"workerId": workerId,
			"task":     task.ID,
		}).Info("Processing task")

		// TODO: result is currently only dummy
		result := mapping.MapReadPairDev(task.ReadPair, gtaIndex)

		fmt.Println("result: ", result)

		taskQueueWriter <- result
	}

	// signal that the worker has finished
	wg.Done()

	logrus.WithFields(logrus.Fields{
		"workerId": workerId,
	}).Info("Finished mapReadPairWorker")
}

func main() {

	logrus.SetFormatter(&logrus.TextFormatter{
		FullTimestamp: true,
	})
	logrus.SetLevel(logrus.InfoLevel)

	//buildAndSerializeIndex()
	testMapping()

	//buildAndSerializeIndex()
	////testMapping()
	//
	//gtaIndex := deserializeIndex()
	//
	//leafs := gtaIndex.SuffixTree.FindLeafNodesRecursive(gtaIndex.SuffixTree.Root)
	//
	//numErrors := 0
	//
	//for _, leaf := range leafs {
	//
	//	//fmt.Println("leaf: ", leaf)
	//
	//	for sequenceIndex, locations := range leaf.Locations {
	//		lenTranscript := len(gtaIndex.GetTranscriptSequenceDna(sequenceIndex))
	//
	//		for _, location := range locations {
	//			if location > lenTranscript {
	//				fmt.Println("location: ", location)
	//				fmt.Println("lenTranscript: ", lenTranscript)
	//				numErrors++
	//				//panic("location > lenTranscript")
	//				fmt.Println(leaf)
	//			}
	//		}
	//	}
	//}
	//
	//gtaIndex.SuffixTree.ToEdgeList()
	//
	//fmt.Println("numErrors: ", numErrors)

	//allSequences := make([]string, 2)
	//
	//allSequences[0] = "abcdefabxybcdmnabcdex1"
	//allSequences[1] = "abcdefgh2"
	//
	//tree := datastructure.CreateNewTree()
	//
	//for i, sequence := range allSequences {
	//	tree.AddSequence(sequence, i)
	//}

	//logrus.Debug()
	//logrus.Debug("final tree:")
	//tree.ToEdgeList()

	//pattern := "e"
	//result := tree.FindPatternExact(&pattern)
	//fmt.Println(result)

	//datastructure.BuildSuffixTree(allSequences)

	//var suffixTree *datastructure.SuffixTree = datastructure.BuildSuffixTree(allSequences)

	//leafs := suffixTree.FindLeafNodesRecursive(suffixTree.Root)

	//for _, leaf := range leafs {
	//
	//	fmt.Println("leaf: ", leaf)
	//
	//	for sequenceIndex, locations := range leaf.Locations {
	//		lenTranscript := len(allSequences[sequenceIndex])
	//
	//		for _, location := range locations {
	//			if location > lenTranscript {
	//				fmt.Println("location: ", location)
	//				fmt.Println("lenTranscript: ", lenTranscript)
	//				panic("location > lenTranscript")
	//			}
	//		}
	//	}
	//}

	//pattern := "CA"
	//
	//result := suffixTree.FindPatternExact(&pattern)
	//
	//suffixTree.ToEdgeList()
	//
	//fmt.Println(result)
}
