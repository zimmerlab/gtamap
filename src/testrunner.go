package main

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core"
	"github.com/KleinSamuel/gtamap/src/core/datastructure"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapping"
	"github.com/KleinSamuel/gtamap/src/dataloader"
	"github.com/KleinSamuel/gtamap/src/dataloader/fastq"
	"github.com/KleinSamuel/gtamap/src/dataloader/gtf"
	"github.com/akamensky/argparse"
	"github.com/sirupsen/logrus"
	"os"
	"time"
)

func inspectAnnotation() {
	pathGtfZeroed := "../resources/ENSG00000173585.zeroed.gtf"
	pathFastaZeroed := "../resources/ENSG00000173585.fasta"
	pathFastaIndex := "../resources/ENSG00000173585.fasta.fai"

	var annotation *gtf.Annotation = dataloader.GenerateInputForIndexFromFile(pathGtfZeroed, pathFastaZeroed, pathFastaIndex)

	for _, trans := range annotation.Genes[0].Transcripts {
		fmt.Println(trans.SequenceDna)
	}
}

/*
func buildAndSerializeIndex() {
	pathGtfZeroed := "../resources/ENSG00000173585.zeroed.gtf"
	pathFastaZeroed := "../resources/ENSG00000173585.fasta"
	pathFastaIndex := "../resources/ENSG00000173585.fasta.fai"

	var annotation *gtf.Annotation = dataloader.GenerateInputForIndexFromFile(pathGtfZeroed, pathFastaZeroed, pathFastaIndex)

	var sequences []string = []string{annotation.Genes[0].Transcripts[0].SequenceDna}

	tree := datastructure.BuildSuffixTree(sequences)

	datastructure.SerializeSuffixTreeFromFile(tree, "../resources/ENSG00000173585.gtai")
}
*/

func deserializeIndex() *index.GtaIndex {
	return index.DezerializeFromPath("../resources/ENSG00000173585.gtai")
}

func buildTestTree() {
	sequences := []string{"ABCABCD", "CABCDA"}

	tree := datastructure.BuildSuffixTree(sequences)

	//tree.ToEdgeList()

	timerStart := time.Now()

	pattern := "C"

	var result *core.PatternSearchResult = tree.Search(&pattern)

	fmt.Println("duration: ", time.Since(timerStart))
	fmt.Println("result: ", result)
}

func search() {
	gtaIndex := deserializeIndex()

	timerStart := time.Now()

	pattern := "ATGA"

	var result *core.PatternSearchResult = gtaIndex.SuffixTreeForwardStrandForwardDirection.Search(&pattern)

	fmt.Println("duration: ", time.Since(timerStart))
	fmt.Println("result: ", result)
}

func testArgparse() {

	parser := argparse.NewParser("gtamap", "Gene-based Trancript-Aware readMAPping")

	var cmdIndex *argparse.Command = parser.NewCommand("index", "Build the GTAMap index (.gtai).")

	var gtfFile *os.File = cmdIndex.File("", "gtf", os.O_RDONLY, 0600, &argparse.Options{
		Required: true,
		Help:     "Gene transfer format (GTF) file.",
	})
	var fastaFile *os.File = cmdIndex.File("", "fasta", os.O_RDONLY, 0600, &argparse.Options{
		Required: true,
		Help:     "Nucleotide sequences (FASTA) file. The corresponding index (.fai) file must be present.",
	})
	var outputFile *os.File = cmdIndex.File("o", "output", os.O_WRONLY|os.O_CREATE, 0600, &argparse.Options{
		Required: true,
		Help:     "Output file (.gtai).",
	})
	var logLevelIndex *string = cmdIndex.Selector("", "loglevel", []string{"ERROR", "INFO", "DEBUG"}, &argparse.Options{
		Required: false,
		Help:     "Log output level.",
		Default:  "ERROR",
	})

	var cmdMap *argparse.Command = parser.NewCommand("map", "Map reads to the GTAMap index.")
	var indexFile *os.File = cmdMap.File("", "index", os.O_RDONLY, 0600, &argparse.Options{
		Required: true,
		Help:     "GTAMap index (.gtai) file",
	})
	var fastqFwFile *os.File = cmdMap.File("", "r1", os.O_RDONLY, 0600, &argparse.Options{
		Required: true,
		Help:     "FASTQ file containing the forward reads.",
	})
	var fastqRwFile *os.File = cmdMap.File("", "r2", os.O_RDONLY, 0600, &argparse.Options{
		Required: false,
		Help:     "FASTQ file containing the reverse reads.",
	})
	var logLevelMap *string = cmdMap.Selector("", "loglevel", []string{"ERROR", "INFO", "DEBUG"}, &argparse.Options{
		Required: false,
		Help:     "log output level",
		Default:  "ERROR",
	})

	err := parser.Parse(os.Args)
	if err != nil {
		fmt.Print(parser.Usage(err))
		return
	}

	if cmdIndex.Happened() {

		fmt.Println("Indexing..")

		fmt.Println("GTF file: ", *gtfFile)
		fmt.Println("FASTA file: ", *fastaFile)
		fmt.Println("Output file: ", *outputFile)
		fmt.Println("loglevel: ", *logLevelIndex)

	} else if cmdMap.Happened() {

		fmt.Println("Mapping..")

		fmt.Println("loglevel: ", *logLevelMap)
		fmt.Println("index: ", *indexFile)
		fmt.Println("fastq fw: ", *fastqFwFile)
		fmt.Println("fastq rw: ", *fastqRwFile)

	} else {
		fmt.Println(parser.Usage("no valid command supplied"))
		return
	}
}

func testFastqReader() {

	pathReadsFw := "../resources/reads_first_1.1.fq"
	//pathReadsRv := "../resources/reads_first_1.2.fq"

	reader := fastq.InitFromPaths(pathReadsFw, "")

	read := reader.NextRead()

	fmt.Println(read.ReadR1)
	fmt.Println(read.ReadR2)
}

func testMapping() {

	timerStart := time.Now()

	gtaIndex := deserializeIndex()

	fmt.Println("read tree")
	fmt.Println("duration: ", time.Since(timerStart))

	timerStart = time.Now()

	pathReadsR1 := "../resources/reads_ccr9.1.fq"
	pathReadsR2 := "../resources/reads_ccr9.2.fq"

	reader := fastq.InitFromPaths(pathReadsR1, pathReadsR2)

	for read := reader.NextRead(); read != nil; read = reader.NextRead() {

		timerStart = time.Now()

		mapping.MapReadPair(read, gtaIndex)

		fmt.Println("map read pair duration: ", time.Since(timerStart))

		break
	}

}

func main() {

	logrus.SetFormatter(&logrus.TextFormatter{
		FullTimestamp: true,
	})
	logrus.SetLevel(logrus.InfoLevel)

	testMapping()
}
