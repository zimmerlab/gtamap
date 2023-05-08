package main

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core"
	"github.com/KleinSamuel/gtamap/src/core/datastructure"
	"github.com/KleinSamuel/gtamap/src/dataloader"
	"github.com/KleinSamuel/gtamap/src/dataloader/gtf"
	"github.com/sirupsen/logrus"
	"time"
)

func inspectAnnotation() {
	pathGtfZeroed := "../resources/ENSG00000173585.zeroed.gtf"
	pathFastaZeroed := "../resources/ENSG00000173585.fasta"

	var annotation *gtf.Annotation = dataloader.GenerateInputForIndex(pathGtfZeroed, pathFastaZeroed)

	for _, trans := range annotation.Genes[0].Transcripts {
		fmt.Println(trans.SequenceDna)
	}
}

func buildAndSerializeSuffixTree() {
	pathGtfZeroed := "../resources/ENSG00000173585.zeroed.gtf"
	pathFastaZeroed := "../resources/ENSG00000173585.fasta"

	var annotation *gtf.Annotation = dataloader.GenerateInputForIndex(pathGtfZeroed, pathFastaZeroed)
	sequences := []string{annotation.Genes[0].Transcripts[0].SequenceDna}

	tree := datastructure.BuildSuffixTree(sequences)

	datastructure.SerializeSuffixTree(tree, "../resources/ENSG00000173585.gtai")
}

func deserializeSuffixTree() *datastructure.SuffixTree {
	return datastructure.DezerializeSuffixTree("../resources/ENSG00000173585.gtai")
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
	tree := deserializeSuffixTree()

	timerStart := time.Now()

	pattern := "ATGA"

	var result *core.PatternSearchResult = tree.Search(&pattern)

	fmt.Println("duration: ", time.Since(timerStart))
	fmt.Println("result: ", result)
}

func main() {

	logrus.SetFormatter(&logrus.TextFormatter{
		FullTimestamp: true,
	})
	logrus.SetLevel(logrus.InfoLevel)

	buildAndSerializeSuffixTree()
}
