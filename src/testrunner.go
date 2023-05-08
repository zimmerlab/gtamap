package main

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core"
	"github.com/KleinSamuel/gtamap/src/core/datastructure"
	"github.com/KleinSamuel/gtamap/src/dataloader"
	"github.com/KleinSamuel/gtamap/src/dataloader/gtf"
	"github.com/sirupsen/logrus"
)

func serializeExample() {
	sequences := []string{"ABCABCABCD"}

	tree := datastructure.BuildSuffixTree(sequences)

	fmt.Println("serializing tree..")

	datastructure.SerializeSuffixTree(tree, "../resources/tree")

	fmt.Println("done serializing tree")

	fmt.Println("reading tree..")

	var newTree *datastructure.SuffixTree = datastructure.DezerializeSuffixTree("../resources/tree")

	pattern := "ABCD"

	var result *core.PatternSearchResult = newTree.Search(&pattern)

	fmt.Println("result: ", result)
}

func main() {

	logrus.SetFormatter(&logrus.TextFormatter{
		FullTimestamp: true,
	})
	logrus.SetLevel(logrus.InfoLevel)

	pathGtfZeroed := "../resources/ENSG00000173585.zeroed.gtf"
	pathFastaZeroed := "../resources/ENSG00000173585.fasta"

	var annotation *gtf.Annotation = dataloader.GenerateInputForIndex(pathGtfZeroed, pathFastaZeroed)

	sequences := []string{annotation.Genes[0].Transcripts[0].SequenceDna}

	// ATGACACCCACAGACTT
	// ABCXABCY
	//sequences := []string{"ABCABCABCD"}

	tree := datastructure.BuildSuffixTree(sequences)
	//tree.PrintEdgeList()

	datastructure.SerializeSuffixTree(tree, "../resources/ENSG00000173585.gtai")

	/*
		pattern := "ABCD"

		var result *core.PatternSearchResult = tree.Search(&pattern)

		fmt.Println("result: ", result)
	*/
}
