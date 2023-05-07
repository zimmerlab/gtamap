package main

import (
	"github.com/KleinSamuel/gtamap/src/dataloader"
	"github.com/KleinSamuel/gtamap/src/datastructure"
	"github.com/KleinSamuel/gtamap/src/gtf"
	"github.com/sirupsen/logrus"
)

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
	//sequences := []string{"AAAABAAAABAAC"}

	tree := datastructure.BuildSuffixTree(sequences)

	//tree.PrintEdgeList()

	tree.Search("ATGACACCCACAGACTT")
}
