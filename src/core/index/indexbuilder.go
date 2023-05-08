package index

import (
	"github.com/KleinSamuel/gtamap/src/core/datastructure"
	"github.com/KleinSamuel/gtamap/src/dataloader"
	"github.com/KleinSamuel/gtamap/src/dataloader/gtf"
	"github.com/sirupsen/logrus"
	"log"
	"os"
	"time"
)

func BuildAndSerializeIndex(gtfFile *os.File, fastaFile *os.File, outputFile *os.File) {

	timerStart := time.Now()

	var fastIndexFilePath string = fastaFile.Name() + ".fai"
	dataloader.ExitIfFastaIndexIsMissing(fastIndexFilePath)
	fastaIndexFile, err := os.Open(fastIndexFilePath)
	if err != nil {
		log.Fatal("Error reading fasta index file (.fai)", err)
	}

	var annotation *gtf.Annotation = dataloader.GenerateInputForIndex(gtfFile, fastaFile, fastaIndexFile)

	var sequences []string = make([]string, len(annotation.Genes[0].Transcripts))
	for i, transcript := range annotation.Genes[0].Transcripts {
		sequences[i] = transcript.SequenceDna
	}

	var suffixTree *datastructure.SuffixTree = datastructure.BuildSuffixTree(sequences)

	datastructure.SerializeSuffixTree(suffixTree, outputFile)

	logrus.WithFields(logrus.Fields{
		"total": time.Since(timerStart).String(),
	}).Info("Successfully built GTAMap index.")
}
