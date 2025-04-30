package index

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/dataloader"
	"github.com/KleinSamuel/gtamap/src/formats/fasta"
	"github.com/KleinSamuel/gtamap/src/formats/gtf"
	"github.com/sirupsen/logrus"
	"os"
	"path/filepath"
)

func ExtractGeneSequenceFromGtfAndFastaForIndex(gtfPath string, fastaPath string, outputPath string,
	geneIds map[string]struct{}) {

	fastaIndexPath := fastaPath + ".fai"

	logrus.WithFields(logrus.Fields{
		"gtfPath":        gtfPath,
		"fastaPath":      fastaPath,
		"fastaIndexPath": fastaIndexPath,
		"outputPath":     outputPath,
		"numGenes":       len(geneIds),
	}).Info("Extracting gene sequence")

	// check if output directory exists
	if _, err := os.Stat(outputPath); os.IsNotExist(err) {
		err := os.MkdirAll(outputPath, os.ModePerm)
		if err != nil {
			logrus.Fatal("Error creating output directory", err)
		}
	}

	geneInfo := gtf.ReadGenesFromGtfUsingPath(gtfPath, geneIds)

	fastaIndex := fasta.ReadFastaIndexUsingPath(fastaIndexPath)
	if fastaIndex == nil {
		panic("Could not read fasta index")
	}

	fastaFile, fastaFileErr := os.Open(fastaPath)
	if fastaFileErr != nil {
		panic(fastaFileErr)
	}
	defer fastaFile.Close()

	outFilePath := filepath.Join(outputPath, "sequences.fa")
	outFile, errOpen := os.Create(outFilePath)
	defer outFile.Close()

	if errOpen != nil {
		logrus.Fatal("Error creating output file (.fai)", errOpen)
	}

	for _, gene := range geneInfo {

		logrus.WithFields(logrus.Fields{
			"geneId": gene.GeneId,
		}).Info("Extracting sequence information")

		seq := dataloader.ExtractSequenceAsStringFromFasta(fastaFile, fastaIndex, gene.Contig,
			gene.StartGenomic, gene.EndGenomic)

		strand := "+"
		if !gene.IsForwardStrand {
			strand = "-"
		}

		_, errWrite := outFile.WriteString(fmt.Sprintf(">%s\t%s\t%s\t%d\t%d\n%s\n",
			gene.GeneId, gene.Contig, strand, gene.StartGenomic, gene.EndGenomic, seq))

		if errWrite != nil {
			logrus.Fatal("Error writing to output file", errWrite)
		}

	}
}
