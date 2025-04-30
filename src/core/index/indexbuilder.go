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
	geneIds map[string]struct{}, separateExtraction bool) {

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

	if separateExtraction {
		// write each gene to separate fa file
		for _, gene := range geneInfo {
			outFilePath := filepath.Join(outputPath, gene.GeneId+".fa")
			outFile, errOpen := os.Create(outFilePath)

			if errOpen != nil {
				logrus.Fatal("Error creating separate output file (.fa)", errOpen)
			}

			logrus.WithFields(logrus.Fields{
				"separated fasta": outputPath,
			}).Info("Creating a separate fasta file")

			writeToFasta(gene, outFile, fastaFile, fastaIndex)
			outFile.Close()
		}
	} else {
		// create shared fasta file
		outFilePath := filepath.Join(outputPath, "sequences.fa")
		outFile, errOpen := os.Create(outFilePath)
		if errOpen != nil {
			logrus.Fatal("Error creating shared output file (.fa)", errOpen)
		}

		logrus.WithFields(logrus.Fields{
			"shared fasta": outputPath,
		}).Info("Creating a shared fasta file")

		for _, gene := range geneInfo {
			writeToFasta(gene, outFile, fastaFile, fastaIndex)
		}
		outFile.Close()
	}
}

func writeToFasta(gene *gtf.GeneBasic, outFile *os.File, fastaFile *os.File, fastaIndex *fasta.Index) {
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
