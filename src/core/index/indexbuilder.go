package index

import (
	"fmt"
	"os"
	"path/filepath"

	"github.com/KleinSamuel/gtamap/src/dataloader"
	"github.com/KleinSamuel/gtamap/src/formats/fasta"
	"github.com/KleinSamuel/gtamap/src/formats/gtf"
	"github.com/sirupsen/logrus"
)

func ExtractSequenceFromFastaForIndex(fastaPath string, chromosome string, start int, end int, outputPath string) {
	fastaIndexPath := fastaPath + ".fai"

	logrus.WithFields(logrus.Fields{
		"fastaPath":      fastaPath,
		"fastaIndexPath": fastaIndexPath,
		"chromosome":     chromosome,
		"start":          start,
		"end":            end,
		"outputPath":     outputPath,
	}).Info("Extracting sequence from fasta")

	// check if output directory exists
	if _, err := os.Stat(outputPath); os.IsNotExist(err) {
		err := os.MkdirAll(outputPath, os.ModePerm)
		if err != nil {
			logrus.Fatal("Error creating output directory", err)
		}
	}

	fastaIndex := fasta.ReadFastaIndexUsingPath(fastaIndexPath)
	if fastaIndex == nil {
		panic("Could not read fasta index: " + fastaIndexPath)
	}

	fastaFile, fastaFileErr := os.Open(fastaPath)
	if fastaFileErr != nil {
		panic(fastaFileErr)
	}
	defer fastaFile.Close()

	name := fmt.Sprintf("%s:%d-%d", chromosome, start, end)

	outFilePath := filepath.Join(outputPath, name+".fa")

	outFile, errOpen := os.Create(outFilePath)
	if errOpen != nil {
		logrus.Fatal("Error creating output file (.fa)", errOpen)
	}
	defer outFile.Close()

	logrus.WithFields(logrus.Fields{
		"outputFile": outFilePath,
	}).Info("Creating output file")

	seq := dataloader.ExtractSequenceAsStringFromFasta(fastaFile, fastaIndex, chromosome, uint32(start), uint32(end))

	AppendSequenceToFastaFile(name, chromosome, true, start, end, []byte(seq), 60, outFile)
}

func ExtractGeneSequenceFromGtfAndFastaForIndex(gtfPath string, fastaPath string, outputPath string,
	geneIds map[string]struct{}, upstreamBases int, downstreamBases int, separateExtraction bool,
) {
	fastaIndexPath := fastaPath + ".fai"

	logrus.WithFields(logrus.Fields{
		"gtfPath":         gtfPath,
		"fastaPath":       fastaPath,
		"fastaIndexPath":  fastaIndexPath,
		"outputPath":      outputPath,
		"numGenes":        len(geneIds),
		"upstreamBases":   upstreamBases,
		"downstreamBases": downstreamBases,
	}).Info("Extracting gene sequence")

	// check if output directory exists
	if _, err := os.Stat(outputPath); os.IsNotExist(err) {
		err := os.MkdirAll(outputPath, os.ModePerm)
		if err != nil {
			logrus.Fatal("Error creating output directory", err)
		}
	}

	geneInfo := gtf.ReadGenesFromGtfUsingPath(gtfPath, geneIds)
	if len(geneInfo) == 0 {
		logrus.WithFields(logrus.Fields{
			"provided geneIds": geneIds,
			"gtf":              gtfPath,
		}).Warnf("Non of the target geneIds were found in provided gtf")
		return
	}

	fastaIndex := fasta.ReadFastaIndexUsingPath(fastaIndexPath)
	if fastaIndex == nil {
		panic("Could not read fasta index: " + fastaIndexPath)
	}

	fastaFile, fastaFileErr := os.Open(fastaPath)
	if fastaFileErr != nil {
		panic(fastaFileErr)
	}
	defer fastaFile.Close()

	if separateExtraction {
		// write each gene to separate fa file
		for _, gene := range geneInfo {
			if gene == nil {
				logrus.Warn("Skipping gene since it was not found in provided gtf.")
				continue
			}

			outFilePath := filepath.Join(outputPath, gene.GeneId+".fa")
			outFile, errOpen := os.Create(outFilePath)
			if errOpen != nil {
				logrus.Fatal("Error creating separate output file (.fa)", errOpen)
			}

			logrus.WithFields(logrus.Fields{
				"separated fasta": outputPath,
			}).Info("Creating a separate fasta file")

			geneStart := int(gene.StartGenomic) - upstreamBases
			if geneStart < 0 {
				geneStart = 0

				logrus.WithFields(logrus.Fields{
					"geneId":   gene.GeneId,
					"start":    gene.StartGenomic,
					"upstream": upstreamBases,
				}).Warn("Genomic region is negative, setting to 0")
			}

			geneEnd := int(gene.EndGenomic) + downstreamBases

			seq := dataloader.ExtractSequenceAsStringFromFasta(fastaFile, fastaIndex, gene.Contig,
				uint32(geneStart), uint32(geneEnd))

			AppendSequenceToFastaFile(gene.GeneId, gene.Contig, gene.IsForwardStrand, geneStart, geneEnd,
				[]byte(seq), 60, outFile)

			outFile.Close()
		}
	} else {
		// create shared fasta file
		outFilePath := ""
		if len(geneInfo) > 1 {
			// if several gene ids are in geneInfo, name shared .fa file genes.fa
			outFilePath = filepath.Join(outputPath, "genes.fa")
		} else {
			// if there's only one element in geneInfo and --splitgenes was not set,
			// use the gene id as file name
			outFilePath = filepath.Join(outputPath, geneInfo[0].GeneId+".fa")
		}
		outFile, errOpen := os.Create(outFilePath)
		if errOpen != nil {
			logrus.Fatal("Error creating shared output file (.fa)", errOpen)
		}

		logrus.WithFields(logrus.Fields{
			"shared fasta": outputPath,
		}).Info("Creating a shared fasta file")

		for _, gene := range geneInfo {

			geneStart := int(gene.StartGenomic) - upstreamBases
			if geneStart < 0 {
				geneStart = 0

				logrus.WithFields(logrus.Fields{
					"geneId":   gene.GeneId,
					"start":    gene.StartGenomic,
					"upstream": upstreamBases,
				}).Warn("Genomic region is negative, setting to 0")
			}

			geneEnd := int(gene.EndGenomic) + downstreamBases

			seq := dataloader.ExtractSequenceAsStringFromFasta(fastaFile, fastaIndex, gene.Contig,
				uint32(geneStart), uint32(geneEnd))

			AppendSequenceToFastaFile(gene.GeneId, gene.Contig, gene.IsForwardStrand, geneStart, geneEnd,
				[]byte(seq), 60, outFile)
		}
		outFile.Close()
	}
	return
}

func AppendSequenceToFastaFile(name string, contig string, isForwardStrand bool, startGenomic int, endGenomic int,
	sequence []byte, charPerLine int, outFile *os.File,
) {
	strand := "+"
	if !isForwardStrand {
		strand = "-"
	}

	_, errWrite := outFile.WriteString(fmt.Sprintf(">%s\t%s\t%s\t%d\t%d\n",
		contig, name, strand, startGenomic, endGenomic))

	if errWrite != nil {
		logrus.Fatal("Error writing to output file", errWrite)
	}

	// write sequence in lines of charPerLine characters
	for i := 0; i < len(sequence); i += charPerLine {
		end := i + charPerLine
		if end > len(sequence) {
			end = len(sequence)
		}
		_, errWrite := outFile.WriteString(string(sequence[i:end]) + "\n")
		if errWrite != nil {
			logrus.Fatal("Error writing to output file", errWrite)
		}
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
		gene.Contig, gene.GeneId, strand, gene.StartGenomic, gene.EndGenomic, seq))

	if errWrite != nil {
		logrus.Fatal("Error writing to output file", errWrite)
	}
}
