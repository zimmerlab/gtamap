package main

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper"
	"github.com/KleinSamuel/gtamap/src/datawriter"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/akamensky/argparse"
	"github.com/sirupsen/logrus"
	"os"
	"strings"
)

func printBanner() {
	logrus.Info("┌─┐┌┬┐┌─┐┌┬┐┌─┐┌─┐")
	logrus.Info("│ ┬ │ ├─┤│││├─┤├─┘")
	logrus.Info("└─┘ ┴ ┴ ┴┴ ┴┴ ┴┴")
	logrus.Info("gtamap v" + config.ToolVersion() + " (S. Klein, M. Weyrich, 2025)")
	logrus.Info("Fast and memory efficient spliced read mapping to a single gene reference.")
	logrus.Info("")
}

func main() {

	parser := argparse.NewParser("gtamap", "Gene-centric spliced read mapping")

	var cmdIndexPre *argparse.Command = parser.NewCommand("index-pre", "Extract gene sequences from genome.")
	var fastaFileIndexPre *string = cmdIndexPre.String("", "fasta", &argparse.Options{
		Required: true,
		Help:     "Nucleotide sequences (FASTA) file (currently only non-compressed).",
	})
	var gtfFileIndexPre *string = cmdIndexPre.String("", "gtf", &argparse.Options{
		Required: true,
		Help:     "Genome annotation (GTF) file (currently only non-compressed).",
	})
	var outputDirIndexPre *string = cmdIndexPre.String("", "output", &argparse.Options{
		Required: true,
		Help:     "Output directory for extracted gene sequences.",
	})
	var geneIdsIndexPre *string = cmdIndexPre.String("", "geneids", &argparse.Options{
		Required: false,
		Help:     "Gene IDs to extract (comma-separated).",
	})
	var logLevelIndexPre *string = cmdIndexPre.Selector("", "loglevel", []string{"ERROR", "INFO", "DEBUG"}, &argparse.Options{
		Required: false,
		Help:     "Log output level.",
		Default:  "INFO",
	})

	var cmdIndex *argparse.Command = parser.NewCommand("index", "Build the GTAMap index (.gtai).")
	var fastaFile *os.File = cmdIndex.File("", "fasta", os.O_RDONLY, 0600, &argparse.Options{
		Required: true,
		Help:     "Nucleotide sequences (FASTA) file.",
	})
	var outputFileIndex *os.File = cmdIndex.File("", "output", os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0600, &argparse.Options{
		Required: true,
		Help:     "Output file (.gtai).",
	})
	var logLevelIndex *string = cmdIndex.Selector("", "loglevel", []string{"ERROR", "INFO", "DEBUG"}, &argparse.Options{
		Required: false,
		Help:     "Log output level.",
		Default:  "INFO",
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
	var outputFileMap *os.File = cmdMap.File("", "output", os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0600, &argparse.Options{
		Required: true,
		Help:     "Output file (.gtai).",
	})
	var logLevelMap *string = cmdMap.Selector("", "loglevel", []string{"ERROR", "INFO", "DEBUG"}, &argparse.Options{
		Required: false,
		Help:     "log output level",
		Default:  "INFO",
	})
	var numThreads *int = cmdMap.Int("", "threads", &argparse.Options{
		Required: false,
		Help:     "Number of threads to use (default: all)",
		Default:  -1,
	})

	err := parser.Parse(os.Args)
	if err != nil {
		fmt.Print(parser.Usage(err))
		return
	}

	logrus.SetFormatter(&logrus.TextFormatter{
		FullTimestamp: true,
	})

	if cmdIndexPre.Happened() {

		level, _ := logrus.ParseLevel(*logLevelIndexPre)
		logrus.SetLevel(level)

		printBanner()
		logrus.Info("Extracting gene sequences from genome")

		// parsing gene ids
		geneIds := make(map[string]struct{})
		if *geneIdsIndexPre != "" {
			genes := strings.Split(*geneIdsIndexPre, ",")
			for _, gene := range genes {
				geneIds[gene] = struct{}{}
			}
		}

		index.ExtractGeneSequenceFromGtfAndFastaForIndex(*gtfFileIndexPre, *fastaFileIndexPre,
			*outputDirIndexPre, geneIds)

		logrus.Info("Done")

	} else if cmdIndex.Happened() {

		level, _ := logrus.ParseLevel(*logLevelIndex)
		logrus.SetLevel(level)

		printBanner()
		logrus.Info("Building GTAMap index (.gtai)")

		index.BuildAndSerializeGenomeIndex(fastaFile, outputFileIndex)

	} else if cmdMap.Happened() {

		level, _ := logrus.ParseLevel(*logLevelMap)
		logrus.SetLevel(level)

		printBanner()
		logrus.Info("Mapping reads to given index")

		genomeIndex := index.ReadGenomeIndexByFile(indexFile)

		reader := fastq.InitFromFiles(fastqFwFile, fastqRwFile)

		writer := datawriter.InitFromFile(outputFileMap)

		config.OutSAMFile = outputFileMap.Name()

		mapper.MapAll(genomeIndex, reader, writer, numThreads)

	} else {
		fmt.Println(parser.Usage("no valid command supplied"))
		return
	}
}
