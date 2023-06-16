package main

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/akamensky/argparse"
	"github.com/sirupsen/logrus"
	"os"
)

func printBanner() {
	logrus.Info(" _____ _____ ___   ___  ___  ___  ______")
	logrus.Info("|  __ \\_   _/ _ \\  |  \\/  | / _ \\ | ___ \\")
	logrus.Info("| |  \\/ | |/ /_\\ \\ | .  . |/ /_\\ \\| |_/ /")
	logrus.Info("| | __  | ||  _  | | |\\/| ||  _  ||  __/")
	logrus.Info("| |_\\ \\ | || | | | | |  | || | | || |")
	logrus.Info(" \\____/ \\_/\\_| |_/ \\_|  |_/\\_| |_/\\_|")
	logrus.Info("GTAMap v" + config.ToolVersion() + " (A. Hadziahmetovic, S. Klein, 2023)")
	logrus.Info("Fast (and memory efficient) RNA-seq read mapping to transcripts of a single gene.")
	logrus.Info("")
}

func main() {

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

	logrus.SetFormatter(&logrus.TextFormatter{
		FullTimestamp: true,
	})

	if cmdIndex.Happened() {

		level, _ := logrus.ParseLevel(*logLevelIndex)
		logrus.SetLevel(level)

		printBanner()
		logrus.Info("Building GTAMap index (.gtai) ..")

		index.BuildAndSerializeIndex(gtfFile, fastaFile, outputFile)

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
