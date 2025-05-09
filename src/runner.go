package main

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/extraction"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper"
	"github.com/KleinSamuel/gtamap/src/datawriter"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/akamensky/argparse"
	"github.com/sirupsen/logrus"
	"os"
	"path/filepath"
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
	var separateExtraction *bool = cmdIndexPre.Flag("", "splitgenes", &argparse.Options{
		Help: "Extract gene sequences into separate fasta files, if more than one gene id is specified in --geneids",
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
	var paralogeFilePathMap *string = cmdMap.String("", "paraloges", &argparse.Options{
		Required: false,
		Help:     "Paraloge region meta file for target regions.",
		Default:  "",
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

	var cmdParalogePre *argparse.Command = parser.NewCommand("paraloge", "Extract known paraloge genes from ENSEMBL Database and prepare paraloge.csv for main target index extension.")
	var geneIdsParalogePre *string = cmdParalogePre.String("", "geneids", &argparse.Options{
		Required: false,
		Help:     "Query Gene IDs for extracting paraloge genes from DB (comma-separated).",
	})
	var indexDirParalogePre *string = cmdParalogePre.String("", "indexdir", &argparse.Options{
		Required: false,
		Help:     "Target directory for .gai index files of paraloge genes.",
		Default:  "index",
	})
	var fastaDirParalogePre *string = cmdParalogePre.String("", "fastadir", &argparse.Options{
		Required: false,
		Help:     "Target directory for .fa files of paraloge genes.",
		Default:  "fasta_in",
	})
	var speciesParalogePre *string = cmdParalogePre.String("", "species", &argparse.Options{
		Required: false,
		Help:     "Corresponding species of target genes. Species name needs to exist in https://www.ensembl.org/index.html",
		Default:  "human",
	})
	var fastaFileParalogePre *string = cmdParalogePre.String("", "fasta", &argparse.Options{
		Required: true,
		Help:     "Nucleotide sequences (FASTA) file (currently only non-compressed). Note: A corresponding .fai file needs to be in the same dir as the genome.fa file.",
	})
	var gtfFileParalogePre *string = cmdParalogePre.String("", "gtf", &argparse.Options{
		Required: true,
		Help:     "Genome annotation (GTF) file (currently only non-compressed).",
	})
	var paralogeMetaOutParalogePre *string = cmdParalogePre.String("", "meta", &argparse.Options{
		Required: false,
		Help:     "Genome annotation (GTF) file (currently only non-compressed).",
		Default:  "paraloges.csv",
	})
	var logLevelParalogePre *string = cmdParalogePre.Selector("", "loglevel", []string{"ERROR", "INFO", "DEBUG"}, &argparse.Options{
		Required: false,
		Help:     "Log output level.",
		Default:  "INFO",
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
			*outputDirIndexPre, geneIds, *separateExtraction)

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
		if *paralogeFilePathMap != "" {
			paralogeFileMap, err := os.Open(*paralogeFilePathMap)
			if err != nil {
				panic("Error reading provided paraloge file. Make sure it exists")
			}
			genomeIndex.LoadParaloges(paralogeFileMap)
		}

		reader := fastq.InitFromFiles(fastqFwFile, fastqRwFile)

		writer := datawriter.InitFromFile(outputFileMap)

		mapper.MapAll(genomeIndex, reader, writer, numThreads)

	} else if cmdParalogePre.Happened() {
		level, _ := logrus.ParseLevel(*logLevelParalogePre)
		logrus.SetLevel(level)

		printBanner()
		logrus.Info("Scanning DB for paraloges of specified target IDs.")

		// parsing gene ids
		targetGeneIds := make([]string, 0)
		if *geneIdsParalogePre != "" {
			genes := strings.Split(*geneIdsParalogePre, ",")
			for _, gene := range genes {
				targetGeneIds = append(targetGeneIds, gene)
			}
		}
		// I. get paraloge seqs per target gene
		targetParaloges := extraction.GetParaloges(targetGeneIds, *speciesParalogePre)

		// II. extract all target genes into separate .fa files in '--fastadir'
		for target, paraloges := range targetParaloges {
			logrus.Infof("Extracting sequences for paraloges of %s", target)
			index.ExtractGeneSequenceFromGtfAndFastaForIndex(*gtfFileParalogePre, *fastaFileParalogePre,
				*fastaDirParalogePre, paraloges, true)
		}

		// III. read in all seqs in '--fastadir' and serialize index into '--indexdir' for each paraloge seq
		indexPaths := make(map[string][]string)
		for target, paraloges := range targetParaloges {
			for paraloge, _ := range paraloges {
				err := os.MkdirAll(*indexDirParalogePre, 0755)
				if err != nil {
					logrus.Fatalf("Something went wrong creating dir %s: %s", *indexDirParalogePre, err)
				}

				// create and format index file
				paralogeIndexName := fmt.Sprintf("%s.gtai", paraloge)
				logrus.Infof("Creating index for paraloge %s of target region %s in: %s ", paraloge, target, paralogeIndexName)
				paralogeIndexPath := filepath.Join(*indexDirParalogePre, paralogeIndexName)
				paralogeIndexFile, err := os.Create(paralogeIndexPath)
				if err != nil {
					logrus.Fatalf("Error creating paraloge index file: %s: %s", paralogeIndexPath, err)
				}
				indexPaths[target] = append(indexPaths[target], paralogeIndexPath)

				// open fa file
				paralogeFastaName := fmt.Sprintf("%s.fa", paraloge)
				paralogeFastaPath := filepath.Join(*fastaDirParalogePre, paralogeFastaName)
				paralogeFastaFile, err := os.Open(paralogeFastaPath)
				if err != nil {
					fmt.Println(err)
					logrus.Fatalf("Error reading paraloge fa file: %s", paralogeFastaPath)
				}

				index.BuildAndSerializeGenomeIndex(paralogeFastaFile, paralogeIndexFile)
			}
		}

		// IV. write paraloges.csv
		metaOut := fmt.Sprintf("%s.csv", *paralogeMetaOutParalogePre)
		extraction.WriteParalogesPre(metaOut, indexPaths)

	} else {
		fmt.Println(parser.Usage("no valid command supplied"))
		return
	}
}
