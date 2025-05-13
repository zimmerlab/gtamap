package main

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/extraction"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper"
	"github.com/KleinSamuel/gtamap/src/datawriter"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/akamensky/argparse"
	"github.com/sirupsen/logrus"
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
	var fastaFile *os.File = cmdIndex.File("", "fasta", os.O_RDONLY, 0o600, &argparse.Options{
		Required: true,
		Help:     "Nucleotide sequences (FASTA) file.",
	})
	var outputFileIndex *os.File = cmdIndex.File("", "output", os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0o600, &argparse.Options{
		Required: true,
		Help:     "Output file (.gtai).",
	})
	var logLevelIndex *string = cmdIndex.Selector("", "loglevel", []string{"ERROR", "INFO", "DEBUG"}, &argparse.Options{
		Required: false,
		Help:     "Log output level.",
		Default:  "INFO",
	})

	var cmdMap *argparse.Command = parser.NewCommand("map", "Map reads to the GTAMap index.")
	var indexFile *os.File = cmdMap.File("", "index", os.O_RDONLY, 0o600, &argparse.Options{
		Required: true,
		Help:     "GTAMap index (.gtai) file",
	})
	var fastqFwFile *os.File = cmdMap.File("", "r1", os.O_RDONLY, 0o600, &argparse.Options{
		Required: true,
		Help:     "FASTQ file containing the forward reads.",
	})
	var fastqRwFile *os.File = cmdMap.File("", "r2", os.O_RDONLY, 0o600, &argparse.Options{
		Required: false,
		Help:     "FASTQ file containing the reverse reads.",
	})
	var outputFileMap *os.File = cmdMap.File("", "output", os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0o600, &argparse.Options{
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
	var paralogFilePathMap *string = cmdMap.String("", "paralogs", &argparse.Options{
		Required: false,
		Help:     "Paralog region meta file for target regions.",
		Default:  "",
	})

	// paralog mode
	var cmdParalogPre *argparse.Command = parser.NewCommand("paralog", "Extract known paralog genes from ENSEMBL Database and prepare paralog.csv for main target index extension.")
	var geneIdsParalogPre *string = cmdParalogPre.String("", "geneids", &argparse.Options{
		Required: false,
		Help:     "Query Gene IDs for extracting paralog genes from DB (comma-separated).",
	})
	var indexDirParalogPre *string = cmdParalogPre.String("", "indexdir", &argparse.Options{
		Required: false,
		Help:     "Target directory for .gai index files of paralog genes.",
		Default:  "index",
	})
	var fastaDirParalogPre *string = cmdParalogPre.String("", "fastadir", &argparse.Options{
		Required: false,
		Help:     "Target directory for .fa files of paralog genes.",
		Default:  "fasta_in",
	})
	var speciesParalogPre *string = cmdParalogPre.String("", "species", &argparse.Options{
		Required: false,
		Help:     "Corresponding species of target genes. Species name needs to exist in https://www.ensembl.org/index.html",
		Default:  "human",
	})
	var fastaFileParalogPre *string = cmdParalogPre.String("", "fasta", &argparse.Options{
		Required: true,
		Help:     "Nucleotide sequences (FASTA) file (currently only non-compressed). Note: A corresponding .fai file needs to be in the same dir as the genome.fa file.",
	})
	var gtfFileParalogPre *string = cmdParalogPre.String("", "gtf", &argparse.Options{
		Required: true,
		Help:     "Genome annotation (GTF) file (currently only non-compressed).",
	})
	var paralogMetaOutParalogPre *string = cmdParalogPre.String("", "meta", &argparse.Options{
		Required: false,
		Help:     "Prefix of paralog csv file.",
		Default:  "paralogs",
	})
	var logLevelParalogPre *string = cmdParalogPre.Selector("", "loglevel", []string{"ERROR", "INFO", "DEBUG"}, &argparse.Options{
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

		reader := fastq.InitFromFiles(fastqFwFile, fastqRwFile)

		writer := datawriter.InitFromFile(outputFileMap)

		if *paralogFilePathMap != "" {
			paralogFileMap, err := os.Open(*paralogFilePathMap)
			if err != nil {
				panic("Error reading provided paralog file. Make sure it exists")
			}
			genomeIndex.LoadParalogs(paralogFileMap)
		}

		mapper.MapAll(genomeIndex, reader, writer, numThreads)

	} else if cmdParalogPre.Happened() {
		level, _ := logrus.ParseLevel(*logLevelParalogPre)
		logrus.SetLevel(level)

		printBanner()
		logrus.Info("Scanning DB for paralogs of specified target IDs.")

		// parsing gene ids
		targetGeneIds := make([]string, 0)
		if *geneIdsParalogPre != "" {
			genes := strings.Split(*geneIdsParalogPre, ",")
			for _, gene := range genes {
				targetGeneIds = append(targetGeneIds, gene)
			}
		}
		// I. get paralog seqs per target gene
		targetParalogs := extraction.GetParaloges(targetGeneIds, *speciesParalogPre)

		// II. extract all target genes into separate .fa files in '--fastadir'
		for target, paralogs := range targetParalogs {
			logrus.Infof("Extracting sequences for paralogs of %s", target)
			index.ExtractGeneSequenceFromGtfAndFastaForIndex(*gtfFileParalogPre, *fastaFileParalogPre,
				*fastaDirParalogPre, paralogs, true)
		}

		// III. read in all seqs in '--fastadir' and serialize index into '--indexdir' for each paralog seq
		indexPaths := make(map[string][]string)
		for target, paralogs := range targetParalogs {
			for paralog := range paralogs {
				err := os.MkdirAll(*indexDirParalogPre, 0o755)
				if err != nil {
					logrus.Fatalf("Something went wrong creating dir %s: %s", *indexDirParalogPre, err)
				}

				// create and format index file
				paralogIndexName := fmt.Sprintf("%s.gtai", paralog)
				logrus.Infof("Creating index for paralog %s of target region %s in: %s ", paralog, target, paralogIndexName)
				paralogIndexPath := filepath.Join(*indexDirParalogPre, paralogIndexName)
				paralogIndexFile, err := os.Create(paralogIndexPath)
				if err != nil {
					logrus.Fatalf("Error creating paralog index file: %s: %s", paralogIndexPath, err)
				}
				indexPaths[target] = append(indexPaths[target], paralogIndexPath)

				// open fa file
				paralogFastaName := fmt.Sprintf("%s.fa", paralog)
				paralogFastaPath := filepath.Join(*fastaDirParalogPre, paralogFastaName)
				paralogFastaFile, err := os.Open(paralogFastaPath)
				if err != nil {
					fmt.Println(err)
					logrus.Fatalf("Error reading paralog fa file: %s", paralogFastaPath)
				}

				index.BuildAndSerializeGenomeIndex(paralogFastaFile, paralogIndexFile)
			}
		}

		// IV. write paralogs.csv
		metaOut := fmt.Sprintf("%s.csv", *paralogMetaOutParalogPre)
		extraction.WriteParalogsPre(metaOut, indexPaths)

	} else {
		fmt.Println(parser.Usage("no valid command supplied"))
		return
	}
}
