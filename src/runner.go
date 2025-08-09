package main

import (
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/KleinSamuel/gtamap/src/config"
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
	// ff, errr := os.Create("cpu.pprof")
	// if errr != nil {
	// 	panic(errr)
	// }
	// defer ff.Close()
	// if err := pprof.StartCPUProfile(ff); err != nil {
	// 	panic(err)
	// }
	// defer pprof.StopCPUProfile()
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
	var upstreamIndexPre *int = cmdIndexPre.Int("", "upstream", &argparse.Options{
		Required: false,
		Help:     "Number of bases to add upstream of the gene start position.",
		Default:  0,
	})
	var downstreamIndexPre *int = cmdIndexPre.Int("", "downstream", &argparse.Options{
		Required: false,
		Help:     "Number of bases to add downstream of the gene end position.",
		Default:  0,
	})
	var logLevelIndexPre *string = cmdIndexPre.Selector("", "loglevel", []string{"ERROR", "INFO", "DEBUG"}, &argparse.Options{
		Required: false,
		Help:     "Log output level.",
		Default:  "INFO",
	})

	var cmdIndexPreRegion *argparse.Command = parser.NewCommand("index-pre-region", "Extract a specific sequence from genome.")
	var fastaFileIndexPreRegion *string = cmdIndexPreRegion.String("", "fasta", &argparse.Options{
		Required: true,
		Help:     "Nucleotide sequences (FASTA) file (currently only non-compressed).",
	})
	var outputDirIndexPreRegion *string = cmdIndexPreRegion.String("", "output", &argparse.Options{
		Required: true,
		Help:     "Output directory for extracted gene sequences.",
	})
	var regionIndexPreRegion *string = cmdIndexPreRegion.String("", "region", &argparse.Options{
		Required: true,
		Help:     "Region to extract (e.g. 1:1000-2000).",
	})
	var logLevelIndexPreRegion *string = cmdIndexPreRegion.Selector("", "loglevel", []string{"ERROR", "INFO", "DEBUG"}, &argparse.Options{
		Required: false,
		Help:     "Log output level.",
		Default:  "INFO",
	})

	var cmdIndex *argparse.Command = parser.NewCommand("index", "Build the GTAMap index (.gtai).")
	var fastaFile *os.File = cmdIndex.File("", "fasta", os.O_RDONLY, 0o600, &argparse.Options{
		Required: true,
		Help:     "Nucleotide sequences (FASTA) file.",
	})
	var blacklistFileName *string = cmdIndex.String("", "blacklist", &argparse.Options{
		Required: true,
		Help:     "Genome wide repeat annotation.",
	})
	var outputFileName *string = cmdIndex.String("", "output", &argparse.Options{
		Required: true,
		Help:     "Output file (file extension: .gtai).",
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
	var outputFileMap *string = cmdMap.String("", "output", &argparse.Options{
		Required: true,
		Help:     "Output SAM file (file extension: .sam).",
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
	// var paralogFilePathMap *string = cmdMap.String("", "paralogs", &argparse.Options{
	// 	Required: false,
	// 	Help:     "Paralog region meta file for target regions.",
	// 	Default:  "",
	// })

	// OLD paralog mode (currently not used anymore)
	// var cmdParalogPre *argparse.Command = parser.NewCommand("paralog", "Extract known paralog genes from ENSEMBL Database and prepare paralog.csv for main target index extension.")
	// var geneIdsParalogPre *string = cmdParalogPre.String("", "geneids", &argparse.Options{
	// 	Required: false,
	// 	Help:     "Query Gene IDs for extracting paralog genes from DB (comma-separated).",
	// })
	// var indexDirParalogPre *string = cmdParalogPre.String("", "indexdir", &argparse.Options{
	// 	Required: false,
	// 	Help:     "Target directory for .gai index files of paralog genes.",
	// 	Default:  "index",
	// })
	// var fastaDirParalogPre *string = cmdParalogPre.String("", "fastadir", &argparse.Options{
	// 	Required: false,
	// 	Help:     "Target directory for .fa files of paralog genes.",
	// 	Default:  "fasta_in",
	// })
	// var speciesParalogPre *string = cmdParalogPre.String("", "species", &argparse.Options{
	// 	Required: false,
	// 	Help:     "Corresponding species of target genes. Species name needs to exist in https://www.ensembl.org/index.html",
	// 	Default:  "human",
	// })
	// var fastaFileParalogPre *string = cmdParalogPre.String("", "fasta", &argparse.Options{
	// 	Required: true,
	// 	Help:     "Nucleotide sequences (FASTA) file (currently only non-compressed). Note: A corresponding .fai file needs to be in the same dir as the genome.fa file.",
	// })
	// var gtfFileParalogPre *string = cmdParalogPre.String("", "gtf", &argparse.Options{
	// 	Required: true,
	// 	Help:     "Genome annotation (GTF) file (currently only non-compressed).",
	// })
	// var paralogMetaOutParalogPre *string = cmdParalogPre.String("", "meta", &argparse.Options{
	// 	Required: false,
	// 	Help:     "Prefix of paralog csv file.",
	// 	Default:  "paralogs",
	// })
	// var logLevelParalogPre *string = cmdParalogPre.Selector("", "loglevel", []string{"ERROR", "INFO", "DEBUG"}, &argparse.Options{
	// 	Required: false,
	// 	Help:     "Log output level.",
	// 	Default:  "INFO",
	// })

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
			*outputDirIndexPre, geneIds, *upstreamIndexPre, *downstreamIndexPre, *separateExtraction)

		logrus.Info("Done")

	} else if cmdIndexPreRegion.Happened() {

		level, _ := logrus.ParseLevel(*logLevelIndexPreRegion)
		logrus.SetLevel(level)

		printBanner()
		logrus.Info("Extracting region sequence from genome")

		// parsing region
		region := strings.Split(*regionIndexPreRegion, ":")
		if len(region) != 2 {
			logrus.Fatal("Invalid region format. Expected format: <chromosome>:<start>-<end>")
		}
		chromosome := region[0]
		startEnd := strings.Split(region[1], "-")
		if len(startEnd) != 2 {
			logrus.Fatal("Invalid region format. Expected format: <chromosome>:<start>-<end>")
		}
		start, err := strconv.Atoi(startEnd[0])
		if err != nil || start < 0 {
			logrus.Fatal("Invalid start position. Expected a positive integer.")
		}
		end, err := strconv.Atoi(startEnd[1])
		if err != nil || end <= start {
			logrus.Fatal("Invalid end position. Expected an integer greater than start position.")
		}

		index.ExtractSequenceFromFastaForIndex(*fastaFileIndexPreRegion, chromosome, start, end,
			*outputDirIndexPreRegion)

		logrus.Info("Done")

	} else if cmdIndex.Happened() {

		level, _ := logrus.ParseLevel(*logLevelIndex)
		logrus.SetLevel(level)

		printBanner()
		logrus.Info("Building GTAMap index (.gtai)")

		// ensure suffix .gtai
		if !strings.HasSuffix(*outputFileName, ".gtai") {
			*outputFileName += ".gtai"
		}

		outputFileIndex, err := os.OpenFile(*outputFileName, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0o600)
		if err != nil {
			logrus.Fatalf("Could not open output sam file: %v", err)
		}

		index.BuildAndSerializeGenomeIndex(fastaFile, *blacklistFileName, outputFileIndex)

	} else if cmdMap.Happened() {

		level, _ := logrus.ParseLevel(*logLevelMap)
		logrus.SetLevel(level)

		printBanner()
		logrus.Info("Mapping reads to given index")

		genomeIndex := index.ReadGenomeIndexByFile(indexFile)

		reader, errFastq := fastq.InitFromFiles(fastqFwFile, fastqRwFile)
		if errFastq != nil {
			logrus.Fatalf("Could not initialize fastq reader: %v", errFastq)
		}

		writer := datawriter.InitFromPath(*outputFileMap)

		// Old paralog logic where main index is extended by using paralog indices
		// if *paralogFilePathMap != "" {
		// 	logrus.Info("Found paralog file. Reading in paralog indices.")
		// 	paralogFileMap, err := os.Open(*paralogFilePathMap)
		// 	if err != nil {
		// 		panic("Error reading provided paralog file. Make sure it exists")
		// 	}
		// 	genomeIndex.LoadParalogs(paralogFileMap)
		// }

		mapper.MapAll(genomeIndex, reader, writer, numThreads)

	} else {
		fmt.Println(parser.Usage("no valid command supplied"))
		return
	}
	// else if cmdParalogPre.Happened() {
	// 	level, _ := logrus.ParseLevel(*logLevelParalogPre)
	// 	logrus.SetLevel(level)
	//
	// 	printBanner()
	// 	logrus.Info("Scanning DB for paralogs of specified target IDs.")
	//
	// 	// make output dirs if non existent
	// 	err := os.MkdirAll(*indexDirParalogPre, 0o755)
	// 	if err != nil {
	// 		logrus.Fatalf("Something went wrong creating dir %s: %s", *indexDirParalogPre, err)
	// 	}
	//
	// 	err = os.MkdirAll(*fastaDirParalogPre, 0o755)
	// 	if err != nil {
	// 		logrus.Fatalf("Something went wrong creating dir %s: %s", *indexDirParalogPre, err)
	// 	}
	//
	// 	// parsing gene ids from cmd parser arg
	// 	targetGeneIds := make([]string, 0)
	// 	if *geneIdsParalogPre != "" {
	// 		genes := strings.Split(*geneIdsParalogPre, ",")
	// 		for _, gene := range genes {
	// 			targetGeneIds = append(targetGeneIds, gene)
	// 		}
	// 	}
	//
	// 	// I. get paralog seqs per target gene (a map where each target region has a set of paralog ids)
	// 	targetParalogs := extraction.GetParalogs(targetGeneIds, *speciesParalogPre)
	//
	// 	// II. check if some of the gene ids are already present in --fastadir
	// 	paralogsToExtractSeq := index.OptimizeFastaExtraction(targetParalogs, fastaDirParalogPre)
	//
	// 	// III. extract all identified non-existent paralogs into separate .fa files in '--fastadir'
	// 	for target, paralogs := range paralogsToExtractSeq {
	// 		if len(paralogs) > 0 {
	// 			logrus.Infof("Extracting %s paralog sequences of target region %s into %s", strconv.Itoa(len(paralogs)), target, *fastaDirParalogPre)
	// 			index.ExtractGeneSequenceFromGtfAndFastaForIndex(*gtfFileParalogPre, *fastaFileParalogPre,
	// 				*fastaDirParalogPre, paralogs, 0, 0, true)
	// 			logrus.Infof("Finished extracting %s paralog sequences of target region %s into %s", strconv.Itoa(len(paralogs)), target, *fastaDirParalogPre)
	// 		}
	// 	}
	//
	// 	// IV. get only ids where index non-existent in --indexdir
	// 	paralogsToSerialize := index.OptimizeIndexSerialisation(targetParalogs, indexDirParalogPre)
	//
	// 	// V. read in all seqs in '--fastadir' and serialize index into '--indexdir' for each paralog seq
	// 	index.BuildAndSerializeAll(paralogsToSerialize, indexDirParalogPre, fastaDirParalogPre)
	//
	// 	// VI. gather abs paths to indices into map[target][]paralogIndexPaths
	// 	indexPaths := extraction.GetAbsPathsPerTarget(targetParalogs, indexDirParalogPre)
	//
	// 	// VI. write paralogs.csv
	// 	metaOut := fmt.Sprintf("%s.csv", *paralogMetaOutParalogPre)
	// 	extraction.WriteParalogsPre(metaOut, indexPaths)
	// 	logrus.Infof("Wrote target to paralog mapping to destination file: %s", metaOut)
}
