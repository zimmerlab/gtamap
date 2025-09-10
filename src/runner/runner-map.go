package runner

import (
	"fmt"
	"log"
	"os"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper"
	"github.com/KleinSamuel/gtamap/src/datawriter"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/akamensky/argparse"
	"github.com/sirupsen/logrus"
)

type ArgsMap struct {
	IndexFile              *os.File
	FastqR1File            *os.File
	FastqR2File            *os.File
	OutputFile             *string
	ReadType               *string
	Threads                *int
	RegionmaskBedFile      *os.File
	RegionmaskPriorityFile *os.File
	RegionmaskAction       *string
}

func AddCommandMap(
	parser *argparse.Parser,
) (
	*argparse.Command,
	*ArgsMap,
) {

	var command *argparse.Command = parser.NewCommand(
		"map",
		"Map reads to the GTAMap index",
	)

	argsObj := &ArgsMap{}

	argsObj.IndexFile = command.File(
		"",
		"index",
		os.O_RDONLY,
		0o600,
		&argparse.Options{
			Required: true,
			Help:     "GTAMap index (.gtai) file",
		},
	)

	argsObj.FastqR1File = command.File(
		"",
		"r1",
		os.O_RDONLY,
		0o600,
		&argparse.Options{
			Required: true,
			Help:     "FASTQ file containing the forward reads",
		},
	)

	argsObj.FastqR2File = command.File(
		"",
		"r2",
		os.O_RDONLY,
		0o600,
		&argparse.Options{
			Required: false,
			Help:     "FASTQ file containing the reverse reads",
		},
	)

	argsObj.OutputFile = command.String(
		"",
		"output",
		&argparse.Options{
			Required: true,
			Help:     "Output SAM file (file extension: .sam)",
			Default:  "",
		},
	)

	// var logFileMap *string = cmdMap.String("", "log", &argparse.Options{
	// 	Help:    "Output log file.",
	// 	Default: "",
	// })

	argsObj.ReadType = command.Selector(
		"",
		"read-type",
		[]string{"DNA", "RNA"},
		&argparse.Options{
			Required: true,
			Help:     "Specify read type",
		},
	)

	argsObj.Threads = command.Int(
		"",
		"threads",
		&argparse.Options{
			Required: false,
			Help:     "Number of threads to use (default: all)",
			Default:  -1,
		},
	)

	argsObj.RegionmaskBedFile = command.File(
		"",
		"regionmaskBed",
		os.O_RDONLY,
		0o600,
		&argparse.Options{
			Required: false,
			Help: "Path to BED file containing regions to mask during " +
				"index construction",
		},
	)

	argsObj.RegionmaskPriorityFile = command.File(
		"",
		"regionmaskPriorities",
		os.O_RDONLY,
		0o600,
		&argparse.Options{
			Required: false,
			Help: "Path to TSV file containing region names and their " +
				"corresponding priority (higher number = higher priority)",
		},
	)

	argsObj.RegionmaskAction = command.Selector(
		"",
		"regionmaskAction",
		[]string{"combine", "ignore-index"},
		&argparse.Options{
			Required: false,
			Help:     "Specify how to handle regionmasks",
			Default:  "combine",
		},
	)

	return command, argsObj
}

func ExecMap(argsObj *ArgsMap) {

	// if *logFileMap != "" {
	// 	config.LogOut = *logFileMap
	// }

	// should already be captured by argument parsers
	if *argsObj.ReadType != "DNA" && *argsObj.ReadType != "RNA" {
		log.Fatalf("Invalid read type: %s. Must be DNA or RNA.", *argsObj.ReadType)
	}

	config.IsOriginRNA = *argsObj.ReadType == "RNA"

	genomeIndex := index.ReadGenomeIndexByFile(argsObj.IndexFile)

	if *argsObj.RegionmaskBedFile == (os.File{}) {
		argsObj.RegionmaskBedFile = nil
	}
	if *argsObj.RegionmaskPriorityFile == (os.File{}) {
		argsObj.RegionmaskPriorityFile = nil
	}
	if (argsObj.RegionmaskBedFile == nil &&
		argsObj.RegionmaskPriorityFile != nil) ||
		(argsObj.RegionmaskBedFile != nil &&
			argsObj.RegionmaskPriorityFile == nil) {
		logrus.Fatal("Both --regionmaskBed and --regionmaskPriorities " +
			"need to be provided to use region masking.")
	}

	// add region mask to index if provided
	if argsObj.RegionmaskBedFile != nil && argsObj.RegionmaskPriorityFile != nil {

		fmt.Println("regionmask action:", *argsObj.RegionmaskAction)

		index.AddRegionmaskToIndex(
			argsObj.RegionmaskBedFile,
			argsObj.RegionmaskPriorityFile,
			genomeIndex,
		)

		logrus.WithFields(logrus.Fields{
			"region mask":   argsObj.RegionmaskBedFile.Name(),
			"priority file": argsObj.RegionmaskPriorityFile.Name(),
		}).Info("Using region mask bed file and priority file")

	}

	reader, errFastq := fastq.InitFromFiles(
		argsObj.FastqR1File,
		argsObj.FastqR2File,
	)

	if errFastq != nil {
		logrus.Fatalf("Could not initialize fastq reader: %v", errFastq)
	}

	writer := datawriter.InitFromPath(*argsObj.OutputFile)

	// Old paralog logic where main index is extended by using paralog indices
	// if *paralogFilePathMap != "" {
	// 	logrus.Info("Found paralog file. Reading in paralog indices.")
	// 	paralogFileMap, err := os.Open(*paralogFilePathMap)
	// 	if err != nil {
	// 		panic("Error reading provided paralog file. Make sure it exists")
	// 	}
	// 	genomeIndex.LoadParalogs(paralogFileMap)
	// }

	// TODO: maybe check threads parameter for valid values

	mapper.MapAll(genomeIndex, reader, writer, argsObj.Threads)
}
