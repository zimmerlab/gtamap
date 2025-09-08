package runner

import (
	"strconv"
	"strings"

	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/akamensky/argparse"
	"github.com/sirupsen/logrus"
)

type ArgsIndexPreRegion struct {
	FastaFile *string
	OutputDir *string
	Region    *string
}

func AddCommandIndexPreRegion(
	parser *argparse.Parser,
) (
	*argparse.Command,
	*ArgsIndexPreRegion,
) {

	var command *argparse.Command = parser.NewCommand(
		"index-pre-region",
		"Extract a specific sequence from genome.",
	)

	argsObj := &ArgsIndexPreRegion{}

	argsObj.FastaFile = command.String(
		"",
		"fasta",
		&argparse.Options{
			Required: true,
			Help: "Nucleotide sequences (FASTA) file (currently only " +
				"non-compressed)",
		},
	)

	argsObj.OutputDir = command.String(
		"",
		"output",
		&argparse.Options{
			Required: true,
			Help:     "Output directory for extracted gene sequences",
		},
	)

	argsObj.Region = command.String(
		"",
		"region",
		&argparse.Options{
			Required: true,
			Help:     "Region to extract (e.g. 1:1000-2000)",
		},
	)

	return command, argsObj
}

func ExecIndexPreRegion(argsObj *ArgsIndexPreRegion) {

	logrus.Info("Extracting region sequence from genome")

	// parse region string to contig, start, end
	region := strings.Split(*argsObj.Region, ":")
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

	index.ExtractSequenceFromFastaForIndex(
		*argsObj.FastaFile,
		chromosome,
		start,
		end,
		*argsObj.OutputDir,
	)
}
