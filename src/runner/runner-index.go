package runner

import (
	"os"
	"strings"

	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/akamensky/argparse"
	"github.com/sirupsen/logrus"
)

type ArgsIndex struct {
	FastaFile              *os.File
	RegionmaskBedFile      *os.File
	RegionmaskPriorityFile *os.File
	OutputFile             *string
}

func AddCommandIndex(
	parser *argparse.Parser,
) (
	*argparse.Command,
	*ArgsIndex,
) {

	var command *argparse.Command = parser.NewCommand(
		"index",
		"Build the GTAMap index (.gtai)",
	)

	argsObj := &ArgsIndex{}

	argsObj.FastaFile = command.File(
		"",
		"fasta",
		os.O_RDONLY,
		0o600,
		&argparse.Options{
			Required: true,
			Help:     "Nucleotide sequences (FASTA) file",
		},
	)

	// repeatMaskerFile := cmdIndex.File(
	// 	"", "repeatmask",
	// 	os.O_RDONLY, 0o600,
	// 	&argparse.Options{
	// 		Required: false,
	// 		Help:     "Path to genome wide repeatmasker file (hg38.fa.out).",
	// 	})

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

	argsObj.OutputFile = command.String(
		"",
		"output",
		&argparse.Options{
			Required: true,
			Help:     "Output file (file extension: .gtai)",
		},
	)

	return command, argsObj
}

func ExecIndex(argsObj *ArgsIndex) {

	// ensure .gtai suffix of output file
	if !strings.HasSuffix(*argsObj.OutputFile, ".gtai") {
		*argsObj.OutputFile += ".gtai"
	}

	outputFileIndex, err := os.OpenFile(
		*argsObj.OutputFile,
		os.O_WRONLY|os.O_CREATE|os.O_TRUNC,
		0o600,
	)

	if err != nil {
		logrus.Fatalf("Could not open output sam file: %v", err)
	}

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

	index.BuildAndSerializeGenomeIndex(
		argsObj.FastaFile,
		outputFileIndex,
		argsObj.RegionmaskBedFile,
		argsObj.RegionmaskPriorityFile,
	)
}
