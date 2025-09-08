package runner

import (
	"fmt"
	"os"

	"github.com/akamensky/argparse"
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
	fmt.Println("exec index")

	fmt.Println(*argsObj.FastaFile)
}
