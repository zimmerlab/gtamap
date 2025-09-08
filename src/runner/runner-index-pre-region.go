package runner

import (
	"fmt"

	"github.com/akamensky/argparse"
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
	fmt.Println("exec index-pre-region")

	fmt.Println(*argsObj.FastaFile)
}
