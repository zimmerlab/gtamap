package runner

import (
	"fmt"
	"os"

	"github.com/akamensky/argparse"
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

	return command, argsObj
}

func ExecMap(argsObj *ArgsMap) {
	fmt.Println("exec map")

	fmt.Println(*argsObj.IndexFile)
}
