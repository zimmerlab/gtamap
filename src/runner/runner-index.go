package runner

import (
	"os"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/sirupsen/logrus"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

func GetCommandIndex() *cobra.Command {

	var fastaFilePath string
	// var regionmaskBedFilePath string
	// var regionmaskPriorityFilePath string
	var outputDirPath string
	var indexFileName string
	var useFastaFileName bool

	indexCmd := &cobra.Command{
		Use:   "index",
		Short: "Build the gtamap index (.gtai)",
		Run: func(cmd *cobra.Command, args []string) {

			if indexFileName != "" && useFastaFileName {
				cmd.PrintErrf("\n%s [-i | --index-file-name] and [-u | --use-fasta-file-name] cannot be used together.\n\n", cmd.ErrPrefix())
				cmd.Usage()
				os.Exit(1)
			}

			SetConfigValue("index.fasta_file_path", fastaFilePath)
			SetConfigValue("general.output_dir", outputDirPath)
			SetConfigValue("index.output.index_file_name", indexFileName)
			SetConfigValue("index.output.use_fasta_file_name", useFastaFileName)

			if err := viper.Unmarshal(config.Mapper); err != nil {
				logrus.Fatalf("Unable to decode config: %v", err)
			}

			ExecIndex()
		},
	}

	flags := indexCmd.Flags()

	flags.StringVarP(
		&fastaFilePath,
		"fasta",
		"f",
		"",
		"Fasta file (required)",
	)
	indexCmd.MarkFlagRequired("fasta")

	flags.StringVarP(
		&outputDirPath,
		"output",
		"o",
		"",
		"Output directory (required)",
	)
	indexCmd.MarkFlagRequired("output")

	flags.StringVarP(
		&indexFileName,
		"index-file-name",
		"i",
		"",
		"Output gtamap index file name (within output directory) (default: index.gtai)",
	)

	flags.BoolVarP(
		&useFastaFileName,
		"use-fasta-file-name",
		"u",
		false,
		"Use the name of the fasta file (without extension) as index file name (within output directory)",
	)

	return indexCmd
}

func ExecIndex() {

	// // ensure .gtai suffix of output file
	// if !strings.HasSuffix(*argsObj.OutputFile, ".gtai") {
	// 	*argsObj.OutputFile += ".gtai"
	// }
	//
	// outputFileIndex, err := os.OpenFile(
	// 	*argsObj.OutputFile,
	// 	os.O_WRONLY|os.O_CREATE|os.O_TRUNC,
	// 	0o600,
	// )
	//
	// if err != nil {
	// 	logrus.Fatalf("Could not open output sam file: %v", err)
	// }
	//
	// if *argsObj.RegionmaskBedFile == (os.File{}) {
	// 	argsObj.RegionmaskBedFile = nil
	// }
	// if *argsObj.RegionmaskPriorityFile == (os.File{}) {
	// 	argsObj.RegionmaskPriorityFile = nil
	// }
	// if (argsObj.RegionmaskBedFile == nil &&
	// 	argsObj.RegionmaskPriorityFile != nil) ||
	// 	(argsObj.RegionmaskBedFile != nil &&
	// 		argsObj.RegionmaskPriorityFile == nil) {
	// 	logrus.Fatal("Both --regionmaskBed and --regionmaskPriorities " +
	// 		"need to be provided to use region masking.")
	// }

	fastaFile, err := os.Open(config.Mapper.Index.FastaFilePath)
	if err != nil {
		logrus.Fatalf("Could not open fasta file: %v", err)
	}

	outputFile, err := config.Mapper.GetIndexOutputFile()
	if err != nil {
		logrus.Fatalf("Could not open output index file: %v", err)
	}

	index.BuildAndSerializeGenomeIndex(
		fastaFile,
		outputFile,
		nil,
		nil,
	)
}

// type ArgsIndex struct {
// 	FastaFile              *os.File
// 	RegionmaskBedFile      *os.File
// 	RegionmaskPriorityFile *os.File
// 	OutputFile             *string
// }
//
// func AddCommandIndex(
// 	parser *argparse.Parser,
// ) (
// 	*argparse.Command,
// 	*ArgsIndex,
// ) {
//
// 	var command *argparse.Command = parser.NewCommand(
// 		"index",
// 		"Build the GTAMap index (.gtai)",
// 	)
//
// 	argsObj := &ArgsIndex{}
//
// 	argsObj.FastaFile = command.File(
// 		"",
// 		"fasta",
// 		os.O_RDONLY,
// 		0o600,
// 		&argparse.Options{
// 			Required: true,
// 			Help:     "Nucleotide sequences (FASTA) file",
// 		},
// 	)
//
// 	// repeatMaskerFile := cmdIndex.File(
// 	// 	"", "repeatmask",
// 	// 	os.O_RDONLY, 0o600,
// 	// 	&argparse.Options{
// 	// 		Required: false,
// 	// 		Help:     "Path to genome wide repeatmasker file (hg38.fa.out).",
// 	// 	})
//
// 	argsObj.RegionmaskBedFile = command.File(
// 		"",
// 		"regionmaskBed",
// 		os.O_RDONLY,
// 		0o600,
// 		&argparse.Options{
// 			Required: false,
// 			Help: "Path to BED file containing regions to mask during " +
// 				"index construction",
// 		},
// 	)
//
// 	argsObj.RegionmaskPriorityFile = command.File(
// 		"",
// 		"regionmaskPriorities",
// 		os.O_RDONLY,
// 		0o600,
// 		&argparse.Options{
// 			Required: false,
// 			Help: "Path to TSV file containing region names and their " +
// 				"corresponding priority (higher number = higher priority)",
// 		},
// 	)
//
// 	argsObj.OutputFile = command.String(
// 		"",
// 		"output",
// 		&argparse.Options{
// 			Required: true,
// 			Help:     "Output file (file extension: .gtai)",
// 		},
// 	)
//
// 	return command, argsObj
// }
//
// func ExecIndex(argsObj *ArgsIndex) {
//
// 	// ensure .gtai suffix of output file
// 	if !strings.HasSuffix(*argsObj.OutputFile, ".gtai") {
// 		*argsObj.OutputFile += ".gtai"
// 	}
//
// 	outputFileIndex, err := os.OpenFile(
// 		*argsObj.OutputFile,
// 		os.O_WRONLY|os.O_CREATE|os.O_TRUNC,
// 		0o600,
// 	)
//
// 	if err != nil {
// 		logrus.Fatalf("Could not open output sam file: %v", err)
// 	}
//
// 	if *argsObj.RegionmaskBedFile == (os.File{}) {
// 		argsObj.RegionmaskBedFile = nil
// 	}
// 	if *argsObj.RegionmaskPriorityFile == (os.File{}) {
// 		argsObj.RegionmaskPriorityFile = nil
// 	}
// 	if (argsObj.RegionmaskBedFile == nil &&
// 		argsObj.RegionmaskPriorityFile != nil) ||
// 		(argsObj.RegionmaskBedFile != nil &&
// 			argsObj.RegionmaskPriorityFile == nil) {
// 		logrus.Fatal("Both --regionmaskBed and --regionmaskPriorities " +
// 			"need to be provided to use region masking.")
// 	}
//
// 	index.BuildAndSerializeGenomeIndex(
// 		argsObj.FastaFile,
// 		outputFileIndex,
// 		argsObj.RegionmaskBedFile,
// 		argsObj.RegionmaskPriorityFile,
// 	)
// }
