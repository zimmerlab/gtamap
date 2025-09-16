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
	var regionmaskFilePath string
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
			SetConfigValue("index.regionmask_file_path", regionmaskFilePath)

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

	flags.StringVarP(
		&regionmaskFilePath,
		"regionmask",
		"m",
		"",
		"Regionmask file (.bed) containing specific mismatch constraints per region",
	)

	return indexCmd
}

func ExecIndex() {

	fastaFile, err := os.Open(config.Mapper.Index.FastaFilePath)
	if err != nil {
		logrus.Fatalf("Could not open fasta file: %v", err)
	}

	outputFile, err := config.Mapper.GetIndexOutputFile()
	if err != nil {
		logrus.Fatalf("Could not open output index file: %v", err)
	}

	var regionmaskFile *os.File = nil
	var errRegionMask error = nil

	if config.Mapper.Index.RegionmaskFilePath != "" {
		regionmaskFile, errRegionMask = os.Open(config.Mapper.Index.RegionmaskFilePath)
		if errRegionMask != nil {
			logrus.Fatalf("Could not open regionmask file: %v", errRegionMask)
		}
		defer regionmaskFile.Close()
	}

	index.BuildAndSerializeGenomeIndex(
		fastaFile,
		outputFile,
		regionmaskFile,
	)
}
