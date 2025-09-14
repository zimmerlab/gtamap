package runner

import (
	"fmt"
	"strconv"
	"strings"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/sirupsen/logrus"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

func GetCommandIndexPreRegion() *cobra.Command {

	var fastaFilePath string
	var fastaIndexFilePath string
	var singleFile bool
	var outputDirPath string
	var fastaFileName string
	var regions []string

	indexPreCmd := &cobra.Command{
		Use:   "index-pre-region",
		Short: "Extract target region sequences from a genome",
		Run: func(cmd *cobra.Command, args []string) {

			fmt.Println("output dir from config:", config.Mapper.General.OutputDir)

			config.Mapper.SetIndexFastaIndex(config.Mapper.Index.FastaIndexFilePath)

			ExecIndexPreRegion()
		},
	}

	flags := indexPreCmd.Flags()

	flags.StringVarP(
		&fastaFilePath,
		"fasta",
		"f",
		"",
		"Fasta file (required)",
	)
	indexPreCmd.MarkFlagRequired("fasta")
	viper.BindPFlag(
		"index.fasta_file_path",
		indexPreCmd.Flags().Lookup("fasta"),
	)

	flags.StringVarP(
		&fastaIndexFilePath,
		"fasta-index",
		"i",
		"",
		"Fasta index file (default: [--fasta].fai)",
	)
	viper.BindPFlag(
		"index.fasta_index_file_path",
		indexPreCmd.Flags().Lookup("fasta-index"),
	)

	flags.StringVarP(
		&outputDirPath,
		"output",
		"o",
		"",
		"Output directory",
	)
	viper.BindPFlag(
		"general.output_dir",
		indexPreCmd.Flags().Lookup("output"),
	)

	flags.BoolVarP(
		&singleFile,
		"single-file",
		"s",
		false,
		"Write all gene sequences to a single fasta file",
	)
	viper.BindPFlag(
		"index.output.single_file",
		indexPreCmd.Flags().Lookup("single-file"),
	)

	flags.StringVar(
		&fastaFileName,
		"fasta-file-name",
		"",
		"Output FASTA file name (within output directory) (only for "+
			"--single-file) (default: genes.fa)",
	)
	viper.BindPFlag(
		"index.output.fasta_file_name",
		indexPreCmd.Flags().Lookup("fasta-file-name"),
	)

	flags.StringSliceVarP(
		&regions,
		"regions",
		"r",
		nil,
		"List of regions to extract (e.g. 1:1000-2000,2:1500-2500) (required)",
	)
	indexPreCmd.MarkFlagRequired("regions")
	viper.BindPFlag(
		"index.regions",
		indexPreCmd.Flags().Lookup("regions"),
	)

	return indexPreCmd
}

type IndexPreRegionArgs struct {
	Contig string
	Start  int
	End    int
}

func parseRegionString(regionStr string) (IndexPreRegionArgs, error) {

	rArr := strings.Split(regionStr, ":")
	if len(rArr) != 2 {
		return IndexPreRegionArgs{}, fmt.Errorf("invalid region format: %s", regionStr)
	}

	contig := rArr[0]

	r2Arr := strings.Split(rArr[1], "-")
	if len(r2Arr) != 2 {
		return IndexPreRegionArgs{}, fmt.Errorf("invalid region format: %s", regionStr)
	}

	start, errStart := strconv.Atoi(r2Arr[0])
	if errStart != nil || start < 0 {
		return IndexPreRegionArgs{}, fmt.Errorf("invalid start position in region: %s", regionStr)
	}

	end, errEnd := strconv.Atoi(r2Arr[1])
	if errEnd != nil || end <= start {
		return IndexPreRegionArgs{}, fmt.Errorf("invalid end position in region: %s", regionStr)
	}

	return IndexPreRegionArgs{
		Contig: contig,
		Start:  start,
		End:    end,
	}, nil
}

func ExecIndexPreRegion() {

	logrus.Info("Extracting region sequence from genome")

	regions := make([]IndexPreRegionArgs, len(config.Mapper.Index.Regions))
	for i, regionStr := range config.Mapper.Index.Regions {
		r, err := parseRegionString(regionStr)
		if err != nil {
			logrus.Fatalf("Error parsing region '%s': %v", regionStr, err)
		}
		regions[i] = r
	}

	index.ExtractSequenceFromFastaForIndex(
		config.Mapper.Index.FastaFilePath,
		config.Mapper.Index.FastaIndexFilePath,
		regions[0].Contig,
		regions[0].Start,
		regions[0].End,
		config.Mapper.General.OutputDir,
	)
}
