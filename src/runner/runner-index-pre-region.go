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

// type ArgsIndexPreRegion struct {
// 	FastaFile *string
// 	OutputDir *string
// 	Region    *string
// }
//
// func AddCommandIndexPreRegion(
// 	parser *argparse.Parser,
// ) (
// 	*argparse.Command,
// 	*ArgsIndexPreRegion,
// ) {
//
// 	var command *argparse.Command = parser.NewCommand(
// 		"index-pre-region",
// 		"Extract a specific sequence from genome.",
// 	)
//
// 	argsObj := &ArgsIndexPreRegion{}
//
// 	argsObj.FastaFile = command.String(
// 		"",
// 		"fasta",
// 		&argparse.Options{
// 			Required: true,
// 			Help: "Nucleotide sequences (FASTA) file (currently only " +
// 				"non-compressed)",
// 		},
// 	)
//
// 	argsObj.OutputDir = command.String(
// 		"",
// 		"output",
// 		&argparse.Options{
// 			Required: true,
// 			Help:     "Output directory for extracted gene sequences",
// 		},
// 	)
//
// 	argsObj.Region = command.String(
// 		"",
// 		"region",
// 		&argparse.Options{
// 			Required: true,
// 			Help:     "Region to extract (e.g. 1:1000-2000)",
// 		},
// 	)
//
// 	return command, argsObj
// }
//
// func ExecIndexPreRegion(argsObj *ArgsIndexPreRegion) {
//
// 	logrus.Info("Extracting region sequence from genome")
//
// 	// parse region string to contig, start, end
// 	region := strings.Split(*argsObj.Region, ":")
// 	if len(region) != 2 {
// 		logrus.Fatal("Invalid region format. Expected format: <chromosome>:<start>-<end>")
// 	}
//
// 	chromosome := region[0]
//
// 	startEnd := strings.Split(region[1], "-")
// 	if len(startEnd) != 2 {
// 		logrus.Fatal("Invalid region format. Expected format: <chromosome>:<start>-<end>")
// 	}
//
// 	start, err := strconv.Atoi(startEnd[0])
// 	if err != nil || start < 0 {
// 		logrus.Fatal("Invalid start position. Expected a positive integer.")
// 	}
//
// 	end, err := strconv.Atoi(startEnd[1])
// 	if err != nil || end <= start {
// 		logrus.Fatal("Invalid end position. Expected an integer greater than start position.")
// 	}
//
// 	index.ExtractSequenceFromFastaForIndex(
// 		*argsObj.FastaFile,
// 		chromosome,
// 		start,
// 		end,
// 		*argsObj.OutputDir,
// 	)
// }

func GetCommandIndexPreRegion(v *viper.Viper) *cobra.Command {

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
	v.BindPFlag(
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
	v.BindPFlag(
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
	v.BindPFlag(
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
	v.BindPFlag(
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
	v.BindPFlag(
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
	v.BindPFlag(
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

	fmt.Println(config.Mapper.Index.Regions)

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

//
// 	logrus.Info("Extracting gene sequences from genome")
//
// 	// convert gene id list to set
// 	geneIds := make(map[string]struct{})
// 	for _, gene := range config.Mapper.Index.GeneIds {
// 		geneIds[gene] = struct{}{}
// 	}
//
// 	index.ExtractGeneSequenceFromGtfAndFastaForIndex(
// 		config.Mapper.Index.GtfFilePath,
// 		config.Mapper.Index.FastaFilePath,
// 		config.Mapper.Index.FastaIndexFilePath,
// 		config.Mapper.General.OutputDir,
// 		geneIds,
// 		config.Mapper.Index.UpstreamBases,
// 		config.Mapper.Index.DownstreamBases,
// 		config.Mapper.Index.Output.SingleFile,
// 	)
// }
