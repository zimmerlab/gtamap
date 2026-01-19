package runner

import (
	"os"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper"
	"github.com/KleinSamuel/gtamap/src/datawriter"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/sirupsen/logrus"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

func GetCommandMap() *cobra.Command {
	var indexFilePath string
	var fastqR1FilePath string
	var fastqR2FilePath string
	var readOrigin string
	var regionmaskFilePath string

	var outputDirPath string
	var samFileName string
	var threads int

	var statsFile string
	var progressFile string
	var progressStageTimeFile string

	mapCmd := &cobra.Command{
		Use:   "map",
		Short: "Run mapping",
		Run: func(cmd *cobra.Command, args []string) {
			SetConfigValue("mapping.index_file_path", indexFilePath)
			SetConfigValue("mapping.fastq_r1_file_path", fastqR1FilePath)
			SetConfigValue("mapping.fastq_r2_file_path", fastqR2FilePath)
			SetConfigValue("mapping.read_origin", readOrigin)
			SetConfigValue("general.output_dir", outputDirPath)
			SetConfigValue("mapping.output.sam_file_name", samFileName)
			SetConfigValue("mapping.threads", threads)
			SetConfigValue("mapping.regionmask_file_path", regionmaskFilePath)
			SetConfigValue("mapping.output.mapping_stats_file_name", statsFile)
			SetConfigValue("mapping.output.mapping_progress_file_name", progressFile)
			SetConfigValue("mapping.output.mapping_progress_stage_time_file_name", progressStageTimeFile)

			if err := viper.Unmarshal(config.Mapper); err != nil {
				logrus.Fatalf("Unable to decode config: %v", err)
			}

			errOrigin := config.Mapper.SetReadOrigin(config.Mapper.Mapping.ReadOrigin)
			if errOrigin != nil {
				cmd.PrintErrf("\n%s %s [-r | --read-origin]\n\n", cmd.ErrPrefix(), errOrigin.Error())
				cmd.Usage()
				os.Exit(1)
			}

			config.Mapper.SetMappingThreads(config.Mapper.Mapping.Threads)

			ExecMap()
		},
	}

	flags := mapCmd.Flags()

	flags.StringVarP(
		&indexFilePath,
		"index",
		"i",
		"",
		"Index file (*.gtai) (required)",
	)
	mapCmd.MarkFlagRequired("index")

	flags.StringVarP(
		&fastqR1FilePath,
		"reads-r1",
		"1",
		"",
		"FASTQ file containing the R1 reads (required)",
	)
	mapCmd.MarkFlagRequired("reads-r1")

	flags.StringVarP(
		&fastqR2FilePath,
		"reads-r2",
		"2",
		"",
		"FASTQ file containing the R2 reads (if paired-end)",
	)

	flags.StringVarP(
		&readOrigin,
		"read-origin",
		"r",
		"",
		"Specify read origin: 'dna' or 'rna' (required)",
	)
	mapCmd.MarkFlagRequired("read-origin")

	flags.StringVarP(
		&outputDirPath,
		"output",
		"o",
		"",
		"Output directory (required)",
	)
	mapCmd.MarkFlagRequired("output")

	flags.StringVarP(
		&samFileName,
		"sam-file-name",
		"f",
		"",
		"Output SAM file name (within output directory) (default: aligned.sam)",
	)

	flags.IntVarP(
		&threads,
		"threads",
		"t",
		0,
		"Number of threads to use (default: all available)",
	)

	flags.StringVarP(
		&statsFile,
		"stats",
		"s",
		"",
		"Mapping statistics file name (within output dir) or full path (default: mapping_stats.tsv)",
	)

	flags.StringVarP(
		&progressFile,
		"progress",
		"p",
		"",
		"Progress file name (within output dir) or full path (default: mapping_progress.tsv)",
	)

	flags.StringVarP(
		&progressStageTimeFile,
		"stagetime",
		"e",
		"",
		"Progress stage time file name (within output dir) or full path (default: mapping_stage.tsv)",
	)

	return mapCmd
}

func ExecMap() {
	indexFile, errIndexFile := os.Open(config.Mapper.Mapping.IndexFilePath)
	if errIndexFile != nil {
		logrus.Fatalf("Could not open index file: %v", errIndexFile)
	}
	defer indexFile.Close()

	genomeIndex := index.ReadGenomeIndexByFile(indexFile)

	// TODO: add regionmask and action support

	// // add region mask to index if provided
	// if argsObj.RegionmaskBedFile != nil && argsObj.RegionmaskPriorityFile != nil {
	//
	// 	fmt.Println("regionmask action:", *argsObj.RegionmaskAction)
	//
	// 	index.AddRegionmaskToIndex(
	// 		argsObj.RegionmaskBedFile,
	// 		argsObj.RegionmaskPriorityFile,
	// 		genomeIndex,
	// 	)
	//
	// 	logrus.WithFields(logrus.Fields{
	// 		"region mask":   argsObj.RegionmaskBedFile.Name(),
	// 		"priority file": argsObj.RegionmaskPriorityFile.Name(),
	// 	}).Info("Using region mask bed file and priority file")
	//
	// }

	fastqReader, errFastqReader := fastq.InitFromPaths(
		config.Mapper.Mapping.FastqR1FilePath,
		config.Mapper.Mapping.FastqR2FilePath,
	)
	if errFastqReader != nil {
		logrus.Fatalf("Could not initialize fastq reader: %v", errFastqReader)
	}

	samFile, errSamFile := config.Mapper.GetMappingOutputSamFile()
	if errSamFile != nil {
		logrus.Fatalf("Could not create/open output SAM file: %v", errSamFile)
	}

	writer := datawriter.InitFromFile(samFile)

	mapper.MapAll(
		genomeIndex,
		fastqReader,
		writer,
		config.Mapper.Mapping.Threads,
	)
}
