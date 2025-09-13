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

// type ArgsMap struct {
// 	IndexFile              *os.File
// 	FastqR1File            *os.File
// 	FastqR2File            *os.File
// 	OutputFile             *string
// 	ReadType               *string
// 	Threads                *int
// 	RegionmaskBedFile      *os.File
// 	RegionmaskPriorityFile *os.File
// 	RegionmaskAction       *string
// }

func GetCommandMap(v *viper.Viper) *cobra.Command {

	var indexFilePath string
	var fastqR1FilePath string
	var fastqR2FilePath string
	var readOrigin string

	var outputDirPath string
	var samFileName string
	var threads int

	mapCmd := &cobra.Command{
		Use:   "map",
		Short: "Run mapping",
		Run: func(cmd *cobra.Command, args []string) {

			errOrigin := config.Mapper.SetReadOrigin(config.Mapper.Mapping.ReadOrigin)
			if errOrigin != nil {
				cmd.PrintErrln("\n", cmd.ErrPrefix(), "[-r | --read-origin]",
					errOrigin.Error(), "\n")
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
	v.BindPFlag(
		"mapping.index_file_path",
		mapCmd.Flags().Lookup("index"),
	)

	flags.StringVarP(
		&fastqR1FilePath,
		"reads-r1",
		"1",
		"",
		"FASTQ file containing the R1 reads (required)",
	)
	mapCmd.MarkFlagRequired("reads-r1")
	v.BindPFlag(
		"mapping.fastq_r1_file_path",
		mapCmd.Flags().Lookup("reads-r1"),
	)

	flags.StringVarP(
		&fastqR2FilePath,
		"reads-r2",
		"2",
		"",
		"FASTQ file containing the R2 reads (if paired-end)",
	)
	v.BindPFlag(
		"mapping.fastq_r2_file_path",
		mapCmd.Flags().Lookup("reads-r2"),
	)

	flags.StringVarP(
		&readOrigin,
		"read-origin",
		"r",
		"",
		"Specify read origin: 'dna' or 'rna' (required)",
	)
	mapCmd.MarkFlagRequired("read-origin")
	v.BindPFlag(
		"mapping.read_origin",
		mapCmd.Flags().Lookup("read-origin"),
	)

	flags.StringVarP(
		&outputDirPath,
		"output",
		"o",
		"",
		"Output directory (required)",
	)
	mapCmd.MarkFlagRequired("output")
	v.BindPFlag(
		"general.output_dir",
		mapCmd.Flags().Lookup("output"),
	)

	flags.StringVarP(
		&samFileName,
		"sam-file-name",
		"f",
		"",
		"Output SAM file name (within output directory) (default: aligned.sam)",
	)
	v.BindPFlag(
		"mapping.output.sam_file_name",
		mapCmd.Flags().Lookup("sam-file-name"),
	)

	flags.IntVarP(
		&threads,
		"threads",
		"t",
		0,
		"Number of threads to use (default: all available)",
	)
	v.BindPFlag(
		"mapping.threads",
		mapCmd.Flags().Lookup("threads"),
	)

	return mapCmd
}

// func AddCommandMap(
// 	parser *argparse.Parser,
// ) (
// 	*argparse.Command,
// 	*ArgsMap,
// ) {
//
// 	var command *argparse.Command = parser.NewCommand(
// 		"map",
// 		"Map reads to the GTAMap index",
// 	)
//
// 	argsObj := &ArgsMap{}
//
// 	argsObj.IndexFile = command.File(
// 		"",
// 		"index",
// 		os.O_RDONLY,
// 		0o600,
// 		&argparse.Options{
// 			Required: true,
// 			Help:     "GTAMap index (.gtai) file",
// 		},
// 	)
//
// 	argsObj.FastqR1File = command.File(
// 		"",
// 		"r1",
// 		os.O_RDONLY,
// 		0o600,
// 		&argparse.Options{
// 			Required: true,
// 			Help:     "FASTQ file containing the forward reads",
// 		},
// 	)
//
// 	argsObj.FastqR2File = command.File(
// 		"",
// 		"r2",
// 		os.O_RDONLY,
// 		0o600,
// 		&argparse.Options{
// 			Required: false,
// 			Help:     "FASTQ file containing the reverse reads",
// 		},
// 	)
//
// 	argsObj.OutputFile = command.String(
// 		"",
// 		"output",
// 		&argparse.Options{
// 			Required: true,
// 			Help:     "Output SAM file (file extension: .sam)",
// 			Default:  "",
// 		},
// 	)
//
// 	// var logFileMap *string = cmdMap.String("", "log", &argparse.Options{
// 	// 	Help:    "Output log file.",
// 	// 	Default: "",
// 	// })
//
// 	argsObj.ReadType = command.Selector(
// 		"",
// 		"read-origin",
// 		[]string{"dna", "rna"},
// 		&argparse.Options{
// 			Required: true,
// 			Help:     "Specify read origin",
// 		},
// 	)
//
// 	argsObj.Threads = command.Int(
// 		"",
// 		"threads",
// 		&argparse.Options{
// 			Required: false,
// 			Help:     "Number of threads to use (default: all)",
// 			Default:  -1,
// 		},
// 	)
//
// 	argsObj.RegionmaskBedFile = command.File(
// 		"",
// 		"regionmask-bed",
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
// 		"regionmask-priorities",
// 		os.O_RDONLY,
// 		0o600,
// 		&argparse.Options{
// 			Required: false,
// 			Help: "Path to TSV file containing region names and their " +
// 				"corresponding priority (higher number = higher priority)",
// 		},
// 	)
//
// 	argsObj.RegionmaskAction = command.Selector(
// 		"",
// 		"regionmask-action",
// 		[]string{"combine", "ignore-index"},
// 		&argparse.Options{
// 			Required: false,
// 			Help:     "Specify how to handle regionmasks",
// 			Default:  "combine",
// 		},
// 	)
//
// 	return command, argsObj
// }

func ExecMap() {

	// if *logFileMap != "" {
	// 	config.LogOut = *logFileMap
	// }

	// should already be captured by argument parser
	// errOrigin := config.Mapper.SetReadOrigin(*argsObj.ReadType)
	// if errOrigin != nil {
	// 	log.Fatal(errOrigin)
	// }

	indexFile, errIndexFile := os.Open(config.Mapper.Mapping.IndexFilePath)
	if errIndexFile != nil {
		logrus.Fatalf("Could not open index file: %v", errIndexFile)
	}

	genomeIndex := index.ReadGenomeIndexByFile(indexFile)

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
	//
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

	// reader, errFastq := fastq.InitFromFiles(
	// 	argsObj.FastqR1File,
	// 	argsObj.FastqR2File,
	// )

	fastqReader, errFastqReader := fastq.InitFromPaths(
		&config.Mapper.Mapping.FastqR1FilePath,
		&config.Mapper.Mapping.FastqR2FilePath,
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
