package config

import (
	"fmt"
	"io/fs"
	"os"
	"path/filepath"
	"runtime"
	"strings"

	"github.com/sirupsen/logrus"
	"github.com/spf13/viper"
)

var Mapper *MapperConfig

type MapperConfig struct {
	General struct {
		OutputDir string `mapstructure:"output_dir" yaml:"output_dir"`
	} `mapstructure:"general" yaml:"general"`

	Index struct {
		FastaFilePath       string   `mapstructure:"fasta_file_path" yaml:"fasta_file_path"`
		GenomeFastaFilePath string   `mapstructure:"genome_fasta_file_path" yaml:"genome_fasta_file_path"`
		FastaIndexFilePath  string   `mapstructure:"fasta_index_file_path" yaml:"fasta_index_file_path"`
		GtfFilePath         string   `mapstructure:"gtf_file_path" yaml:"gtf_file_path"`
		GeneIds             []string `mapstructure:"gene_ids" yaml:"gene_ids"`
		Regions             []string `mapstructure:"regions" yaml:"regions"`
		UpstreamBases       int      `mapstructure:"upstream_bases" yaml:"upstream_bases"`
		DownstreamBases     int      `mapstructure:"downstream_bases" yaml:"downstream_bases"`
		RegionmaskFilePath  string   `mapstructure:"regionmask_file_path" yaml:"regionmask_file_path"`

		Output struct {
			SingleFile       bool   `mapstructure:"single_file" yaml:"single_file"`
			FastaFileName    string `mapstructure:"fasta_file_name" yaml:"fasta_file_name"`
			IndexFileName    string `mapstructure:"index_file_name" yaml:"index_file_name"`
			UseFastaFileName bool   `mapstructure:"use_fasta_file_name" yaml:"use_fasta_file_name"`
		} `mapstructure:"output" yaml:"output"`
	} `mapstructure:"index" yaml:"index"`

	Mapping struct {
		IndexFilePath      string `mapstructure:"index_file_path" yaml:"index_file_path"`
		FastqR1FilePath    string `mapstructure:"fastq_r1_file_path" yaml:"fastq_r1_file_path"`
		FastqR2FilePath    string `mapstructure:"fastq_r2_file_path" yaml:"fastq_r2_file_path"`
		ReadOrigin         string `mapstructure:"read_origin" yaml:"read_origin"`
		IsReadOriginRna    bool   `mapstructure:"is_read_origin_rna" yaml:"is_read_origin_rna"`
		Threads            int    `mapstructure:"threads" yaml:"threads"`
		RegionmaskFilePath string `mapstructure:"regionmask_file_path" yaml:"regionmask_file_path"`
		MaxBranchingDepth  int    `mapstructure:"max_branching_depth" yaml:"max_branching_depth"`

		RnaMode struct {
			FilterMinMatches      int     `mapstructure:"filter_min_matches" yaml:"filter_min_matches"` // the minimum number of exact matches required to keep the read during filtering
			IntronLengthMin       int     `mapstructure:"intron_length_min" yaml:"intron_length_min"`
			MaxGapLength          int     `mapstructure:"max_gap_length" yaml:"max_gap_length"`
			MaxGapLengthTotal     int     `mapstructure:"max_gap_length_total" yaml:"max_gap_length_total"`
			MaxGapCount           int     `mapstructure:"max_gap_count" yaml:"max_gap_count"`
			MaxMismatchCount      int     `mapstructure:"max_mismatch_count" yaml:"max_mismatch_count"`
			MaxMismatchPercentage float64 `mapstructure:"max_mismatch_percentage" yaml:"max_mismatch_percentage"`

			Confident struct {
				MaxMismatchCount          int `mapstructure:"max_mismatch_count" yaml:"max_mismatch_count"`     // how many mm is a conf map allowed to have
				MinAnchorLength           int `mapstructure:"min_anchor_length" yaml:"min_anchor_length"`       // how long does each ali block in a conf map have to be considered conf
				IntronClusterDelta        int `mapstructure:"intron_cluster_delta" yaml:"intron_cluster_delta"` // by default, an intron cluster only absorbes an incoming gap (extending its reach) if the delta of gap.start/stop and cluster.start/stop is less than 100. This allows overlapping introns but also resolves intron coord confilct within close proximity
				IntronClusterRepairWindow int `mapstructure:"intron_cluster_repair_window" yaml:"intron_cluster_repair_window"`
			} `mapstructure:"confident" yaml:"confident"`
		} `mapstructure:"rna_mode" yaml:"rna_mode"`

		DnaMode struct {
			FilterMinMatches                   int     `mapstructure:"filter_min_matches" yaml:"filter_min_matches"` // the minimum number of exact matches required to keep the read during filtering
			MinLengthInitialDiagonal           int     `mapstructure:"min_length_initial_diagonal" yaml:"min_length_initial_diagonal"`
			MinLengthInitialDiagonalPercentage float64 `mapstructure:"min_length_initial_diagonal_percentage" yaml:"min_length_initial_diagonal_percentage"`
			MaxGapLength                       int     `mapstructure:"max_gap_length" yaml:"max_gap_length"`
			MaxGapCount                        int     `mapstructure:"max_gap_count" yaml:"max_gap_count"`
			MaxMismatchCount                   int     `mapstructure:"max_mismatch_count" yaml:"max_mismatch_count"`
			MaxMismatchPercentage              float64 `mapstructure:"max_mismatch_percentage" yaml:"max_mismatch_percentage"`

			Confident struct {
				MaxMismatchCount          int `mapstructure:"max_mismatch_count" yaml:"max_mismatch_count"`     // how many mm is a conf map allowed to have
				MinAnchorLength           int `mapstructure:"min_anchor_length" yaml:"min_anchor_length"`       // how long does each ali block in a conf map have to be considered conf
				IntronClusterDelta        int `mapstructure:"intron_cluster_delta" yaml:"intron_cluster_delta"` // by default, an intron cluster only absorbes an incoming gap (extending its reach) if the delta of gap.start/stop and cluster.start/stop is less than 100. This allows overlapping introns but also resolves intron coord confilct within close proximity
				IntronClusterRepairWindow int `mapstructure:"intron_cluster_repair_window" yaml:"intron_cluster_repair_window"`
			} `mapstructure:"confident" yaml:"confident"`
		} `mapstructure:"dna_mode" yaml:"dna_mode"`

		Output struct {
			MappingProgressFileName string `mapstructure:"mapping_progress_file_name" yaml:"mapping_progress_file_name"` // the name of the mapping progress log file
			MappingStatsFileName    string `mapstructure:"mapping_stats_file_name" yaml:"mapping_stats_file_name"`       // the name of the mapping stats file
			SamFileName             string `mapstructure:"sam_file_name" yaml:"sam_file_name"`                           // the default name of the output sam file
			IncludeMMinSAM          bool   `mapstructure:"sam_use_x" yaml:"sam_use_x"`                                   // if set to true, CIGAR will include "=" and "X" runes instead of only "M"
			IncludeAllPairings      bool   `mapstructure:"sam_include_all_pairings" yaml:"sam_include_all_pairings"`     // if set to true, all vs. all pairings are written to sam output  NOTE: actual pairing is not implemented yet
		} `mapstructure:"output" yaml:"output"`
	} `mapstructure:"mapping" yaml:"mapping"`
}

// SetReadOrigin sets the read origin and the derived boolean flag.
// Valid values are "dna" and "rna".
// Returns an error if the value is invalid.
func (c *MapperConfig) SetReadOrigin(value string) error {
	switch value {
	case "rna":
		c.Mapping.IsReadOriginRna = true
		c.Mapping.ReadOrigin = "rna"
		return nil
	case "dna":
		c.Mapping.IsReadOriginRna = false
		c.Mapping.ReadOrigin = "dna"
		return nil
	default:
		return fmt.Errorf("invalid read origin: %s (must be 'dna' or 'rna')",
			value)
	}
}

// GetMaxMismatches returns the maximum number of mismatches allowed for a read
// uses the max mismatch count if set (>=0), otherwise calculates it from the
// max mismatch percentage relative to the read length
func (c *MapperConfig) GetMaxMismatches(readLength int) int {
	if c.Mapping.IsReadOriginRna {
		if c.Mapping.RnaMode.MaxMismatchCount >= 0 {
			return c.Mapping.RnaMode.MaxMismatchCount
		} else {
			return int(float64(readLength) * c.Mapping.RnaMode.MaxMismatchPercentage)
		}
	} else {
		if c.Mapping.DnaMode.MaxMismatchCount >= 0 {
			return c.Mapping.DnaMode.MaxMismatchCount
		} else {
			return int(float64(readLength) * c.Mapping.DnaMode.MaxMismatchPercentage)
		}
	}
}

// SetMappingThreads sets the number of threads to use for mapping.
// If threads is <= 0 or greater than the number of CPU cores, it defaults to
// the number of CPU cores.
func (c *MapperConfig) SetMappingThreads(threads int) {
	cpuThreads := runtime.NumCPU()
	if threads <= 0 || threads > cpuThreads {
		c.Mapping.Threads = cpuThreads
	} else {
		c.Mapping.Threads = threads
	}
}

func (c *MapperConfig) GetMappingOutputSamFile() (*os.File, error) {
	// combine output dir and sam file name
	outputPath := filepath.Join(
		c.General.OutputDir,
		c.Mapping.Output.SamFileName,
	)

	file, err := os.OpenFile(outputPath, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0o600)
	if err != nil {
		return nil, fmt.Errorf("failed to open output SAM file %s: %w", outputPath, err)
	}

	return file, nil
}

func (c *MapperConfig) GetIndexOutputFile() (*os.File, error) {
	filename := ""

	if c.Index.Output.UseFastaFileName {

		base := filepath.Base(c.Index.FastaFilePath)
		name := strings.TrimSuffix(base, filepath.Ext(base))
		filename = name + ".gtai"

	} else {

		filename = c.Index.Output.IndexFileName

		if filepath.Ext(filename) != ".gtai" {
			filename += ".gtai"
		}

	}

	outputPath := filepath.Join(
		c.General.OutputDir,
		filename,
	)

	file, err := os.OpenFile(outputPath, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0o600)
	if err != nil {
		return nil, fmt.Errorf("failed to open output index file %s: %w", outputPath, err)
	}

	return file, nil
}

func (c *MapperConfig) SetIndexFastaIndex(fastaIndexFilePath string) {
	if fastaIndexFilePath != "" {
		c.Index.FastaIndexFilePath = fastaIndexFilePath
	} else {
		c.Index.FastaIndexFilePath = c.Index.FastaFilePath + ".fai"
	}
}

func (c *MapperConfig) GetConfidentMaxMismatchCount() int {
	if c.Mapping.IsReadOriginRna {
		return c.Mapping.RnaMode.Confident.MaxMismatchCount
	} else {
		return c.Mapping.DnaMode.Confident.MaxMismatchCount
	}
}

func (c *MapperConfig) GetConfidentMinAnchorLength() int {
	if c.Mapping.IsReadOriginRna {
		return c.Mapping.RnaMode.Confident.MinAnchorLength
	} else {
		return c.Mapping.DnaMode.Confident.MinAnchorLength
	}
}

func (c *MapperConfig) GetConfidentIntronClusterDelta() int {
	if c.Mapping.IsReadOriginRna {
		return c.Mapping.RnaMode.Confident.IntronClusterDelta
	} else {
		return c.Mapping.DnaMode.Confident.IntronClusterDelta
	}
}

func (c *MapperConfig) GetConfidentIntronClusterRepairWindow() int {
	if c.Mapping.IsReadOriginRna {
		return c.Mapping.RnaMode.Confident.IntronClusterRepairWindow
	} else {
		return c.Mapping.DnaMode.Confident.IntronClusterRepairWindow
	}
}

func (c *MapperConfig) GetFilterMinMatches() int {
	if c.Mapping.IsReadOriginRna {
		return c.Mapping.RnaMode.FilterMinMatches
	} else {
		return c.Mapping.DnaMode.FilterMinMatches
	}
}

func (c *MapperConfig) GetMapperProgressLogPath() string {
	if filepath.Dir(c.Mapping.Output.MappingProgressFileName) == "." {
		return filepath.Join(c.General.OutputDir, c.Mapping.Output.MappingProgressFileName)
	} else {
		return c.Mapping.Output.MappingProgressFileName
	}
}

func (c *MapperConfig) GetMappingStatsFilePath() string {
	if filepath.Dir(c.Mapping.Output.MappingStatsFileName) == "." {
		return filepath.Join(c.General.OutputDir, c.Mapping.Output.MappingStatsFileName)
	} else {
		return c.Mapping.Output.MappingStatsFileName
	}
}

func setDefaults() {
	// GENERAL
	viper.SetDefault("general.output_dir", "./output")

	// INDEX
	viper.SetDefault("index.gene_ids", []string{})
	viper.SetDefault("index.regions", []string{})
	viper.SetDefault("index.upstream_bases", 0)
	viper.SetDefault("index.downstream_bases", 0)

	// INDEX - OUTPUT
	viper.SetDefault("index.output.single_file", false)
	viper.SetDefault("index.output.fasta_file_name", "genes.fa")
	viper.SetDefault("index.output.index_file_name", "index.gtai")

	// MAPPING
	// should be required by parser (no default) such that the user has to
	// specify it which minimizes the chance of mistakes (hopefully)
	// viper.SetDefault("mapping.read_origin", "dna")
	viper.SetDefault("mapping.threads", -1)
	viper.SetDefault("mapping.max_branching_depth", 30)

	// MAPPING - RNA MODE
	viper.SetDefault("mapping.rna_mode.filter_min_matches", 6)
	viper.SetDefault("mapping.rna_mode.intron_length_min", 20)
	viper.SetDefault("mapping.rna_mode.max_mismatch_count", -1)
	viper.SetDefault("mapping.rna_mode.max_mismatch_percentage", 0.1)

	// MAPPING - RNA MODE - CONFIDENT
	viper.SetDefault("mapping.rna_mode.confident.max_mismatch_count", 6)
	viper.SetDefault("mapping.rna_mode.confident.min_anchor_length", 20)
	viper.SetDefault("mapping.rna_mode.confident.intron_cluster_delta", 100)
	viper.SetDefault("mapping.rna_mode.confident.intron_cluster_repair_window", 10)

	// MAPPING - DNA MODE
	viper.SetDefault("mapping.dna_mode.filter_min_matches", 7)
	viper.SetDefault("mapping.dna_mode.min_length_initial_diagonal", 0.7)
	viper.SetDefault("mapping.dna_mode.max_gap_length", 1000)
	viper.SetDefault("mapping.dna_mode.max_gap_count", 1)
	viper.SetDefault("mapping.dna_mode.max_mismatch_count", -1)
	viper.SetDefault("mapping.dna_mode.max_mismatch_percentage", 0.1)

	// MAPPING - DNA MODE - CONFIDENT
	viper.SetDefault("mapping.dna_mode.confident.max_mismatch_count", 6)
	viper.SetDefault("mapping.dna_mode.confident.min_anchor_length", 50)
	viper.SetDefault("mapping.dna_mode.confident.intron_cluster_delta", 100)
	viper.SetDefault("mapping.dna_mode.confident.intron_cluster_repair_window", 10)

	// MAPPING - OUTPUT
	viper.SetDefault("mapping.output.mapping_progress_file_name", "mapping_progress.log")
	viper.SetDefault("mapping.output.mapping_stats_file_name", "mapping_stats.tsv")
	viper.SetDefault("mapping.output.sam_file_name", "aligned.sam")
	viper.SetDefault("mapping.output.sam_use_x", true)
	viper.SetDefault("mapping.output.sam_include_all_pairings", false)
}

func InitConfig(configFilePath string) {
	setDefaults()

	// default search path is next to executable
	execPath, err := os.Executable()
	if err != nil {
		logrus.Fatalf("failed to get executable path: %v", err)
	}

	binaryDir := filepath.Dir(execPath)

	defaultConfigPath := filepath.Join(binaryDir, "mapperconfig.yaml")
	viper.SetConfigFile(defaultConfigPath)

	if err := viper.ReadInConfig(); err != nil {
		if _, ok := err.(*fs.PathError); !ok {
			logrus.Error("Error reading default configuration", err)
		} else {
			logrus.WithFields(logrus.Fields{
				"config": defaultConfigPath,
			}).Info("No default configuration file found")
		}
	} else {
		logrus.WithFields(logrus.Fields{
			"config": viper.ConfigFileUsed(),
		}).Info("Loaded default configuration file")
	}

	if configFilePath != "" {
		viper.SetConfigFile(configFilePath)
		if err := viper.MergeInConfig(); err != nil {
			logrus.Fatalf("Error reading config file %s: %v", configFilePath, err)
		} else {
			logrus.WithFields(logrus.Fields{
				"config": configFilePath,
			}).Info("Loaded custom configuration file")
		}
	}

	var cfg MapperConfig

	if err := viper.Unmarshal(&cfg); err != nil {
		logrus.Fatalf("Unable to decode config into struct: %v", err)
	}

	Mapper = &cfg
}
