package config

import (
	"fmt"
	"os"
	"path/filepath"

	"github.com/sirupsen/logrus"
	"github.com/spf13/viper"
)

var Mapper *MapperConfig

type MapperConfig struct {
	General struct {
		OutputDir string `mapstructure:"output_dir" yaml:"output_dir"`
	} `mapstructure:"general" yaml:"general"`

	Mapping struct {
		ReadOrigin      string `mapstructure:"read_origin" yaml:"read_origin"` // "dna" or "rna"
		IsReadOriginRna bool   `mapstructure:"-" yaml:"-"`                     // derived, not in config file
		Threads         int    `mapstructure:"threads" yaml:"threads"`

		RnaMode struct {
			IntronLengthMin       int     `mapstructure:"intron_length_min" yaml:"intron_length_min"`
			MaxMismatchCount      int     `mapstructure:"max_mismatch_count" yaml:"max_mismatch_count"`
			MaxMismatchPercentage float64 `mapstructure:"max_mismatch_percentage" yaml:"max_mismatch_percentage"`
		} `mapstructure:"rna_mode" yaml:"rna_mode"`

		DnaMode struct {
			MaxGapLength          int     `mapstructure:"max_gap_length" yaml:"max_gap_length"`
			MaxGapCount           int     `mapstructure:"max_gap_count" yaml:"max_gap_count"`
			MaxMismatchCount      int     `mapstructure:"max_mismatch_count" yaml:"max_mismatch_count"`
			MaxMismatchPercentage float64 `mapstructure:"max_mismatch_percentage" yaml:"max_mismatch_percentage"`
		} `mapstructure:"dna_mode" yaml:"dna_mode"`

		Output struct {
			IncludeMMinSAM     bool `mapstructure:"sam_use_x" yaml:"sam_use_x"`
			IncludeAllPairings bool `mapstructure:"sam_include_all_pairings" yaml:"sam_include_all_pairings"`
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

func setDefaults(v *viper.Viper) {

	// GENERAL
	v.SetDefault("general.output_dir", "./output")

	// MAPPING
	// should be required by parser (no default) such that the user has to
	// specify it which minimizes the chance of mistakes (hopefully)
	// v.SetDefault("mapping.read_origin", "dna")
	v.SetDefault("mapping.threads", -1)

	// MAPPING - RNA MODE
	v.SetDefault("mapping.rna_mode.intron_length_min", 20)
	// v.SetDefault("mapping.rna_mode.max_mismatch_count", 15)
	v.SetDefault("mapping.rna_mode.mismatch_max_percentage", 0.1)

	// MAPPING - DNA MODE
	v.SetDefault("mapping.dna_mode.max_gap_length", 1000)
	v.SetDefault("mapping.dna_mode.max_gap_count", 3)
	// v.SetDefault("mapping.dna_mode.max_mismatch_count", 15)
	v.SetDefault("mapping.dna_mode.max_mismatch_percentage", 0.1)

	// MAPPING - OUTPUT
	v.SetDefault("mapping.output.sam_use_x", true)
	v.SetDefault("mapping.output.sam_include_all_pairings", false)
}

func LoadMapperConfig(path string) (*MapperConfig, error) {

	v := viper.New()
	v.SetConfigName("mapperconfig")
	v.SetConfigType("yaml")

	if path != "" {
		// if a custom path is provided, use it (overrides default)
		v.SetConfigFile(path)
	} else {
		// use default: look for mapperconfig.yaml next to the binary
		execPath, err := os.Executable()
		if err != nil {
			return nil, fmt.Errorf("failed to get executable path: %w", err)
		}
		binaryDir := filepath.Dir(execPath)
		v.AddConfigPath(binaryDir)
		// fallback to current directory
		v.AddConfigPath("./")
	}

	setDefaults(v)

	if err := v.ReadInConfig(); err != nil {
		if path != "" {
			return nil, fmt.Errorf("failed to read config from %s: %w", path, err)
		}
		logrus.Info("No mapperconfig.yaml found, using default configuration")
	} else {
		logrus.WithFields(logrus.Fields{
			"config": v.ConfigFileUsed(),
		}).Info("Using mapper configuration")
	}

	var cfg MapperConfig
	if err := v.Unmarshal(&cfg); err != nil {
		logrus.Error(err)
		return nil, err
	}

	fmt.Println("no error")
	fmt.Println(cfg)

	Mapper = &cfg

	return Mapper, nil
}
