package config

import (
	"github.com/sirupsen/logrus"
	"github.com/spf13/viper"
)

type MapperConfig struct {
	General struct {
		OutputDir string `yaml:"output_dir"`
	} `yaml:"general"`

	Mapping struct {
		ReadOrigin string `yaml:"read_origin"` // "dna" or "rna"

		IntronLengthMin int `yaml:"intron_length_min"`

		MismatchHardCutoff    int     `yaml:"mismatch_hard_cutoff"`
		MismatchMaxPercentage float64 `yaml:"mismatch_max_percentage"`

		Output struct {
			IncludeMMinSAM     bool `yaml:"sam_use_x"`
			IncludeAllPairings bool `yaml:"sam_include_all_pairings"`
		} `yaml:"output"`
	} `yaml:"mapping"`
}

func setDefaults(v *viper.Viper) {
	// GENERAL
	v.SetDefault("general.output_dir", "./output")
	// MAPPING
	v.SetDefault("mapping.read_origin", "dna")
	v.SetDefault("mapping.intron_length_min", 20)
	v.SetDefault("mapping.mismatch_hard_cutoff", 5)
	v.SetDefault("mapping.mismatch_max_percentage", 0.2)
	// MAPPING OUTPUT
	v.SetDefault("mapping.output.sam_use_x", true)
	v.SetDefault("mapping.output.sam_include_all_pairings", false)
}

func LoadMapperConfig() (*MapperConfig, error) {

	v := viper.New()
	v.SetConfigName("mapperconfig")
	v.SetConfigType("yaml")
	v.AddConfigPath("./")

	setDefaults(v)

	if err := v.ReadInConfig(); err != nil {
		logrus.Info("No mapperconfig.yaml found, using default configuration")
	} else {
		logrus.WithFields(logrus.Fields{
			"config": v.ConfigFileUsed(),
		}).Info("Using mapper configuration")
	}

	var cfg MapperConfig
	if err := v.Unmarshal(&cfg); err != nil {
		return nil, err
	}

	return &cfg, nil
}
