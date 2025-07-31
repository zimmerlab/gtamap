package config

import (
	"github.com/sirupsen/logrus"
	"github.com/spf13/viper"
	"log"
)

type Config struct {
	Data struct {
		RunDir           string `mapstructure:"run_dir"`
		TargetName       string `mapstructure:"target_name"`
		TargetRegion     string `mapstructure:"target_region"`
		TargetContig     string `mapstructure:"target_contig"`
		TargetStart      int    `mapstructure:"target_start"`
		TargetEnd        int    `mapstructure:"target_end"`
		TargetStrand     string `mapstructure:"target_strand"`
		TargetFasta      string `mapstructure:"target_fasta"`
		TargetFastaIndex string `mapstructure:"target_fasta_index"`
	} `mapstructure:"data"`
}

var ServerConfig *Config

func LoadConfig() {
	viper.SetConfigName("config.server")
	viper.SetConfigType("toml")

	//viper.AddConfigPath(".")
	//viper.AddConfigPath("./config")

	viper.AddConfigPath("/home/sam/Data/gtamap/output/3c2c2757/92f24301/664d462a/")

	// environment variables
	//viper.AutomaticEnv()
	//viper.SetEnvPrefix("APP")

	// defaults
	viper.SetDefault("server.port", 8080)

	if err := viper.ReadInConfig(); err != nil {
		logrus.Error("No config file found")
		logrus.Error(err)
	}

	ServerConfig = &Config{}
	if err := viper.Unmarshal(ServerConfig); err != nil {
		log.Fatalf("Unable to decode config: %v", err)
	}
}

func GetRunDir() string {
	return ServerConfig.Data.RunDir
}

func GetTargetName() string {
	return ServerConfig.Data.TargetName
}

func GetTargetRegion() string {
	return ServerConfig.Data.TargetRegion
}

func GetTargetContig() string {
	return ServerConfig.Data.TargetContig
}

func GetTargetStart() int {
	return ServerConfig.Data.TargetStart
}

func GetTargetEnd() int {
	return ServerConfig.Data.TargetEnd
}

func GetTargetStrand() string {
	return ServerConfig.Data.TargetStrand
}

func GetTargetFasta() string {
	return ServerConfig.Data.TargetFasta
}

func GetTargetFastaIndex() string {
	return ServerConfig.Data.TargetFastaIndex
}
