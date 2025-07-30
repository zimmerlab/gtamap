package server

import (
	"fmt"
	"github.com/sirupsen/logrus"
	"github.com/spf13/viper"
	"log"
)

type Config struct {
	Data struct {
		RunDir           string `mapstructure:"run_dir"`
		TargetName       string `mapstructure:"target_name"`
		TargetRegion     string `mapstructure:"target_region"`
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
	fmt.Printf(ServerConfig.Data.TargetName)
	return ServerConfig.Data.TargetName
}

func GetTargetRegion() string {
	return ServerConfig.Data.TargetRegion
}

func GetTargetFasta() string {
	fmt.Printf(ServerConfig.Data.TargetFasta)
	return ServerConfig.Data.TargetFasta
}

func GetTargetFastaIndex() string {
	return ServerConfig.Data.TargetFastaIndex
}
