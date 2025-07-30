package server

import (
	"github.com/spf13/viper"
	"log"
)

type Config struct {
	Server struct {
		Port int `mapstructure:"port"`
	} `mapstructure:"server"`
	Data struct {
		OutputDir string `mapstructure:"outputDir"`
	} `mapstructure:"data"`
}

var ServerConfig *Config

func LoadConfig() {
	viper.SetConfigName("config.server")
	viper.SetConfigType("toml")

	viper.AddConfigPath(".")
	//viper.AddConfigPath("./config")

	// environment variables
	//viper.AutomaticEnv()
	//viper.SetEnvPrefix("APP")

	// defaults
	viper.SetDefault("server.port", 8080)

	if err := viper.ReadInConfig(); err != nil {
		log.Printf("No config file found: %v", err)
	}

	ServerConfig = &Config{}
	if err := viper.Unmarshal(ServerConfig); err != nil {
		log.Fatalf("Unable to decode config: %v", err)
	}
}

func GetServerPort() int {
	return ServerConfig.Server.Port
}

func GetOutputDirectory() string {
	return ServerConfig.Data.OutputDir
}
