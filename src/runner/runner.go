package runner

import (
	"fmt"
	"os"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/sirupsen/logrus"
	"github.com/spf13/cobra"
)

func PrintBanner() {
	logrus.Info("┌─┐┌┬┐┌─┐┌┬┐┌─┐┌─┐")
	logrus.Info("│ ┬ │ ├─┤│││├─┤├─┘")
	logrus.Info("└─┘ ┴ ┴ ┴┴ ┴┴ ┴┴")
	logrus.Info("gtamap v" + config.ToolVersion() + " (S. Klein, M. Weyrich, 2025)")
	logrus.Info("Sensitive read (DNA / RNA) mapping to a target region")
	logrus.Info("")
}

func Execute() {

	logrus.SetFormatter(
		&logrus.TextFormatter{
			FullTimestamp: true,
		},
	)

	PrintBanner()

	var configFile string
	var logLevel string

	rootCmd := &cobra.Command{
		Use: "gtamap",
		// Short: "Sensitive read mapping to a target region",
		PersistentPreRunE: func(cmd *cobra.Command, args []string) error {
			// global setup before any subcommand runs

			level, err := logrus.ParseLevel(logLevel)
			if err != nil {
				return fmt.Errorf("invalid log level: %w", err)
			}
			logrus.SetLevel(level)

			config.InitConfig(configFile)

			return nil
		},
	}

	// global flags (available to all subcommands)
	rootCmd.PersistentFlags().StringVar(
		&configFile,
		"config",
		"",
		"Path to config YAML",
	)
	rootCmd.PersistentFlags().StringVar(
		&logLevel,
		"loglevel",
		"INFO",
		"Log output level (ERROR, INFO, DEBUG)",
	)

	indexPreCmd := GetCommandIndexPre()
	indexPreRegionCmd := GetCommandIndexPreRegion()
	mapCmd := GetCommandMap()

	rootCmd.AddCommand(
		indexPreCmd,
		indexPreRegionCmd,
		mapCmd,
	)

	if err := rootCmd.Execute(); err != nil {
		os.Exit(1)
	}
}

// func Run() {
//
// 	logrus.SetFormatter(&logrus.TextFormatter{
// 		FullTimestamp: true,
// 	})
//
// 	PrintBanner()
//
// 	parser := argparse.NewParser(
// 		"gtamap",
// 		"Sensitive read mapping to a target region",
// 	)
//
// 	var configFile *string = parser.String(
// 		"",
// 		"config",
// 		&argparse.Options{
// 			Required: false,
// 			Help:     "Path to configuration file (YAML)",
// 			Default:  "",
// 		},
// 	)
//
// 	var logLevel *string = parser.Selector(
// 		"",
// 		"loglevel",
// 		[]string{"ERROR", "INFO", "DEBUG"},
// 		&argparse.Options{
// 			Required: false,
// 			Help:     "Log output level",
// 			Default:  "INFO",
// 		},
// 	)
//
// 	cmdIndexPre, argsIndexPre := AddCommandIndexPre(parser)
// 	cmdIndexPreRegion, argsIndexPreRegion := AddCommandIndexPreRegion(parser)
// 	cmdIndex, argsIndex := AddCommandIndex(parser)
// 	cmdMap, argsMap := AddCommandMap(parser)
//
// 	err := parser.Parse(os.Args)
// 	if err != nil {
// 		fmt.Printf("\n%s\n\n%s", err, parser.Usage(nil))
// 		os.Exit(1)
// 	}
//
// 	level, errLevel := logrus.ParseLevel(*logLevel)
// 	if errLevel != nil {
// 		logrus.Fatalf("Invalid log level: %v", errLevel)
// 	}
// 	logrus.SetLevel(level)
//
// 	_, errConfig := config.LoadMapperConfig(*configFile)
// 	if errConfig != nil {
// 		logrus.Fatalf("Failed to load configuration: %v", errConfig)
// 	}
//
// 	if cmdIndexPre.Happened() {
//
// 		ExecIndexPre(argsIndexPre)
//
// 	} else if cmdIndexPreRegion.Happened() {
//
// 		ExecIndexPreRegion(argsIndexPreRegion)
//
// 	} else if cmdIndex.Happened() {
//
// 		ExecIndex(argsIndex)
//
// 	} else if cmdMap.Happened() {
//
// 		ExecMap(argsMap)
//
// 	}
// }
