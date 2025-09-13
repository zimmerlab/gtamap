package runner

import (
	"fmt"
	"os"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/sirupsen/logrus"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
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

	viper := viper.New()

	// TODO: prescan os.Args for --key=value and set in viper
	// which allows to supply any argument even when not defined as flag

	// for _, arg := range os.Args[1:] {
	// 	if strings.HasPrefix(arg, "--") {
	// 		parts := strings.SplitN(arg[2:], "=", 2)
	// 		// map --foo-bar -> foo.bar
	// 		key := strings.ReplaceAll(parts[0], "-", ".")
	// 		value := ""
	// 		if len(parts) > 1 {
	// 			value = parts[1]
	// 		}
	// 		fmt.Println(key, " -> ", value)
	// 		// viper.Set(key, value)
	// 	}
	// }

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

			_, errConfig := config.LoadMapperConfig(viper, configFile)
			if errConfig != nil {
				return fmt.Errorf("failed to load config: %w", errConfig)
			}

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

	// subcommand: index-pre
	indexPreCmd := &cobra.Command{
		Use:   "index-pre",
		Short: "Prepare index (pre step)",
		Run: func(cmd *cobra.Command, args []string) {
			fmt.Println(args)
			// ExecIndexPre(args)
		},
	}

	// Subcommand: index-pre-region
	indexPreRegionCmd := &cobra.Command{
		Use:   "index-pre-region",
		Short: "Prepare index with region",
		Run: func(cmd *cobra.Command, args []string) {
			fmt.Println(args)
			// ExecIndexPreRegion(args)
		},
	}

	// Subcommand: index
	indexCmd := &cobra.Command{
		Use:   "index",
		Short: "Build index",
		Run: func(cmd *cobra.Command, args []string) {
			fmt.Println(args)
			// ExecIndex(args)
		},
	}

	mapCmd := GetCommandMap(viper)

	rootCmd.AddCommand(
		indexPreCmd,
		indexPreRegionCmd,
		indexCmd,
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
