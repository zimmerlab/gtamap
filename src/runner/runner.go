package runner

import (
	"fmt"
	"os"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/akamensky/argparse"
	"github.com/sirupsen/logrus"
)

func PrintBanner() {
	logrus.Info("┌─┐┌┬┐┌─┐┌┬┐┌─┐┌─┐")
	logrus.Info("│ ┬ │ ├─┤│││├─┤├─┘")
	logrus.Info("└─┘ ┴ ┴ ┴┴ ┴┴ ┴┴")
	logrus.Info("gtamap v" + config.ToolVersion() + " (S. Klein, M. Weyrich, 2025)")
	logrus.Info("Sensitive read (DNA / RNA) mapping to a target region")
	logrus.Info("")
}

func Run() {

	logrus.SetFormatter(&logrus.TextFormatter{
		FullTimestamp: true,
	})

	PrintBanner()

	parser := argparse.NewParser(
		"gtamap",
		"Sensitive read mapping to a target region",
	)

	var configFile *string = parser.String(
		"",
		"config",
		&argparse.Options{
			Required: false,
			Help:     "Path to configuration file (YAML)",
			Default:  "",
		},
	)

	var logLevel *string = parser.Selector(
		"",
		"loglevel",
		[]string{"ERROR", "INFO", "DEBUG"},
		&argparse.Options{
			Required: false,
			Help:     "Log output level",
			Default:  "INFO",
		},
	)

	cmdIndexPre, argsIndexPre := AddCommandIndexPre(parser)
	cmdIndexPreRegion, argsIndexPreRegion := AddCommandIndexPreRegion(parser)
	cmdIndex, argsIndex := AddCommandIndex(parser)
	cmdMap, argsMap := AddCommandMap(parser)

	err := parser.Parse(os.Args)
	if err != nil {
		fmt.Printf("\n%s\n\n%s", err, parser.Usage(nil))
		os.Exit(1)
	}

	level, errLevel := logrus.ParseLevel(*logLevel)
	if errLevel != nil {
		logrus.Fatalf("Invalid log level: %v", errLevel)
	}
	logrus.SetLevel(level)

	_, errConfig := config.LoadMapperConfig(*configFile)
	if errConfig != nil {
		logrus.Fatalf("Failed to load configuration: %v", errConfig)
	}

	if cmdIndexPre.Happened() {

		ExecIndexPre(argsIndexPre)

	} else if cmdIndexPreRegion.Happened() {

		ExecIndexPreRegion(argsIndexPreRegion)

	} else if cmdIndex.Happened() {

		ExecIndex(argsIndex)

	} else if cmdMap.Happened() {

		ExecMap(argsMap)

	}
}
