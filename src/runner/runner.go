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
	logrus.Info("Fast and memory efficient spliced read mapping to a single gene reference.")
	logrus.Info("")
}

func Run() {

	logrus.SetFormatter(&logrus.TextFormatter{
		FullTimestamp: true,
	})

	PrintBanner()

	parser := argparse.NewParser(
		"gtamap",
		"Gene-centric spliced read mapping",
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
