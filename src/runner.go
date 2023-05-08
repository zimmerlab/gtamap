package main

import (
	"github.com/sirupsen/logrus"
)

func printBanner() {
	logrus.Info("_____ _____ ___   ___  ___  ___  ______")
	logrus.Info("|  __ \\_   _/ _ \\  |  \\/  | / _ \\ | ___ \\")
	logrus.Info("| |  \\/ | |/ /_\\ \\ | .  . |/ /_\\ \\| |_/ /")
	logrus.Info("| | __  | ||  _  | | |\\/| ||  _  ||  __/")
	logrus.Info("| |_\\ \\ | || | | | | |  | || | | || |")
	logrus.Info(" \\____/ \\_/\\_| |_/ \\_|  |_/\\_| |_/\\_|")
	logrus.Info("Starting GTAMap v0.1.0 (A. Hadziahmetovic, S. Klein, 2023)")
}

func main() {

	logrus.SetFormatter(&logrus.TextFormatter{
		FullTimestamp: true,
	})
	logrus.SetLevel(logrus.DebugLevel)

	printBanner()

	logrus.Info("Runner does nothing for now..")

}
