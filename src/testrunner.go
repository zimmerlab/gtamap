package main

import (
	"github.com/KleinSamuel/gtamap/src/datastructure"
	"github.com/sirupsen/logrus"
)

func main() {

	logrus.SetFormatter(&logrus.TextFormatter{
		FullTimestamp: true,
	})
	logrus.SetLevel(logrus.InfoLevel)

	sequences := []string{"abacad", "abd"}

	tree := datastructure.BuildSuffixTree(sequences)

	//tree.PrintEdgeList()

	tree.Search("a")
}
