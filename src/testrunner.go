package main

import "github.com/KleinSamuel/gtamap/src/datastructure"

func main() {

	sequences := []string{"abcabd", "abd"}

	tree := datastructure.BuildSuffixTree(sequences)
	tree.PrintEdgeList()
}
