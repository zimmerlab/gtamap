package main

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/dataloader"
	"github.com/KleinSamuel/gtamap/src/gtf"
)

func main() {

	pathFastaCcr9Zeroed := "../resources/ENSG00000173585.fasta"
	pathGtfCcr9Zeroed := "../resources/ENSG00000173585.zeroed.gtf"

	var annotation *gtf.Annotation = dataloader.GenerateInputForIndex(pathGtfCcr9Zeroed, pathFastaCcr9Zeroed)

	fmt.Println("annotation: ", annotation)

}
