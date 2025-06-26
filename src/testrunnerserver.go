package main

import (
	"github.com/KleinSamuel/gtamap/src/server"
	"log"
)

func main() {

	s := server.NewServer()

	gtamapTargetPath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.gtamap.target.sam"

	if err := s.AnalysisService.AddMapperInfo("gtamap", gtamapTargetPath); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	minimap2GenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.minimap2.genome.sam"

	if err := s.AnalysisService.AddMapperInfo("minimap2", minimap2GenomePath); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	s.Start()

}
