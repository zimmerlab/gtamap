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

	starGenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.star.genome.sam"
	if err := s.AnalysisService.AddMapperInfo("star", starGenomePath); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	hisat2GenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.hisat2.genome.sam"
	if err := s.AnalysisService.AddMapperInfo("hisat2", hisat2GenomePath); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	bbmapGenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.bbmap.genome.sam"
	if err := s.AnalysisService.AddMapperInfo("bbmap", bbmapGenomePath); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	bowtie2GenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.bowtie2.genome.sam"
	if err := s.AnalysisService.AddMapperInfo("bowtie2", bowtie2GenomePath); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	bwaMemGenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.bwa-mem.genome.sam"
	if err := s.AnalysisService.AddMapperInfo("bwa-mem", bwaMemGenomePath); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	bwaMem2GenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.bwa-mem2.genome.sam"
	if err := s.AnalysisService.AddMapperInfo("bwa-mem2", bwaMem2GenomePath); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	contextmap2bowtieGenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.contextmap2-bowtie2.genome.sam"
	if err := s.AnalysisService.AddMapperInfo("contextmap2-bowtie2", contextmap2bowtieGenomePath); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	contextmap2bwaGenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.contextmap2-bwa.genome.sam"
	if err := s.AnalysisService.AddMapperInfo("contextmap2-bwa", contextmap2bwaGenomePath); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	s.Start()
}
