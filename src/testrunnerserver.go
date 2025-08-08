package main

import (
	"log"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/server"
)

func loadCcr9(s *server.Server) {

	gtamapTargetPath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.gtamap.target.sam"
	if err := s.AnalysisService.AddMapperInfo("gtamap", gtamapTargetPath, "target"); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	minimap2GenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.minimap2.genome.sam"
	if err := s.AnalysisService.AddMapperInfo("minimap2", minimap2GenomePath, "genome"); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	//starGenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.star.genome.sam"
	//if err := s.AnalysisService.AddMapperInfo("star", starGenomePath); err != nil {
	//	log.Fatalf("Error initializing analysis service: %v", err)
	//}
	//
	//hisat2GenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.hisat2.genome.sam"
	//if err := s.AnalysisService.AddMapperInfo("hisat2", hisat2GenomePath); err != nil {
	//	log.Fatalf("Error initializing analysis service: %v", err)
	//}
	//
	//bbmapGenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.bbmap.genome.sam"
	//if err := s.AnalysisService.AddMapperInfo("bbmap", bbmapGenomePath); err != nil {
	//	log.Fatalf("Error initializing analysis service: %v", err)
	//}
	//
	//bowtie2GenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.bowtie2.genome.sam"
	//if err := s.AnalysisService.AddMapperInfo("bowtie2", bowtie2GenomePath); err != nil {
	//	log.Fatalf("Error initializing analysis service: %v", err)
	//}
	//
	//bwaMemGenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.bwa-mem.genome.sam"
	//if err := s.AnalysisService.AddMapperInfo("bwa-mem", bwaMemGenomePath); err != nil {
	//	log.Fatalf("Error initializing analysis service: %v", err)
	//}
	//
	//bwaMem2GenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.bwa-mem2.genome.sam"
	//if err := s.AnalysisService.AddMapperInfo("bwa-mem2", bwaMem2GenomePath); err != nil {
	//	log.Fatalf("Error initializing analysis service: %v", err)
	//}
	//
	//contextmap2bowtieGenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.contextmap2-bowtie2.genome.sam"
	//if err := s.AnalysisService.AddMapperInfo("contextmap2-bowtie2", contextmap2bowtieGenomePath); err != nil {
	//	log.Fatalf("Error initializing analysis service: %v", err)
	//}
	//
	//contextmap2bwaGenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.contextmap2-bwa.genome.sam"
	//if err := s.AnalysisService.AddMapperInfo("contextmap2-bwa", contextmap2bwaGenomePath); err != nil {
	//	log.Fatalf("Error initializing analysis service: %v", err)
	//}
}

func loadNsun5(s *server.Server) {

	gtamapTargetPath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000130305/ENSG00000130305.gtamap.target.sam"
	if err := s.AnalysisService.AddMapperInfo("gtamap", gtamapTargetPath, "target"); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	minimap2GenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000130305/ENSG00000130305.minimap2.genome.sam"
	if err := s.AnalysisService.AddMapperInfo("minimap2", minimap2GenomePath, "genome"); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	starGenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000130305/ENSG00000130305.star.genome.sam"
	if err := s.AnalysisService.AddMapperInfo("star", starGenomePath, "genome"); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	hisat2GenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000130305/ENSG00000130305.hisat2.genome.sam"
	if err := s.AnalysisService.AddMapperInfo("hisat2", hisat2GenomePath, "genome"); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	bbmapGenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000130305/ENSG00000130305.bbmap.genome.sam"
	if err := s.AnalysisService.AddMapperInfo("bbmap", bbmapGenomePath, "genome"); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	bowtie2GenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000130305/ENSG00000130305.bowtie2.genome.sam"
	if err := s.AnalysisService.AddMapperInfo("bowtie2", bowtie2GenomePath, "genome"); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	bwaMemGenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000130305/ENSG00000130305.bwa-mem.genome.sam"
	if err := s.AnalysisService.AddMapperInfo("bwa-mem", bwaMemGenomePath, "genome"); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	bwaMem2GenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000130305/ENSG00000130305.bwa-mem2.genome.sam"
	if err := s.AnalysisService.AddMapperInfo("bwa-mem2", bwaMem2GenomePath, "genome"); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	contextmap2bowtieGenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000130305/ENSG00000130305.contextmap2-bowtie2.genome.sam"
	if err := s.AnalysisService.AddMapperInfo("contextmap2-bowtie2", contextmap2bowtieGenomePath, "genome"); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}

	contextmap2bwaGenomePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000130305/ENSG00000130305.contextmap2-bwa.genome.sam"
	if err := s.AnalysisService.AddMapperInfo("contextmap2-bwa", contextmap2bwaGenomePath, "genome"); err != nil {
		log.Fatalf("Error initializing analysis service: %v", err)
	}
}

func main() {

	config.LoadConfig()

	// s := server.NewServer()
	// s.LoadMapperResults()
	//
	// s.Start()

	fastaFilePath := "/home/sam/Data/reference-genomes/ensembl/113/homo_sapiens/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
	fastaIndexFilePath := "/home/sam/Data/reference-genomes/ensembl/113/homo_sapiens/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
	// gtamapSamFilePath := "/home/sam/Data/gtamap/output/3c2c2757/92f24301/664d462a/ENSG00000173585.gtamap.target.sam"

	s := server.NewServer()
	s.InitMappingDataHandler(fastaFilePath, fastaIndexFilePath)

	s.Handler.LoadMapperResults()

	s.Start()
}
