package main

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper"
	"github.com/KleinSamuel/gtamap/src/dataloader"
	"github.com/KleinSamuel/gtamap/src/datawriter"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/sirupsen/logrus"
	"os"
	"strconv"
)

func buildAndSerializeIndexGenome(pathFasta string, pathOutput string) {

	fastaFile, errFasta := os.Open(pathFasta)
	if errFasta != nil {
		logrus.Fatal("Error reading fasta file", errFasta)
	}

	outputFile, errOutput := os.Create(pathOutput)
	if errOutput != nil {
		logrus.Fatal("Error creating output file (.fai)", errOutput)
	}

	fastaEntries, err := dataloader.ReadFasta(fastaFile)

	if err != nil {
		logrus.Fatal("Error extracting sequence from fasta file", err)
	}

	logrus.WithFields(logrus.Fields{
		"NumSequences": len(fastaEntries),
	}).Info("Read sequence(s) from fasta")

	for i, entry := range fastaEntries {
		logrus.WithFields(logrus.Fields{
			"Header": entry.Header,
			"Length": len(entry.Sequence),
		}).Info("Added sequence #" + strconv.Itoa(i+1))
	}

	genomeIndex := index.BuildGenomeIndex(fastaEntries)

	index.WriteGenomeIndex(genomeIndex, outputFile)
}

func deserializeGenomeIndex() *index.GenomeIndex {
	return index.ReadGenomeIndexByPath("./resources/ENSG00000173585.genome.gtai")
}

func testFastqReader2() {

	pathReadsFw := "../resources/reads/manual/reads_ccr9.1.fq"
	pathReadsRw := "../resources/reads/manual/reads_ccr9.2.fq"

	reader := fastq.InitFromPaths(&pathReadsFw, &pathReadsRw)

	read := reader.NextRead()
	fmt.Println(read.ReadR1)
	fmt.Println(read.ReadR2)
}

func extractGeneSequenceFromGtfAndFastaForIndex() {

	geneIds := make(map[string]struct{})
	geneIds["ENSG00000173585"] = struct{}{}

	gtfPath := "/home/sam/Data/gtamap/Homo_sapiens.GRCh38.113.chr.gtf"
	fastaPath := "/home/sam/Data/gtamap/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
	//fastaIndexPath := "/home/sam/Data/gtamap/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
	outputPath := "/home/sam/Data/gtamap/"

	index.ExtractGeneSequenceFromGtfAndFastaForIndex(gtfPath, fastaPath, outputPath, geneIds,
		0, 0, false)
}

func testTas2Read() {

	genomeIndexPath := "/home/sam/Data/gtamap/tas2/tas2r4/index/ENSG00000127364.gtai"

	//readsFwPath := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_1.tas2r4.fastq"
	//readsRvPath := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_2.tas2r4.fastq"

	readsFwPath := "/home/sam/Data/gtamap/tas2/tas2r4/bugs/018/r1.fastq"
	readsRvPath := "/home/sam/Data/gtamap/tas2/tas2r4/bugs/018/r2.fastq"

	outputPath := "/home/sam/Data/gtamap/tas2/tas2r4/bugs/018/aligned.sam"

	genomeIndex := index.ReadGenomeIndexByPath(genomeIndexPath)
	reader := fastq.InitFromPaths(&readsFwPath, &readsRvPath)
	writer := datawriter.InitFromPath(outputPath)

	numThreads := 1

	mapper.MapAll(genomeIndex, reader, writer, &numThreads)
}

func testTas2ReadsAll() {

	genomeIndexPath := "/home/sam/Data/gtamap/tas2/tas2r4/index/ENSG00000127364.gtai"

	readsFwPath := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_1.fastq.gz"
	readsRvPath := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_2.fastq.gz"

	outputPath := "/home/sam/Data/gtamap/tas2/tas2r4/aligned.sam"

	genomeIndex := index.ReadGenomeIndexByPath(genomeIndexPath)
	reader := fastq.InitFromPaths(&readsFwPath, &readsRvPath)
	writer := datawriter.InitFromPath(outputPath)

	numThreads := 1

	mapper.MapAll(genomeIndex, reader, writer, &numThreads)
}

func testTas2r4DeletionReads() {

	genomeIndexPath := "/home/sam/Data/gtamap/tas2/tas2r4/index/ENSG00000127364.gtai"

	readsFwPath := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_1.with-del.fastq"
	readsRvPath := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_2.with-del.fastq"

	outputPath := "/home/sam/Data/gtamap/tas2/tas2r4/aligned.with-del.sam"

	genomeIndex := index.ReadGenomeIndexByPath(genomeIndexPath)
	reader := fastq.InitFromPaths(&readsFwPath, &readsRvPath)
	writer := datawriter.InitFromPath(outputPath)

	numThreads := 1

	mapper.MapAll(genomeIndex, reader, writer, &numThreads)
}

func testIndex() {

	pattern := "CAAATAAATCTTAAAATTTTATAAATTACATGACTTTTCTCATT"

	var p1 [10]byte
	copy(p1[:], pattern[0:10])

	genomeIndexPath := "/home/sam/Data/gtamap/tas2/tas2r4/index/ENSG00000127364.gtai"

	genomeIndex := index.ReadGenomeIndexByPath(genomeIndexPath)

	matches := genomeIndex.GetKeywordFromMap(p1)

	for _, match := range matches {
		fmt.Println("Match found: ", match)
	}
}

func main() {

	//f, _ := os.Create("cpu_profile.prof")
	//err := pprof.StartCPUProfile(f)
	//if err != nil {
	//	return
	//}
	//defer pprof.StopCPUProfile()

	logrus.SetFormatter(&logrus.TextFormatter{
		FullTimestamp: true,
	})
	logrus.SetLevel(logrus.DebugLevel)

	//m := deserializeGenomeIndex()

	//genomeSeqPath := "./resources/ENSG00000173585.zeroed.fasta"
	//genomeIndexPath := "./resources/ENSG00000173585.genome.gtai"

	//genomeSeqPath := ""
	//genomeIndexPath := "./tas2r39.gtai"

	//outputPath := "./out/test.sam"

	//readsFwPath := "./resources/reads/manual/reads_ccr9.1.fq"
	//readsRwPath := "./resources/reads/manual/reads_ccr9.2.fq"
	//readsFwPath := "/home/sam/Data/reads/sim_ccr9/10mio/fw.fastq"
	//readsRvPath := "/home/sam/Data/reads/sim_ccr9/10mio/rw.fastq"

	//readsFwPath := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_1.10k.fastq"
	//readsRvPath := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_2.10k.fastq"

	//readsFwPath := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_1.fastq.gz"
	//readsRvPath := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_2.fastq.gz"

	//readsFwPath := "/home/sam/Data/genomes/test.r1.fastq"
	//readsRvPath := "/home/sam/Data/genomes/test.r2.fastq"

	//buildAndSerializeIndexGenome(genomeSeqPath, genomeIndexPath)

	//genomeIndex := index.ReadGenomeIndexByPath(genomeIndexPath)
	//reader := fastq.InitFromPaths(&readsFwPath, &readsRvPath)
	//writer := datawriter.InitFromPath(outputPath)
	//
	//numThreads := 1
	//
	//mapper.MapAll(genomeIndex, reader, writer, &numThreads)

	//extractGeneSequenceFromGtfAndFastaForIndex()

	testTas2Read()
	//testTas2r4DeletionReads()
	//testTas2ReadsAll()

	//testIndex()
}
