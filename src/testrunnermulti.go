package main

import (
	"os"
	"runtime/pprof"
	"strconv"

	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper"
	"github.com/KleinSamuel/gtamap/src/dataloader"
	"github.com/KleinSamuel/gtamap/src/datawriter"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/sirupsen/logrus"
)

func buildMultiFastaIndex(pathFasta string, pathOutput string, pathRepeatmasker string) {

	fastaFile, errFasta := os.Open(pathFasta)
	if errFasta != nil {
		logrus.Fatal("Error reading fasta file", errFasta)
	}

	outputFile, errOutput := os.Create(pathOutput)
	if errOutput != nil {
		logrus.Fatal("Error creating output file (.fai)", errOutput)
	}

	repeatmaskFile, errRm := os.Open(pathRepeatmasker)
	if errRm != nil {
		logrus.Fatal("Error reading repeatmasker file", errRm)
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

	genomeIndex := index.BuildGenomeIndex(fastaEntries, repeatmaskFile)
	// genomeIndex := index.BuildGenomeIndex(fastaEntries, nil)

	// if 1 == 2 {
	// 	fmt.Println(repeatmaskFile)
	// }

	index.WriteGenomeIndex(genomeIndex, outputFile)
}

func mapTas2R4All(genomeIndexPath string, outputPath string) {

	// readsFwPath := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_1.fastq.gz"
	// readsRvPath := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_2.fastq.gz"
	readsFwPath := "/home/sam/Data/genome-sam/30m.fw.fq"
	readsRvPath := "/home/sam/Data/genome-sam/30m.rv.fq"

	genomeIndex := index.ReadGenomeIndexByPath(genomeIndexPath)
	reader, _ := fastq.InitFromPaths(&readsFwPath, &readsRvPath)
	writer := datawriter.InitFromPath(outputPath)

	numThreads := 18

	mapper.MapAll(genomeIndex, reader, writer, &numThreads)
}

func multiBuildSingle() {
	pathRepeatmasker := "/home/sam/Data/repeatmasker/hg38.fa.out"
	pathFastaMulti := "/home/sam/Data/gtamap/multi-index/ENSG00000075624.fa"
	pathOutputMulti := "/home/sam/Data/gtamap/multi-index/ENSG00000075624.gtai"

	buildMultiFastaIndex(pathFastaMulti, pathOutputMulti, pathRepeatmasker)
}

func multiMapSingle() {
	pathOutputMulti := "/home/sam/Data/gtamap/multi-index/ENSG00000075624.gtai"
	pathOutputMapping := "/home/sam/Data/gtamap/multi-index/ENSG00000075624.sam"
	mapTas2R4All(pathOutputMulti, pathOutputMapping)
}

func multiBuildAll() {
	pathRepeatmasker := "/home/sam/Data/repeatmasker/hg38.fa.out"
	pathFastaMulti := "/home/sam/Data/gtamap/multi-index/genes.fa"
	pathOutputMulti := "/home/sam/Data/gtamap/multi-index/genes.gtai"
	buildMultiFastaIndex(pathFastaMulti, pathOutputMulti, pathRepeatmasker)
}

func multiMapAll() {
	pathOutputMulti := "/home/sam/Data/gtamap/multi-index/genes.gtai"
	pathOutputMapping := "/home/sam/Data/gtamap/multi-index/genes.sam"
	mapTas2R4All(pathOutputMulti, pathOutputMapping)
}

func main() {
	logrus.SetFormatter(&logrus.TextFormatter{
		FullTimestamp: true,
	})
	logrus.SetLevel(logrus.InfoLevel)

	ff, errr := os.Create("cpu.pprof")
	if errr != nil {
		panic(errr)
	}
	defer ff.Close()
	if err := pprof.StartCPUProfile(ff); err != nil {
		panic(err)
	}
	defer pprof.StopCPUProfile()

	// multiBuildSingle()
	multiMapSingle()
	// multiBuildAll()
	// multiMapAll()
}
