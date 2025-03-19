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

func buildAndSerializeIndexGenome() {

	pathFasta := "./resources/ENSG00000173585.zeroed.fasta"
	//pathFasta := "./resources/test.2.fasta"
	pathOutput := "./resources/ENSG00000173585.genome.gtai"

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
	logrus.SetLevel(logrus.ErrorLevel)

	//buildAndSerializeIndexGenome()
	//m := deserializeGenomeIndex()

	genomeIndexPath := "./resources/ENSG00000173585.genome.gtai"
	outputPath := "./out/test.sam"
	//readsFwPath := "./resources/reads/manual/reads_ccr9.1.fq"
	//readsRwPath := "./resources/reads/manual/reads_ccr9.2.fq"
	readsFwPath := "/home/sam/Data/reads/sim_ccr9/10mio/fw.fastq"
	readsRvPath := "/home/sam/Data/reads/sim_ccr9/10mio/rw.fastq"

	genomeIndex := index.ReadGenomeIndexByPath(genomeIndexPath)
	reader := fastq.InitFromPaths(&readsFwPath, &readsRvPath)
	writer := datawriter.InitFromPath(outputPath)

	mapper.MapAll(genomeIndex, reader, writer)

}
