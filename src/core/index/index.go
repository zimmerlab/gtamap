package index

import (
	"encoding/gob"
	"github.com/KleinSamuel/gtamap/src/core/datastructure"
	"github.com/KleinSamuel/gtamap/src/dataloader"
	"github.com/KleinSamuel/gtamap/src/formats/gtf"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
	"log"
	"os"
	"time"
)

type GtaIndex struct {
	Gene                                    *Gene
	Transcripts                             []*Transcript
	SuffixTreeForwardStrandForwardDirection *datastructure.SuffixTree
	SuffixTreeForwardStrandReverseDirection *datastructure.SuffixTree
	SuffixTreeReverseStrandForwardDirection *datastructure.SuffixTree
	SuffixTreeReverseStrandReverseDirection *datastructure.SuffixTree
}

type Gene struct {
	GeneIdEnsembl  string // e.g. "ENSG00000173585"
	Chromosome     string // e.g. "1", "2", "X", "Y", "MT"
	IsFowardStrand bool   // true if on forward strand
	StartGenomic   uint32 // 0-based genomic start location
	EndGenomic     uint32 // exclusive genomic end location
}

type Transcript struct {
	TranscriptIdEnsembl                      string // e.g. "ENST00000342992"
	SequenceDnaForwardStrandForwardDirection string
	SequenceDnaForwardStrandReverseDirection string
	SequenceDnaReverseStrandForwardDirection string
	SequenceDnaReverseStrandReverseDirection string
	SequenceLength                           int
	Exons                                    []*Exon
}

type Exon struct {
	StartRelative uint32 // 0-based start location relative to genomic location of parent gene
	EndRelative   uint32 // exclusive end location relative to genomic location of parent gene
}

func BuildAndSerializeIndex(gtfFile *os.File, fastaFile *os.File, outputFile *os.File) {

	timerStart := time.Now()

	var fastIndexFilePath string = fastaFile.Name() + ".fai"
	dataloader.ExitIfFastaIndexIsMissing(fastIndexFilePath)
	fastaIndexFile, err := os.Open(fastIndexFilePath)
	if err != nil {
		log.Fatal("Error reading fasta index file (.fai)", err)
	}

	var annotation *gtf.Annotation = dataloader.GenerateInputForIndex(gtfFile, fastaFile, fastaIndexFile)

	var gene *Gene = &Gene{
		GeneIdEnsembl:  annotation.Genes[0].GeneIdEnsembl,
		Chromosome:     annotation.Genes[0].Chromosome,
		IsFowardStrand: annotation.Genes[0].IsForwardStrand,
	}

	var sequencesFwFw []string = make([]string, len(annotation.Genes[0].Transcripts))
	var sequencesFwRv []string = make([]string, len(annotation.Genes[0].Transcripts))
	var sequencesRvFw []string = make([]string, len(annotation.Genes[0].Transcripts))
	var sequencesRvRv []string = make([]string, len(annotation.Genes[0].Transcripts))

	var transcripts []*Transcript = make([]*Transcript, len(annotation.Genes[0].Transcripts))

	for i, transcript := range annotation.Genes[0].Transcripts {
		sequencesFwFw[i] = transcript.SequenceDna
		sequencesFwRv[i] = utils.ReverseString(transcript.SequenceDna)
		sequencesRvFw[i] = utils.ReverseComplementDNA(transcript.SequenceDna)
		sequencesRvRv[i] = utils.ComplementDNA(transcript.SequenceDna)

		transcripts[i] = &Transcript{
			TranscriptIdEnsembl:                      transcript.TranscriptIdEnsembl,
			SequenceLength:                           len(transcript.SequenceDna),
			SequenceDnaForwardStrandForwardDirection: sequencesFwFw[i],
			SequenceDnaForwardStrandReverseDirection: sequencesFwRv[i],
			SequenceDnaReverseStrandForwardDirection: sequencesRvFw[i],
			SequenceDnaReverseStrandReverseDirection: sequencesRvRv[i],
		}

		transcripts[i].Exons = make([]*Exon, len(transcript.Exons))

		for j, exon := range transcript.Exons {
			transcripts[i].Exons[j] = &Exon{
				StartRelative: exon.StartRelative,
				EndRelative:   exon.EndRelative,
			}
		}
	}

	// builds the suffix tree for the forward strand
	var suffixTreeFwFw *datastructure.SuffixTree = datastructure.BuildSuffixTree(sequencesFwFw)
	var suffixTreeFwRv *datastructure.SuffixTree = datastructure.BuildSuffixTree(sequencesFwRv)

	// builds the suffix tree for the reverse strand
	var suffixTreeRvFw *datastructure.SuffixTree = datastructure.BuildSuffixTree(sequencesRvFw)
	var suffixTreeRvRv *datastructure.SuffixTree = datastructure.BuildSuffixTree(sequencesRvRv)

	var gtaIndex GtaIndex = GtaIndex{
		Gene:                                    gene,
		Transcripts:                             transcripts,
		SuffixTreeForwardStrandForwardDirection: suffixTreeFwFw,
		SuffixTreeForwardStrandReverseDirection: suffixTreeFwRv,
		SuffixTreeReverseStrandForwardDirection: suffixTreeRvFw,
		SuffixTreeReverseStrandReverseDirection: suffixTreeRvRv,
	}

	SerializeFromFile(&gtaIndex, outputFile)

	logrus.WithFields(logrus.Fields{
		"total": time.Since(timerStart).String(),
	}).Info("Successfully built and serialized GTAMap index")
}

func SerializeFromPath(gtaIndex *GtaIndex, filePath string) {
	file, err := os.Create(filePath)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	SerializeFromFile(gtaIndex, file)
}

func SerializeFromFile(gtaIndex *GtaIndex, outputFile *os.File) {

	timerStart := time.Now()

	enc := gob.NewEncoder(outputFile)

	if err := enc.Encode(gtaIndex); err != nil {
		log.Fatal("encode error:", err)
	}

	logrus.WithFields(logrus.Fields{
		"duration": time.Since(timerStart),
		"output":   outputFile.Name(),
	}).Info("Serialized index")
}

func DezerializeFromPath(indexFilePath string) *GtaIndex {
	// Open the file for reading
	file, err := os.Open(indexFilePath)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	return DezerializeFromFile(file)
}

func DezerializeFromFile(indexFile *os.File) *GtaIndex {

	timerStart := time.Now()

	// Create a decoder and deserialize the person struct from the file
	decoder := gob.NewDecoder(indexFile)

	var gtaIndex GtaIndex

	if err := decoder.Decode(&gtaIndex); err != nil {
		panic(err)
	}

	logrus.WithFields(logrus.Fields{
		"duration": time.Since(timerStart),
	}).Info("Deserialized index")

	return &gtaIndex
}
