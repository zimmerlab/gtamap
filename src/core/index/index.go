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
	Gene         *Gene                     // information about the gene
	Transcripts  []*Transcript             // all transcripts of the gene
	SuffixTree   *datastructure.SuffixTree // a single suffix tree containing all transcripts
	NumSequences int
}

func (i GtaIndex) GetTranscriptSequenceDna(sequenceIndex int) string {
	transcriptIndex := sequenceIndex / 4
	sequenceIndexInTranscript := sequenceIndex % 4

	return i.Transcripts[transcriptIndex].GetSequenceDna(sequenceIndexInTranscript)
}

type Gene struct {
	GeneIdEnsembl  string // e.g. "ENSG00000173585"
	Chromosome     string // e.g. "1", "2", "X", "Y", "MT"
	IsFowardStrand bool   // true if on forward strand
	StartGenomic   uint32 // 0-based genomic start location
	EndGenomic     uint32 // exclusive genomic end location
}

type Transcript struct {
	TranscriptIdEnsembl  string // e.g. "ENST00000342992"
	SequenceDnaForward53 string // the DNA sequence on the forward strand in 5'-3' direction
	SequenceDnaForward35 string // the DNA sequence on the forward strand in 3'-5' direction
	SequenceDnaReverse53 string // the DNA sequence on the reverse strand in 5'-3' direction (rev-comp of forward)
	SequenceDnaReverse35 string // the DNA sequence on the reverse strand in 3'-5' direction (rev-comp of forward)
	SequenceLength       int
	Exons                []*Exon
}

func (t Transcript) GetSequenceDna(index int) string {
	switch index {
	case 0:
		return t.SequenceDnaForward53
	case 1:
		return t.SequenceDnaForward35
	case 2:
		return t.SequenceDnaReverse53
	case 3:
		return t.SequenceDnaReverse35
	}
	return ""
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

	var sequencesForward []string = make([]string, len(annotation.Genes[0].Transcripts))
	var sequencesReverseComplement []string = make([]string, len(annotation.Genes[0].Transcripts))

	var transcripts []*Transcript = make([]*Transcript, len(annotation.Genes[0].Transcripts))

	for i, transcript := range annotation.Genes[0].Transcripts {

		sequencesForward[i] = transcript.SequenceDna
		sequencesReverseComplement[i] = utils.ReverseComplementDNA(transcript.SequenceDna)

		transcripts[i] = &Transcript{
			TranscriptIdEnsembl:  transcript.TranscriptIdEnsembl,
			SequenceLength:       len(transcript.SequenceDna),
			SequenceDnaForward53: sequencesForward[i],
			SequenceDnaReverse53: sequencesReverseComplement[i],
		}

		transcripts[i].Exons = make([]*Exon, len(transcript.Exons))

		for j, exon := range transcript.Exons {
			transcripts[i].Exons[j] = &Exon{
				StartRelative: exon.StartRelative,
				EndRelative:   exon.EndRelative,
			}
		}
	}

	// contains all forward sequences of the transcripts and then all rev-comp sequences
	var allSequences []string = make([]string, 0)
	allSequences = append(allSequences, sequencesForward...)
	allSequences = append(allSequences, sequencesReverseComplement...)

	timerBuildTree := time.Now()

	var suffixTree = datastructure.CreateTree()

	// add each sequence to the suffix tree
	for sequenceIndex, sequence := range allSequences {
		suffixTree.AddSequence(sequence, sequenceIndex)
	}
	// suffix links are only required for tree construction but not for the search
	suffixTree.RemoveAllSuffixLinks()

	totalBuildTree := time.Since(timerBuildTree)

	var gtaIndex GtaIndex = GtaIndex{
		Gene:         gene,
		Transcripts:  transcripts,
		SuffixTree:   suffixTree,
		NumSequences: len(allSequences),
	}

	SerializeFromFile(&gtaIndex, outputFile)

	logrus.WithFields(logrus.Fields{
		"buildTree": totalBuildTree.String(),
		"total":     time.Since(timerStart).String(),
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
