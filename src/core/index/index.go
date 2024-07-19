package index

import (
	"encoding/gob"
	"github.com/KleinSamuel/gtamap/src/core/datastructure"
	"github.com/KleinSamuel/gtamap/src/dataloader"
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

// GetTranscriptSequenceDna returns the DNA sequence of the transcript for the given sequence index
// sequenceIndex is the index of the sequence in the suffix tree
// the index of the transcript is sequenceIndex / 2
// if sequenceIndex is even, the sequence is on the forward strand, otherwise on the reverse strand
func (i GtaIndex) GetTranscriptSequenceDna(sequenceIndex int) string {
	//// the index of the transcript
	//var transcriptIndex uint8 = uint8(sequenceIndex / 2)
	//// true if forward strand, false if reverse strand
	//var isForward bool = sequenceIndex%2 == 0
	//// return the sequence of the transcript
	//return i.Transcripts[transcriptIndex].GetSequenceDna(isForward)
	return i.SuffixTree.Sequences[sequenceIndex]
}

type Gene struct {
	GeneIdEnsembl  string // e.g. "ENSG00000173585"
	Chromosome     string // e.g. "1", "2", "X", "Y", "MT"
	IsFowardStrand bool   // true if on forward strand
	StartGenomic   uint32 // 0-based genomic start location
	EndGenomic     uint32 // exclusive genomic end location
}

type Transcript struct {
	TranscriptIdEnsembl       string // e.g. "ENST00000342992"
	SequenceDnaForward53Index int    // the index of the forward sequence in the sequence list of the suffix tree
	SequenceDnaReverse53Index int    // the index in the reverse sequence in the sequence list of the suffix tree
	SequenceLength            int    // the length of the dna sequence
	Exons                     []*Exon
}

type Exon struct {
	StartRelative uint32 // 0-based start location relative to genomic location of parent gene
	EndRelative   uint32 // exclusive end location relative to genomic location of parent gene
}

func BuildAndSerializeIndex(gtfFile *os.File, fastaFile *os.File, outputFile *os.File) {

	timerStart := time.Now()

	timerReadReference := time.Now()
	durationReadReference := time.Duration(0)

	timerBuildTree := time.Now()
	durationBuildTree := time.Duration(0)

	timerReadReference = time.Now()
	var fastIndexFilePath string = fastaFile.Name() + ".fai"
	dataloader.ExitIfFastaIndexIsMissing(fastIndexFilePath)
	fastaIndexFile, err := os.Open(fastIndexFilePath)
	if err != nil {
		log.Fatal("Error reading fasta index file (.fai)", err)
	}
	durationReadReference += time.Since(timerReadReference)

	var gtaIndex = GtaIndex{
		Gene:         nil,
		Transcripts:  nil,
		SuffixTree:   nil,
		NumSequences: 0,
	}

	// read the reference files (gtf + fasta) and generate the annotation
	timerReadReference = time.Now()
	var annotation = dataloader.GenerateInputForIndex(gtfFile, fastaFile, fastaIndexFile)
	durationReadReference += time.Since(timerReadReference)

	gtaIndex.Gene = &Gene{
		GeneIdEnsembl:  annotation.Genes[0].GeneIdEnsembl,
		Chromosome:     annotation.Genes[0].Chromosome,
		IsFowardStrand: annotation.Genes[0].IsForwardStrand,
	}

	gtaIndex.Transcripts = make([]*Transcript, len(annotation.Genes[0].Transcripts))

	gtaIndex.SuffixTree = datastructure.CreateTree()

	for i, transcript := range annotation.Genes[0].Transcripts {

		//the index of the sequence in the suffix tree (used for retrieval)
		sequenceIndexForward := len(gtaIndex.SuffixTree.Sequences)
		timerBuildTree = time.Now()
		gtaIndex.SuffixTree.AddSequence(transcript.SequenceDna, sequenceIndexForward)
		durationBuildTree += time.Since(timerBuildTree)
		gtaIndex.NumSequences++

		// the index of the sequence in the suffix tree (used for retrieval)
		sequenceIndexReverse := len(gtaIndex.SuffixTree.Sequences)
		timerBuildTree = time.Now()
		gtaIndex.SuffixTree.AddSequence(utils.ReverseComplementDNA(transcript.SequenceDna), sequenceIndexReverse)
		durationBuildTree += time.Since(timerBuildTree)
		gtaIndex.NumSequences++

		gtaIndex.Transcripts[i] = &Transcript{
			TranscriptIdEnsembl:       transcript.TranscriptIdEnsembl,
			SequenceDnaForward53Index: sequenceIndexForward,
			SequenceDnaReverse53Index: sequenceIndexReverse,
			SequenceLength:            len(transcript.SequenceDna),
			Exons:                     make([]*Exon, len(transcript.Exons)),
		}

		for j, exon := range transcript.Exons {
			gtaIndex.Transcripts[i].Exons[j] = &Exon{
				StartRelative: exon.StartRelative,
				EndRelative:   exon.EndRelative,
			}
		}
	}

	// suffix links are only required for tree construction but not for the search
	timerBuildTree = time.Now()
	gtaIndex.SuffixTree.RemoveAllSuffixLinks()
	durationBuildTree += time.Since(timerBuildTree)

	// write the index to file as a serialized gob
	SerializeFromFile(&gtaIndex, outputFile)

	logrus.WithFields(logrus.Fields{
		"build": durationBuildTree.String(),
		"read":  durationReadReference.String(),
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

	// create a decoder and deserialize the person struct from the file
	decoder := gob.NewDecoder(indexFile)

	var gtaIndex GtaIndex

	if err := decoder.Decode(&gtaIndex); err != nil {
		panic(err)
	}

	logrus.WithFields(logrus.Fields{
		"duration": time.Since(timerStart),
	}).Info("Deserialized index")

	logrus.WithFields(logrus.Fields{
		"geneId":         gtaIndex.Gene.GeneIdEnsembl,
		"numTranscripts": len(gtaIndex.Transcripts),
	}).Info("Index info")

	// information about the transcripts contained in the index
	for i, transcript := range gtaIndex.Transcripts {
		logrus.WithFields(logrus.Fields{
			"sequenceIndex":  i,
			"transcriptId":   transcript.TranscriptIdEnsembl,
			"sequenceLength": transcript.SequenceLength,
		}).Info("Contained transcript")
	}

	return &gtaIndex
}
