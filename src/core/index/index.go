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
	SuffixTree   *datastructure.SuffixTree // a single suffix tree containing all transcript sequences forward and reverse complemented
	NumSequences int
}

func (i GtaIndex) TransIndexToTransName(transcriptIndex uint32) string {
	return i.Transcripts[transcriptIndex].TranscriptIdEnsembl
}

func (i GtaIndex) SequenceIndexToTranscriptIndex(sequenceIndex uint32) uint32 {
	return sequenceIndex / 2
}

func (i GtaIndex) SequenceIndexIsForward(sequenceIndex uint32) bool {
	return sequenceIndex%2 == 0
}

// GetTranscriptSequenceDna returns the DNA sequence of the sequence at the given index.
// This returns the forward sequence if the index is even and the reverse sequence if the index is odd.
func (i GtaIndex) GetTranscriptSequenceDna(sequenceIndex int) string {
	return i.SuffixTree.Sequences[sequenceIndex]
}

// TranslateRwTransPosToFwTransPos translates a relative position revTransPos within a reverse
// transcript to a relative position within the forward transcript for the transcript at index transcriptIndex.
func (i GtaIndex) TranslateRwTransPosToFwTransPos(transcriptIndex uint32,
	revTransPos uint32) uint32 {

	transcript := i.Transcripts[transcriptIndex]

	if revTransPos >= uint32(transcript.SequenceLength) {
		panic("Relative position is out of bounds")
	}

	return uint32(transcript.SequenceLength) - revTransPos
}

// TranslateRelativeTranscriptPositionToRelativeGenePosition translates a relative position within a transcript to the
// relative position within the gene. Can be used to check whether two transcript positions relate to the same
// position on the gene.
// TODO: untested function
func (i GtaIndex) TranslateRelativeTranscriptPositionToRelativeGenePosition(transcriptIndex uint32,
	relativeTranscriptPosition uint32) uint32 {

	transcript := i.Transcripts[transcriptIndex]

	if relativeTranscriptPosition > uint32(transcript.SequenceLength) {
		panic("Relative position is out of bounds")
	}

	// Find the exon that contains the relative position by iterating over all exons and subtracting the
	// length of the exons that are before the relative position until the relative position is within the exon
	// then return the genomic location by adding the relative position to the genomic start location of the exon.
	for _, exon := range transcript.Exons {

		exonLength := exon.EndRelative - exon.StartRelative

		if relativeTranscriptPosition >= exonLength {
			relativeTranscriptPosition -= exonLength
			continue
		}

		return exon.StartRelative + relativeTranscriptPosition
	}

	// panic if the relative position could not be translated to a genomic location
	// should never happen because every relative position should be within the exons
	panic("Could not translate to genomic location")
}

// TransIntervalToGenomicIntervals translates a relative interval on a transcript to genomic intervals.
// The transcript interval is one continuous interval on the transcript sequence. This interval can span across
// multiple exons and can therefore consist of multiple genomic intervals.
// Idea:
// Iterate over the exons of the transcript and check if the interval intersects with the current exon.
// The length of each exon is subtracted from the relative transcript start position and until the relative start
// is larger than the length of the current exon, they do not overlap.
// When an exon overlaps the transcript interval, then the interval starts either at the start or within the exon.
// The interval starts at the start of the exon if the relative start position is larger than 0.
// The interval ends either at the end of the exon or within the exon.
// The interval ends at the end of the exon if the relative end position is larger than the length of the exon.
// Otherwise the offset (transcript start or transcript end) is added to the start of the exon to get the position
// within the exon.
// TODO: untested function
func (i GtaIndex) TransIntervalToGenomicIntervals(tranIndex uint32, startTransPos int, endTransPos int) []uint32 {

	transcript := i.Transcripts[tranIndex]

	//intervals := make([]*core.Interval, 0)
	intervals := make([]uint32, 0)

	for _, exon := range transcript.Exons {

		exonLength := exon.EndRelative - exon.StartRelative

		if startTransPos >= int(exonLength) {
			startTransPos -= int(exonLength)
			endTransPos -= int(exonLength)
			continue
		}

		var startGenomic uint32
		var endGenomic uint32

		if startTransPos >= 0 {
			startGenomic = exon.StartRelative + uint32(startTransPos)
			startTransPos -= int(exonLength)
		} else {
			startGenomic = exon.StartRelative
		}

		if endTransPos >= int(exonLength) {
			endGenomic = exon.EndRelative
			endTransPos -= int(exonLength)
		} else {
			endGenomic = exon.StartRelative + uint32(endTransPos)
		}

		//intervals = append(intervals, &core.Interval{
		//	Start: int(startGenomic),
		//	End:   int(endGenomic),
		//})
		intervals = append(intervals, startGenomic)
		intervals = append(intervals, endGenomic)
	}

	return intervals
}

// AddGenomicLocationsToNodes adds the genomic locations to the nodes of the suffix tree.
// Each node in the suffix tree has a list of positions which are the start positions of the suffix
// on each of the transcript sequences (forward and reverse). This function adds the genomic locations
// (the start positions of the suffix on the gene) to the positions.
func (i GtaIndex) AddGenomicLocationsToNodes() {
	for _, node := range i.SuffixTree.Nodes {
		for _, pos := range node.Positions {

			transcriptIndex := i.SequenceIndexToTranscriptIndex(uint32(pos.Index))

			transPos := uint32(pos.Start)
			if !i.SequenceIndexIsForward(uint32(pos.Index)) {
				transPos = i.TranslateRwTransPosToFwTransPos(transcriptIndex, uint32(pos.Start))
			}

			pos.StartGenomic = int(i.TranslateRelativeTranscriptPositionToRelativeGenePosition(transcriptIndex, transPos))
		}
	}
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

	// add the regular transcripts to the suffix tree and the index
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
	//gtaIndex.AddGenomicLocationsToNodes()
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
