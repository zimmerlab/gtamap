package index

import (
	"encoding/gob"
	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/datastructure"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/genemodel"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/keywordtree"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/keywordtreebyte"
	"github.com/KleinSamuel/gtamap/src/core/interval"
	"github.com/KleinSamuel/gtamap/src/dataloader"
	"github.com/KleinSamuel/gtamap/src/formats/gtf"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
	"log"
	"math"
	"os"
	"strconv"
	"strings"
	"time"
)

type GtaIndex struct {
	Gene               *genemodel.Gene           // information about the gene
	Transcripts        []*genemodel.Transcript   // all transcripts of the gene
	SuffixTree         *datastructure.SuffixTree // a single suffix tree containing all transcript sequences forward and reverse complemented
	NumSequences       int
	Sequences          []string
	KeywordTree        *keywordtree.KeywordTree
	EquivalenceClasses []*genemodel.EquivalenceClass
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

// GetSequenceByIndex returns the DNA sequence of the sequence at the given index.
// This returns the forward sequence if the index is even and the reverse sequence if the index is odd.
func (i GtaIndex) GetSequenceByIndex(sequenceIndex int) *string {
	if i.SequenceIndexIsForward(uint32(sequenceIndex)) {
		return &i.Transcripts[i.SequenceIndexToTranscriptIndex(uint32(sequenceIndex))].SequenceDnaForward
	} else {
		return &i.Transcripts[i.SequenceIndexToTranscriptIndex(uint32(sequenceIndex))].SequenceDnaReverse
	}
}

// TranslateRwTransPosToFwTransPos translates a relative position revTransPos within a reverse
// transcript to a relative position within the forward transcript for the transcript at index transcriptIndex.
func (i GtaIndex) TranslateRwTransPosToFwTransPos(transcriptIndex uint32,
	revTransPos uint32) uint32 {

	transcript := i.Transcripts[transcriptIndex]

	if revTransPos >= uint32(transcript.SequenceLength) {
		panic("Relative position is out of bounds")
	}

	return uint32(transcript.SequenceLength) - revTransPos - 1
}

// TranslateRelativeTranscriptPositionToRelativeGenePosition translates a relative position within a transcript to the
// relative position within the gene. Can be used to check whether two transcript positions relate to the same
// position on the gene.
// TODO: untested function
func (i GtaIndex) TranslateRelativeTranscriptPositionToRelativeGenePosition(transcriptIndex uint32,
	relativeTranscriptPosition uint32) uint32 {

	transcript := i.Transcripts[transcriptIndex]

	if relativeTranscriptPosition < 0 {
		panic("Relative position is out of bounds (less than 0)")
	}
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

func (i GtaIndex) FindEquivalenceClassByGenomicLocation(genomicLocation uint32) *genemodel.EquivalenceClass {
	for _, ec := range i.EquivalenceClasses {
		if genomicLocation >= ec.FromGenomic && genomicLocation < ec.ToGenomic {
			return ec
		}
	}
	return nil
}

// GetCoveredEquivalenceClassIds returns the ids of the equivalence classes that are covered by
// the given interval.
// The given startRelative and endRelative position is relative to the sequence (transcript) and not the gene.
// Idea:
// Each position of a transcript must be contained within an equivalence class because these are built upon
// the exons of all transcripts.
// The start position is translated to a genomic location and the equivalence class that contains this location
// is obtained as the initial and current equivalence class.
// While the given interval spans across the current equivalence class, the next equivalence class is retrieved.
// TODO: not tested
func (i GtaIndex) GetCoveredEquivalenceClassIds(sequenceIndex int, startRelative int, endRelative int) []uint32 {

	ecIds := make([]uint32, 0)

	transcriptIndex := i.SequenceIndexToTranscriptIndex(uint32(sequenceIndex))

	if !i.SequenceIndexIsForward(uint32(sequenceIndex)) {
		startRelative = int(i.TranslateRwTransPosToFwTransPos(transcriptIndex, uint32(startRelative)))
		endRelative = int(i.TranslateRwTransPosToFwTransPos(transcriptIndex, uint32(endRelative)))
	}

	start := i.TranslateRelativeTranscriptPositionToRelativeGenePosition(transcriptIndex, uint32(startRelative))

	ec := i.FindEquivalenceClassByGenomicLocation(start)

	if ec == nil {
		// should not happen because the equivalence class are based on the positions covered by exons
		panic("Could not find equivalence class for start position")
	}

	posRelative := startRelative + int(ec.FromGenomic-ec.ToGenomic)
	ecIds = append(ecIds, ec.Id)

	for posRelative < endRelative {

		// the given interval spans across the current equivalence class
		// find the next equivalence class
		ec := i.EquivalenceClasses[ec.Id+1]
		ecIds = append(ecIds, ec.Id)

		posRelative += int(ec.ToGenomic - ec.FromGenomic)
	}

	return ecIds
}

// AddSequenceToKeywordTree adds all keywords of a sequence to the keyword tree.
// The length of each keyword is defined by the keyword length of the keyword tree.
// Each keyword contained in the sequence is added to the keyword tree and the position of the keyword
// within the sequence is added to the keyword node as well as all equivalence classes that are covered by the
// keyword.
// TODO: not tested
func (i GtaIndex) AddSequenceToKeywordTree(sequence *string, sequenceIndex uint32) {

	for kStart := 0; kStart < len(*sequence)-int(i.KeywordTree.KeywordLength); kStart++ {

		kEnd := kStart + int(i.KeywordTree.KeywordLength)

		ecIds := i.GetCoveredEquivalenceClassIds(int(sequenceIndex), kStart, kEnd)

		node := i.KeywordTree.AddKeyword((*sequence)[kStart:kEnd])

		node.Positions = append(node.Positions, keywordtree.Position{
			SequenceIndex:       sequenceIndex,
			Position:            uint32(kStart),
			EquivalenceClassIds: ecIds,
		})

		//break
	}
}

func (i GtaIndex) EquivalenceClassIdsMatch(ecIds1 []uint32, ecIds2 []uint32) bool {

	if len(ecIds1) != len(ecIds2) {
		return false
	}

	for i, ecId := range ecIds1 {
		if ecId != ecIds2[i] {
			return false
		}
	}

	return true
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

	gtaIndex.Gene = &genemodel.Gene{
		GeneIdEnsembl:  annotation.Genes[0].GeneIdEnsembl,
		Chromosome:     annotation.Genes[0].Chromosome,
		IsFowardStrand: annotation.Genes[0].IsForwardStrand,
	}

	gtaIndex.Transcripts = make([]*genemodel.Transcript, len(annotation.Genes[0].Transcripts))

	// add the transcripts and its information of the first gene to the index
	for i, transcript := range annotation.Genes[0].Transcripts {

		gtaIndex.Transcripts[i] = &genemodel.Transcript{
			TranscriptIdEnsembl:       transcript.TranscriptIdEnsembl,
			SequenceDnaForward53Index: i * 2,
			SequenceDnaReverse53Index: i*2 + 1,
			SequenceDnaForward:        transcript.SequenceDna,
			SequenceDnaReverse:        utils.ReverseComplementDNA(transcript.SequenceDna),
			SequenceLength:            len(transcript.SequenceDna),
			Exons:                     make([]*genemodel.Exon, len(transcript.Exons)),
		}

		// the exons of the transcript are added to the transcript in the index
		for j, exon := range transcript.Exons {
			gtaIndex.Transcripts[i].Exons[j] = &genemodel.Exon{
				StartRelative: exon.StartRelative,
				EndRelative:   exon.EndRelative,
			}
		}
	}

	// generate equivalence classes
	exonIntervals := make([]*interval.Interval, 0)

	for _, transcript := range gtaIndex.Transcripts {
		for _, exon := range transcript.Exons {
			exonIntervals = append(exonIntervals, &interval.Interval{
				Start: int(exon.StartRelative),
				End:   int(exon.EndRelative),
			})
		}
	}

	equivalenceClasses := interval.GenerateEquivalenceClasses(exonIntervals)

	gtaIndex.EquivalenceClasses = make([]*genemodel.EquivalenceClass, len(equivalenceClasses))

	for i, ec := range equivalenceClasses {
		gtaIndex.EquivalenceClasses[i] = &genemodel.EquivalenceClass{
			Id:          uint32(i),
			FromGenomic: uint32(ec.Start),
			ToGenomic:   uint32(ec.End),
		}
	}

	//gtaIndex.SuffixTree = datastructure.CreateTree()

	gtaIndex.KeywordTree = keywordtree.NewKeywordTree(config.KmerLength())

	for _, transcript := range gtaIndex.Transcripts {

		timerBuildTree = time.Now()
		gtaIndex.AddSequenceToKeywordTree(&transcript.SequenceDnaForward, uint32(transcript.SequenceDnaForward53Index))
		durationBuildTree += time.Since(timerBuildTree)
		gtaIndex.NumSequences++

		timerBuildTree = time.Now()
		gtaIndex.AddSequenceToKeywordTree(&transcript.SequenceDnaReverse, uint32(transcript.SequenceDnaReverse53Index))
		durationBuildTree += time.Since(timerBuildTree)
		gtaIndex.NumSequences++

		//// the index of the forward sequence in the suffix tree (used for retrieval)
		//sequenceIndexForward := len(gtaIndex.SuffixTree.Sequences)
		//timerBuildTree = time.Now()
		//gtaIndex.SuffixTree.AddSequence(transcript.SequenceDna, sequenceIndexForward)
		////gtaIndex.SuffixTree.AddSequence(seq, sequenceIndexForward)
		//durationBuildTree += time.Since(timerBuildTree)
		//gtaIndex.NumSequences++
		//
		//// the index of the reverse complement sequence in the suffix tree (used for retrieval)
		//sequenceIndexReverse := len(gtaIndex.SuffixTree.Sequences)
		//timerBuildTree = time.Now()
		//gtaIndex.SuffixTree.AddSequence(utils.ReverseComplementDNA(transcript.SequenceDna), sequenceIndexReverse)
		////gtaIndex.SuffixTree.AddSequence(utils.ReverseComplementDNA(seq), sequenceIndexReverse)
		//durationBuildTree += time.Since(timerBuildTree)
		//gtaIndex.NumSequences++

		//break
	}

	// suffix links are only required for tree construction but not for the search
	//timerBuildTree = time.Now()
	//gtaIndex.SuffixTree.RemoveAllSuffixLinks()
	//durationBuildTree += time.Since(timerBuildTree)

	//posPerSeq := make(map[int][]int)
	//
	//for _, node := range gtaIndex.SuffixTree.Nodes {
	//	//fmt.Println("node id: ", node.Id)
	//
	//	for _, pos := range node.Positions {
	//		//fmt.Println("pos: ", pos)
	//		posPerSeq[pos.Index] = append(posPerSeq[pos.Index], pos.Start)
	//	}
	//}
	//
	//fmt.Println(len(gtaIndex.SuffixTree.Sequences))
	//
	//for i := 0; i < len(gtaIndex.SuffixTree.Sequences); i++ {
	//	fmt.Println("seqIndex:\t", i)
	//	fmt.Println("length:\t\t", len(gtaIndex.SuffixTree.Sequences[i]))
	//
	//	seen := make(map[int]bool) // Map to store seen integers
	//	var unique []int
	//
	//	for _, num := range posPerSeq[i] {
	//		if !seen[num] {
	//			seen[num] = true
	//			unique = append(unique, num)
	//		}
	//	}
	//
	//	sort.Ints(unique)
	//
	//	fmt.Println("length pos\t", len(unique))
	//
	//	for j := 0; j < len(unique); j++ {
	//		if unique[j] != j {
	//			fmt.Println("missing: ", j)
	//		}
	//	}
	//
	//	fmt.Println("")
	//}

	/*
		for seqIndex, positions := range posPerSeq {
			fmt.Println("seqIndex: ", seqIndex)
			fmt.Println("positions: ", positions)
		}
	*/

	// write the index to file as a serialized gob
	SerializeFromFile(&gtaIndex, outputFile)

	logrus.WithFields(logrus.Fields{
		"build": durationBuildTree.String(),
		"read":  durationReadReference.String(),
		"total": time.Since(timerStart).String(),
	}).Info("Successfully built and serialized the GTAMap index")
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

func SerializeFromFileTest(tree *keywordtree.KeywordTree, outputFile *os.File) {

	timerStart := time.Now()

	enc := gob.NewEncoder(outputFile)

	if err := enc.Encode(tree); err != nil {
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

/* GENOME INDEX */

type GenomeIndex struct {
	SequenceHeaders []string                     // the headers of the sequences at indices 0(+1), 2(+3), 4(+5), etc.
	SequenceInfo    []*gtf.GeneBasic             // information about the sequences 0(+1) 2(+3) ..
	Sequences       []*[]byte                    // the genome sequences (0: forward, 1: reverse complement, etc.)
	KeywordTree     *keywordtreebyte.KeywordTree // the keyword tree containing all kmers of both genome sequences
	KeywordMap      map[[10]byte][]*keywordtreebyte.Position
	KeywordMapSmall map[[5]byte][]*keywordtreebyte.Position
}

func (i GenomeIndex) AddKeywordToMap(keyword [10]byte, sequenceIndex uint8, position uint32) {
	if _, ok := i.KeywordMap[keyword]; !ok {
		i.KeywordMap[keyword] = make([]*keywordtreebyte.Position, 0)
	}
	i.KeywordMap[keyword] = append(i.KeywordMap[keyword], &keywordtreebyte.Position{
		SequenceIndex: sequenceIndex,
		Position:      position,
	})
}

func (i GenomeIndex) AddKeywordToMapSmall(keyword [5]byte, sequenceIndex uint8, position uint32) {
	if _, ok := i.KeywordMapSmall[keyword]; !ok {
		i.KeywordMapSmall[keyword] = make([]*keywordtreebyte.Position, 0)
	}
	i.KeywordMapSmall[keyword] = append(i.KeywordMapSmall[keyword], &keywordtreebyte.Position{
		SequenceIndex: sequenceIndex,
		Position:      position,
	})
}

func (i GenomeIndex) GetKeywordFromMap(keyword [10]byte) []*keywordtreebyte.Position {
	return i.KeywordMap[keyword]
}

func (i GenomeIndex) GetKeywordFromMapSmall(keyword [5]byte) []*keywordtreebyte.Position {
	return i.KeywordMapSmall[keyword]
}

func (i GenomeIndex) AddSequenceToKeywordTree(sequence *[]byte, sequenceIndex uint8) {

	for kStart := 0; kStart <= len(*sequence)-int(i.KeywordTree.KeywordLength); kStart++ {

		kEnd := kStart + int(i.KeywordTree.KeywordLength)

		node := i.KeywordTree.AddKeyword((*sequence)[kStart:kEnd])

		node.Positions = append(node.Positions, keywordtreebyte.Position{
			SequenceIndex: sequenceIndex,
			Position:      uint32(kStart),
		})

		seqBytes := *(*[10]byte)((*sequence)[kStart:kEnd])

		i.AddKeywordToMap(seqBytes, sequenceIndex, uint32(kStart))
	}
}

func (i GenomeIndex) AddSequenceToMapSmall(sequence *[]byte, sequenceIndex uint8) {

	for kStart := 0; kStart <= len(*sequence)-5; kStart++ {

		kEnd := kStart + 5

		seqBytes := *(*[5]byte)((*sequence)[kStart:kEnd])

		i.AddKeywordToMapSmall(seqBytes, sequenceIndex, uint32(kStart))
	}
}

func (i GenomeIndex) GetSequenceInfos() []sam.SequenceInfo {
	infos := make([]sam.SequenceInfo, len(i.SequenceInfo))

	for seqIndex, seqInfo := range i.SequenceInfo {
		infos[seqIndex] = sam.SequenceInfo{
			Name:   seqInfo.Contig,
			Length: int(seqInfo.EndGenomic - seqInfo.StartGenomic),
		}
	}

	return infos
}

func (i GenomeIndex) IsSequenceForward(sequenceIndex int) bool {
	return sequenceIndex%2 == 0
}

func (i GenomeIndex) GetRevCompIndex(sequenceIndex int) int {
	if i.IsSequenceForward(sequenceIndex) {
		return sequenceIndex + 1
	} else {
		return sequenceIndex - 1
	}
}

func (i GenomeIndex) GetSequenceHeader(sequenceIndex int) string {
	if i.IsSequenceForward(sequenceIndex) {
		return i.SequenceHeaders[sequenceIndex]
	} else {
		return i.SequenceHeaders[sequenceIndex-1]
	}
}

func (i GenomeIndex) GetSequenceInfo(sequenceIndex int) *gtf.GeneBasic {
	if i.IsSequenceForward(sequenceIndex) {
		return i.SequenceInfo[sequenceIndex]
	} else {
		return i.SequenceInfo[sequenceIndex-1]
	}
}

func (i GenomeIndex) GetSequenceContig(sequenceIndex int) string {
	if i.IsSequenceForward(sequenceIndex) {
		return i.SequenceInfo[sequenceIndex].Contig
	} else {
		return i.SequenceInfo[sequenceIndex-1].Contig
	}
}

func BuildGenomeIndex(fastaEntries []*dataloader.FastaEntry) *GenomeIndex {

	var index = GenomeIndex{
		SequenceHeaders: make([]string, len(fastaEntries)),
		SequenceInfo:    make([]*gtf.GeneBasic, len(fastaEntries)),
		Sequences:       make([]*[]byte, len(fastaEntries)*2),
		KeywordTree:     keywordtreebyte.NewKeywordTree(config.KmerLength()),
		KeywordMap:      make(map[[10]byte][]*keywordtreebyte.Position, int(math.Pow(4, 10))),
		KeywordMapSmall: make(map[[5]byte][]*keywordtreebyte.Position, int(math.Pow(4, 5))),
	}

	for i, entry := range fastaEntries {
		sequence := entry.Sequence
		sequenceRevComp := utils.ReverseComplementDnaBytes(sequence)

		info := parseFastaHeader(entry.Header)
		index.SequenceInfo[i] = info

		index.SequenceHeaders[i] = entry.Header

		index.Sequences[i*2] = &sequence
		index.Sequences[i*2+1] = &sequenceRevComp
	}

	for i, seq := range index.Sequences {
		index.AddSequenceToKeywordTree(seq, uint8(i))
		index.AddSequenceToMapSmall(seq, uint8(i))
	}

	index.KeywordTree.NumSequences = uint8(len(index.Sequences))

	return &index
}

func parseFastaHeader(header string) *gtf.GeneBasic {

	headerParts := strings.Split(header, "\t")

	if len(headerParts) < 4 {
		logrus.Fatal("invalid header:", header)
	}

	geneId := headerParts[0]
	contig := headerParts[1]
	isForwardStrand := headerParts[2] == "+"
	startGenomicInt, err := strconv.Atoi(headerParts[3])
	if err != nil {
		logrus.Fatal("invalid start genomic position:", header)
	}
	startGenomic := uint32(startGenomicInt)
	endGenomicInt, err := strconv.Atoi(headerParts[4])
	if err != nil {
		logrus.Fatal("invalid end genomic position:", header)
	}
	endGenomic := uint32(endGenomicInt)

	return &gtf.GeneBasic{
		GeneId:          geneId,
		Contig:          contig,
		IsForwardStrand: isForwardStrand,
		StartGenomic:    startGenomic,
		EndGenomic:      endGenomic,
	}
}

func WriteGenomeIndex(genomeIndex *GenomeIndex, outputFile *os.File) {

	timerStart := time.Now()

	enc := gob.NewEncoder(outputFile)

	if err := enc.Encode(genomeIndex); err != nil {
		log.Fatal("encode error:", err)
	}

	logrus.WithFields(logrus.Fields{
		"duration": time.Since(timerStart),
		"output":   outputFile.Name(),
	}).Info("Serialized index")
}

func ReadGenomeIndexByPath(indexFilePath string) *GenomeIndex {
	// Open the file for reading
	file, err := os.Open(indexFilePath)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	return ReadGenomeIndexByFile(file)
}

func ReadGenomeIndexByFile(indexFile *os.File) *GenomeIndex {

	timerStart := time.Now()

	// create a decoder and deserialize the person struct from the file
	decoder := gob.NewDecoder(indexFile)

	var genomeIndex GenomeIndex

	if err := decoder.Decode(&genomeIndex); err != nil {
		panic(err)
	}

	logrus.WithFields(logrus.Fields{
		"duration": time.Since(timerStart),
	}).Info("Deserialized index")

	return &genomeIndex
}

func BuildAndSerializeGenomeIndex(fastaFile *os.File, outputFile *os.File) {

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

	genomeIndex := BuildGenomeIndex(fastaEntries)

	WriteGenomeIndex(genomeIndex, outputFile)
}
