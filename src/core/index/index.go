package index

import (
	"encoding/gob"
	"fmt"
	"io/fs"
	"log"
	"math"
	"os"
	"path/filepath"
	"slices"
	"strconv"
	"strings"
	"time"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/datastructure"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/keywordtreebyte"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/dataloader"
	"github.com/KleinSamuel/gtamap/src/formats/gtf"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
)

type GenomeIndex struct {
	SequenceHeaders []string                     // the headers of the sequences at indices 0(+1), 2(+3), 4(+5), etc.
	SequenceInfo    []*gtf.GeneBasic             // information about the sequences 0(+1) 2(+3) ..
	Sequences       []*[]byte                    // the genome sequences (0: forward, 1: reverse complement, etc.)
	KeywordTree     *keywordtreebyte.KeywordTree // the keyword tree containing all kmers of both genome sequences
	KeywordMap      map[[10]byte][]*keywordtreebyte.Position
	KeywordMapSmall map[[5]byte][]*keywordtreebyte.Position
	// ParalogRegions  map[string]*GenomeIndex // additional index per target region of paralog regions
	// Blacklist        map[int]*datastructure.INTree  // interval tree used to store repeat regions. Tree at 0 -> for target region 0 etc...
	// RepeatRegions    map[int][]datastructure.Bounds // interval tree used to store repeat regions. Tree at 0 -> for target region 0 etc...
	// ContigRepeatmask      map[string]*ContigRepeatmask
	RegionMask            *RegionMask
	ContigToTargetRegions map[string]*regionvector.RegionVector // map to get the target region for a given contig
}

func (i *GenomeIndex) SetRegionMask(mask *RegionMask) {
	i.RegionMask = mask
}

func (i *GenomeIndex) AddKeywordToMap(keyword [10]byte, sequenceIndex uint8, position uint32) {
	if _, ok := i.KeywordMap[keyword]; !ok {
		i.KeywordMap[keyword] = make([]*keywordtreebyte.Position, 0)
	}
	i.KeywordMap[keyword] = append(i.KeywordMap[keyword], &keywordtreebyte.Position{
		SequenceIndex: sequenceIndex,
		Position:      position,
	})
}

// func (i *GenomeIndex) LoadParalogs(paralogFile *os.File) {
// 	// init map to store paralog regions per target sequence
// 	// let's say we have target sequence X and Y
// 	// then i.ParalogRegions will have X and Y as key and store
// 	// one paralogIndex per target seq
// 	// i.ParalogRegions[X] -> *genomeIndex (which will contain all kmers of paralog seqs of X in one map)
// 	i.ParalogRegions = make(map[string]*GenomeIndex)
// 	scanner := bufio.NewScanner(paralogFile)
// 	for scanner.Scan() {
// 		line := scanner.Text()
// 		contentParts := strings.Split(line, ",")
// 		if len(contentParts) != 2 {
// 			logrus.Fatal("Error parsing paralog file. Wrong format!")
// 		}
// 		targetRegion := contentParts[0]
// 		pathToParalogIndex := contentParts[1]
// 		paralogIndex := ReadGenomeIndexByPath(pathToParalogIndex)
// 		i.AddParalogRegionIndex(targetRegion, paralogIndex)
// 	}
// }

// func (i *GenomeIndex) AddParalogRegionIndex(seqId string, index *GenomeIndex) {
// 	_, exists := i.ParalogRegions[seqId]
// 	if !exists {
// 		i.ParalogRegions[seqId] = index
// 	} else {
// 		i.ParalogRegions[seqId].ExtendParalogRegionIndex(seqId, index)
// 	}
// }

func (i *GenomeIndex) ExtendParalogRegionIndex(targetGene string, indexExtension *GenomeIndex) {
	// since the indices we add to our main index i are mainly paralog indices with only
	// one sequence, we panic if there are more than 1 sequenceHeader in indexExtension.
	// indexExtension should look like this
	// Sequences -> fw rv
	// SequenceHeaders -> one id
	// SequenceInfo -> one *gtf.GeneBasic
	// Keyword(Map/Tree)(Small) -> one each
	if len(indexExtension.SequenceHeaders) > 1 {
		logrus.Fatal("Error extending main index with paralog region: paralog region index extension contains more than one sequence: ", indexExtension.SequenceHeaders)
	}
	// append seq header
	i.SequenceHeaders = append(i.SequenceHeaders, indexExtension.SequenceHeaders[0])

	// append fw and rv seq
	for _, sequence := range indexExtension.Sequences {
		i.Sequences = append(i.Sequences, sequence)
	}

	// append gtf.GeneBasic info
	i.SequenceInfo = append(i.SequenceInfo, indexExtension.SequenceInfo[0])

	// append Maps; here we need to make sure that the pos.SequenceIndex matches the
	// current extension. Since the extension is only allowed to hold one
	// sequence, the pos.SequenceIndex is in {0, 1}. By adding an offset based on how many
	// seqs we already added, the kmer positions reference the correct sequenceIndex
	// if we add a seq to an index which holds only one seq:
	// 0, 1 -> 2, 3
	// if we then add another seq
	// 0, 1 -> 4, 5

	// how many sequences does the subindex currently hold (after adding the new seqId)
	numSequences := uint8(len(i.SequenceHeaders))
	for kmer, positions := range indexExtension.KeywordMap {
		for _, pos := range positions {
			var offset uint8
			if pos.SequenceIndex == 0 {
				offset = (numSequences - 1) * 2
			} else {
				offset = (numSequences-1)*2 + 1
			}
			i.AddKeywordToMap(kmer, offset, pos.Position)
		}
	}

	// do the same for the small keyword map
	for kmer, positions := range indexExtension.KeywordMapSmall {
		for _, pos := range positions {
			var offset uint8
			if pos.SequenceIndex == 0 {
				offset = (numSequences - 1) * 2
			} else {
				offset = (numSequences-1)*2 + 1
			}
			i.AddKeywordToMapSmall(kmer, offset, pos.Position)
		}
	}

	// log
	logrus.WithFields(logrus.Fields{
		"Added paralog region":                            indexExtension.SequenceInfo[0].GeneId,
		"to the paralog index belonging to target region": targetGene,
	}).Debug("Extended paralog index of main index")
}

func (i *GenomeIndex) AddKeywordToMapSmall(keyword [5]byte, sequenceIndex uint8, position uint32) {
	if _, ok := i.KeywordMapSmall[keyword]; !ok {
		i.KeywordMapSmall[keyword] = make([]*keywordtreebyte.Position, 0)
	}
	i.KeywordMapSmall[keyword] = append(i.KeywordMapSmall[keyword], &keywordtreebyte.Position{
		SequenceIndex: sequenceIndex,
		Position:      position,
	})
}

func (i *GenomeIndex) GetKeywordFromMap(keyword [10]byte) []*keywordtreebyte.Position {
	return i.KeywordMap[keyword]
}

func (i *GenomeIndex) FindKeywordMatchesInMap(keyword *[10]byte, posInRead int) []*mapperutils.Match {
	positions := i.KeywordMap[*keyword]

	matches := make([]*mapperutils.Match, len(positions))

	for i, pos := range positions {
		matches[i] = &mapperutils.Match{
			SequenceIndex: int(pos.SequenceIndex),
			FromGenome:    int(pos.Position),
			ToGenome:      int(pos.Position + 10),
			FromRead:      posInRead,
			ToRead:        posInRead + 10,
			StartGenome:   int(pos.Position) - posInRead,
		}
	}

	return matches
}

func (i *GenomeIndex) GetKeywordFromMapSmall(keyword [5]byte) []*keywordtreebyte.Position {
	return i.KeywordMapSmall[keyword]
}

func (i *GenomeIndex) AddSequenceToKeywordTree(sequence *[]byte, sequenceIndex uint8) {
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

func (i *GenomeIndex) AddSequenceToMap(sequence *[]byte, sequenceIndex uint8) {
	for kStart := 0; kStart <= len(*sequence)-int(config.KmerLength()); kStart++ {

		kEnd := kStart + int(config.KmerLength())

		seqBytes := *(*[10]byte)((*sequence)[kStart:kEnd])

		i.AddKeywordToMap(seqBytes, sequenceIndex, uint32(kStart))
	}
}

func (i *GenomeIndex) AddSequenceToMapSmall(sequence *[]byte, sequenceIndex uint8) {
	for kStart := 0; kStart <= len(*sequence)-5; kStart++ {

		kEnd := kStart + 5

		seqBytes := *(*[5]byte)((*sequence)[kStart:kEnd])

		i.AddKeywordToMapSmall(seqBytes, sequenceIndex, uint32(kStart))
	}
}

func (i *GenomeIndex) GetSequenceInfos() []sam.SequenceInfo {
	infos := make([]sam.SequenceInfo, len(i.SequenceInfo))

	for seqIndex, seqInfo := range i.SequenceInfo {
		infos[seqIndex] = sam.SequenceInfo{
			Name: seqInfo.Contig,
			// Length: int(seqInfo.EndGenomic - seqInfo.StartGenomic),
			Length: int(seqInfo.EndGenomic),
		}
	}

	return infos
}

// IsSequenceForward checks if the sequence at the given index is forward or reverse complemented.
// Returns true if the given sequence index refers to a 5->3 forward sequence.
func (i *GenomeIndex) IsSequenceForward(sequenceIndex int) bool {
	return sequenceIndex%2 == 0
}

// GetRevCompIndex returns the index of the reverse complement sequence for the given sequence index.
func (i *GenomeIndex) GetRevCompIndex(sequenceIndex int) int {
	if i.IsSequenceForward(sequenceIndex) {
		return sequenceIndex + 1
	} else {
		return sequenceIndex - 1
	}
}

// GetForwardSequence returns the forward sequence for the given sequence index.
func (i *GenomeIndex) GetForwardSequence(sequenceIndex int) *[]byte {
	if i.IsSequenceForward(sequenceIndex) {
		return i.Sequences[sequenceIndex]
	} else {
		return i.Sequences[sequenceIndex-1]
	}
}

func (i *GenomeIndex) GetSequenceHeader(sequenceIndex int) string {
	if i.IsSequenceForward(sequenceIndex) {
		return i.SequenceHeaders[sequenceIndex]
	} else {
		return i.SequenceHeaders[sequenceIndex-1]
	}
}

func (i *GenomeIndex) GetSequenceInfo(sequenceIndex int) *gtf.GeneBasic {
	return i.SequenceInfo[sequenceIndex/2]
}

func (i *GenomeIndex) GetSequenceContig(sequenceIndex int) string {
	// if i.IsSequenceForward(sequenceIndex) {
	// 	return i.SequenceInfo[sequenceIndex].Contig
	// } else {
	// 	return i.SequenceInfo[sequenceIndex-1].Contig
	// }
	return i.SequenceInfo[sequenceIndex/2].Contig
}

func (i *GenomeIndex) NumSequences() int {
	return len(i.SequenceHeaders)
}

func (i *GenomeIndex) IsResultValid(result *mapperutils.ReadMatchResult) bool {
	sequenceInfo := i.GetSequenceInfo(result.SequenceIndex)

	if _, found := i.RegionMask.ContigMasks[sequenceInfo.Contig]; !found {
		return true
	}

	length := int(sequenceInfo.EndGenomic - sequenceInfo.StartGenomic)
	isForward := i.IsSequenceForward(result.SequenceIndex)

	fmt.Println("--")
	// DEBUG: print result
	fmt.Println(result.String())

	result.MergeRegionsIfBothConsecutive()

	fmt.Println(result.String())

	for _, genomicRegion := range result.MatchedGenome.Regions {

		var startGenomic int
		if isForward {
			startGenomic = int(sequenceInfo.StartGenomic) + genomicRegion.Start
		} else {
			startGenomic = int(sequenceInfo.StartGenomic) + (length - genomicRegion.End)
		}

		fmt.Printf("region:\nrel:\t%d - %d\nglb:\t%d - %d\n\n", genomicRegion.Start, genomicRegion.End, startGenomic, startGenomic+genomicRegion.Length())

		i.RegionMask.ContigMasks[sequenceInfo.Contig].ApplyMaskToRegion(startGenomic, startGenomic+genomicRegion.Length())

		//
		// 	found, name, size := i.RegionMask.ContigMasks[sequenceInfo.Contig].
		// 		GetMostImportantItemThatOverlaps(
		// 			startGenomic,
		// 			startGenomic+genomicRegion.Length(),
		// 		)
		//
		// 	if !found {
		// 		continue
		// 	}
		//
		// 	// percentOverlap := float32(size) / float32(genomicRegion.Length())
		// 	fmt.Println()
		// 	fmt.Println(genomicRegion)
		// 	fmt.Println(startGenomic, startGenomic+genomicRegion.Length())
		// 	fmt.Println(name, size)
		// 	fmt.Println(float32(size) / float32(genomicRegion.Length()))
		//
	}

	return true
}

// TODO: adjust to region mask
// func (i *GenomeIndex) IsPartOfRepeat(result *mapperutils.ReadMatchResult) bool {
//
// 	seqInfo := i.GetSequenceInfo(result.SequenceIndex)
//
// 	// no repeatmask for this contig
// 	if _, found := i.ContigRepeatmask[seqInfo.Contig]; !found {
// 		return false
// 	}
//
// 	length := int(seqInfo.EndGenomic - seqInfo.StartGenomic)
// 	isForward := i.IsSequenceForward(result.SequenceIndex)
//
// 	for _, genomicRegion := range result.MatchedGenome.Regions {
//
// 		var startGenomic int
// 		if isForward {
// 			startGenomic = int(seqInfo.StartGenomic) + genomicRegion.Start
// 		} else {
// 			startGenomic = int(seqInfo.StartGenomic) + (length - genomicRegion.End)
// 		}
//
// 		repeats := i.ContigRepeatmask[seqInfo.Contig].TreeFw.Including(startGenomic)
//
// 		if len(repeats) != 0 && genomicRegion.Length() > 70 {
// 			return true
// 		}
// 	}
//
// 	return false
// }

// CleanResults removes results that do not match contraints.
// Removes results that are part of a repeat and have more than the allowed
// number of mismatches.
func (i *GenomeIndex) CleanResults(
	results []*mapperutils.ReadMatchResult,
) []*mapperutils.ReadMatchResult {
	cleaned := make([]*mapperutils.ReadMatchResult, 0)

	for _, result := range results {
		// skip results that are part of a repeat region and have more than 4 mismatches
		// TODO: adjust to region mask
		// if i.IsPartOfRepeat(result) && len(result.MismatchesRead) > 4 {
		// 	continue
		// }
		cleaned = append(cleaned, result)
	}

	return cleaned
}

func BuildGenomeIndex(fastaEntries []*dataloader.FastaEntry) *GenomeIndex {
	timerStart := time.Now()

	index := GenomeIndex{
		SequenceHeaders: make([]string, len(fastaEntries)),
		SequenceInfo:    make([]*gtf.GeneBasic, len(fastaEntries)),
		Sequences:       make([]*[]byte, len(fastaEntries)*2),
		KeywordMap:      make(map[[10]byte][]*keywordtreebyte.Position, int(math.Pow(4, 10))),
		KeywordMapSmall: make(map[[5]byte][]*keywordtreebyte.Position, int(math.Pow(4, 5))),
		// Blacklist:        make(map[int]*datastructure.INTree),
		// RepeatRegions:    make(map[int][]datastructure.Bounds),
		// ContigRepeatmask:      make(map[string]*ContigRepeatmask),
		ContigToTargetRegions: make(map[string]*regionvector.RegionVector),
	}

	// key = contig, value = list of bounds (start+end of targets)
	// used to build interval trees per contig to speed up repeatmasker loading
	// contigMapTargets :=

	for i, entry := range fastaEntries {

		sequence := entry.Sequence
		sequenceRevComp, revCompErr := utils.ReverseComplementDnaBytes(sequence)
		if revCompErr != nil {
			logrus.Fatal("Error reversing and complementing sequence", revCompErr)
		}

		info := parseFastaHeader(entry.Header)
		index.SequenceInfo[i] = info
		index.SequenceHeaders[i] = entry.Header
		index.Sequences[i*2] = &sequence
		index.Sequences[i*2+1] = &sequenceRevComp

		// // populate contig to target region map
		// if _, exists := index.ContigToTargetRegions[info.Contig]; !exists {
		// 	index.ContigToTargetRegions[info.Contig] = regionvector.NewRegionVector()
		// }
		// index.ContigToTargetRegions[info.Contig].AddRegionAndMerge(
		// 	int(info.StartGenomic),
		// 	int(info.EndGenomic),
		// )
	}

	for i, seq := range index.Sequences {
		index.AddSequenceToMap(seq, uint8(i))
		// TODO: do we plan to use this? if not then remove
		index.AddSequenceToMapSmall(seq, uint8(i))
	}

	logrus.WithFields(logrus.Fields{
		"duration": utils.FormatDuration(time.Since(timerStart)),
	}).Info("Built index")

	return &index
}

func parseFastaHeader(header string) *gtf.GeneBasic {
	headerParts := strings.Split(header, "\t")

	if len(headerParts) < 4 {
		logrus.Fatal("invalid header:", header)
	}

	contig := headerParts[0]
	geneId := headerParts[1]
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
	gob.Register(&datastructure.RepeatRegion{})

	if err := enc.Encode(genomeIndex); err != nil {
		log.Fatal("encode error:", err)
	}

	logrus.WithFields(logrus.Fields{
		"duration": utils.FormatDuration(time.Since(timerStart)),
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
	gob.Register(&datastructure.RepeatRegion{})

	var genomeIndex GenomeIndex

	if err := decoder.Decode(&genomeIndex); err != nil {
		panic(err)
	}

	logrus.WithFields(logrus.Fields{
		"duration": utils.FormatDuration(time.Since(timerStart)),
	}).Info("Deserialized index")

	return &genomeIndex
}

func BuildAndSerializeGenomeIndex(
	fastaFile *os.File,
	outputFile *os.File,
	regionmaskFile *os.File,
) {
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

	// add region mask to index if provided
	AddRegionmaskToIndex(regionmaskFile, genomeIndex)

	fmt.Println(genomeIndex.RegionMask)

	WriteGenomeIndex(genomeIndex, outputFile)
}

func AddRegionmaskToIndex(
	regionmaskFile *os.File,
	genomeIndex *GenomeIndex,
) {
	if regionmaskFile != nil {
		logrus.WithFields(logrus.Fields{
			"region mask": regionmaskFile.Name(),
		}).Info("Using region mask bed file and priority file")
	} else {
		logrus.Info("No region mask used")
	}

	mask, errMask := NewRegionMask(
		regionmaskFile,
		// genomeIndex.ContigToTargetRegions,
		genomeIndex.SequenceInfo,
	)

	if errMask != nil {
		logrus.Fatal("Error loading region mask", errMask)
	}

	genomeIndex.SetRegionMask(mask)
}

func BuildBlacklistInTree(start uint32, end uint32, contig string, blackListFile *os.File) (*datastructure.INTree, *datastructure.INTree, []datastructure.Bounds, []datastructure.Bounds, error) {
	repeatRegions := make([]datastructure.Bounds, 0)
	repeatRegionsRev := make([]datastructure.Bounds, 0)
	regionLength := end - start

	scanner, err := dataloader.OpenFileScanner(blackListFile)
	if err != nil {
		logrus.Fatal(err)
	}

	for scanner.Scan() {
		line := scanner.Text()
		chr := extractFieldFromBlacklist(4, line)
		if strings.Contains(chr, "chr") {
			chr = strings.Trim(chr, "chr")
		}
		if chr != contig {
			continue
		}

		regionStart, err := strconv.ParseUint(extractFieldFromBlacklist(5, line), 10, 32)
		if err != nil {
			logrus.Fatal(err)
		}
		regionEnd, err := strconv.ParseUint(extractFieldFromBlacklist(6, line), 10, 32)
		if err != nil {
			logrus.Fatal(err)
		}
		// append relative coords to gene
		regionStartConv := uint32(regionStart)
		regionEndConv := uint32(regionEnd)
		if regionStartConv >= start && regionEndConv <= end {
			repeatRegions = append(repeatRegions, &datastructure.RepeatRegion{Lower: int(regionStartConv - start), Upper: int(regionEndConv - start)})
			repeatRegionsRev = append(repeatRegionsRev, &datastructure.RepeatRegion{Lower: int(regionLength - (regionEndConv - start)), Upper: int(regionLength - (regionStartConv - start))})
		}
		if regionStartConv > end && regionEndConv > end {
			// TODO: should repeatRegionsRev be reversed here too?
			return datastructure.NewINTree(repeatRegions), datastructure.NewINTree(repeatRegionsRev), repeatRegions, repeatRegionsRev, nil
		}

	}

	// TODO: why is this reversed here but not in early termination case above?
	slices.Reverse(repeatRegionsRev)
	return datastructure.NewINTree(repeatRegions), datastructure.NewINTree(repeatRegionsRev), repeatRegions, repeatRegionsRev, nil
}

type ContigRepeatmask struct {
	Contig    string
	TreeFw    *datastructure.INTree
	TreeRv    *datastructure.INTree
	RegionsFw []datastructure.Bounds
	RegionsRv []datastructure.Bounds
}

func BuildContigRepeatmaskTrees(repeatmaskFile *os.File, targetTrees map[string]*datastructure.INTree) map[string]*ContigRepeatmask {
	scanner, err := dataloader.OpenFileScanner(repeatmaskFile)
	if err != nil {
		logrus.Fatal(err)
	}

	contigMap := make(map[string]*ContigRepeatmask)

	for scanner.Scan() {

		line := scanner.Text()
		lineParts := strings.Fields(line)

		if len(lineParts) < 7 {
			continue
		}

		chr := lineParts[4]
		if strings.Contains(chr, "chr") {
			chr = strings.Trim(chr, "chr")
		}

		if targetTrees != nil {
			if _, exists := targetTrees[chr]; !exists {
				continue
			}
		}

		regionStart, errStart := strconv.ParseUint(lineParts[5], 10, 32)
		if errStart != nil {
			// logrus.Fatal(errStart)
			continue
		}
		regionEnd, errEnd := strconv.ParseUint(lineParts[6], 10, 32)
		if errEnd != nil {
			// logrus.Fatal(errEnd)
			continue
		}

		regionStartInt := int(regionStart)
		regionEndInt := int(regionEnd)

		targetTree := targetTrees[chr]

		if len(targetTree.Including(regionStartInt)) == 0 {
			continue
		}

		if _, exists := contigMap[chr]; !exists {
			contigMap[chr] = &ContigRepeatmask{
				Contig:    chr,
				RegionsFw: make([]datastructure.Bounds, 0),
			}
		}

		contigMap[chr].RegionsFw = append(contigMap[chr].RegionsFw,
			&datastructure.RepeatRegion{
				Lower: regionStartInt,
				Upper: regionEndInt,
			})
	}

	for _, contigRepeatmask := range contigMap {
		contigRepeatmask.TreeFw = datastructure.NewINTree(contigRepeatmask.RegionsFw)
	}

	return contigMap
}

func extractFieldFromBlacklist(field int, line string) string {
	parts := strings.Fields(line)
	if field >= 0 && field < len(parts) {
		return parts[field]
	}
	return ""
}

func OptimizeFastaExtraction(targetParalogs map[string]map[string]struct{}, fastaDirParalogPre *string) map[string]map[string]struct{} {
	paralogsToExtractSeq := make(map[string]map[string]struct{}, 0)
	foundFiles := 0
	existingFaFiles := make(map[string]struct{}, 0)

	for target, paralogs := range targetParalogs {
		logrus.Infof("Found %s paralog region(s) of target region '%s' in DB.", strconv.Itoa(len(paralogs)), target)
	}

	err := filepath.WalkDir(*fastaDirParalogPre, func(path string, d fs.DirEntry, err error) error {
		if err != nil {
			return err
		}
		if !d.IsDir() {
			faFileName := strings.Split(d.Name(), ".fa")[0]
			var target string
			for t, paralogs := range targetParalogs {
				_, exists := paralogs[faFileName]
				if exists {
					target = t
					logrus.Infof("Found paralog region '%s' of target region '%s' in --fastaout '%s'. No sequence extraction necessary.", faFileName, target, *fastaDirParalogPre)
				}
			}
			existingFaFiles[faFileName] = struct{}{}
			foundFiles++
		}
		return nil
	})
	if err != nil {
		logrus.Fatalf("Error reading %s directory to check for pre-computed fa seqs: %s", *fastaDirParalogPre, err)
	}

	// if no file exists in --fastaout then all seqs need to be extracted
	if foundFiles == 0 {
		logrus.Infof("All paralog fasta sequences of all target regions need to be extracted.")
	}

	for target, paralogs := range targetParalogs {
		paralogsToExtractSeq[target] = make(map[string]struct{})
		for paralog := range paralogs {
			_, exists := existingFaFiles[paralog]
			if !exists {
				logrus.Infof("Will have to extract paralog region '%s' of target region '%s' since .fa non-existent in --fastaout '%s'.", paralog, target, *fastaDirParalogPre)
				paralogsToExtractSeq[target][paralog] = struct{}{}
			}
		}
	}
	return paralogsToExtractSeq
}

func OptimizeIndexSerialisation(targetParalogs map[string]map[string]struct{}, indexDirParalogPre *string) map[string]map[string]struct{} {
	foundIndices := 0
	paralogsToSerialize := make(map[string]map[string]struct{})
	existingIndices := make(map[string]struct{}, 0)
	err := filepath.WalkDir(*indexDirParalogPre, func(path string, d fs.DirEntry, err error) error {
		if err != nil {
			return err
		}

		if !d.IsDir() {
			gtaiName := strings.Split(d.Name(), ".gtai")[0]
			var target string
			for t, paralogs := range targetParalogs {
				_, exists := paralogs[gtaiName]
				if exists {
					target = t
					logrus.Infof("Found serialized index '%s'.gai of target region '%s' in --indexout '%s'.", gtaiName, target, *indexDirParalogPre)
				}
			}
			existingIndices[gtaiName] = struct{}{}
			foundIndices++
		}
		return nil
	})
	if err != nil {
		logrus.Fatalf("Error reading %s directory to check for already existing indices: %s", *indexDirParalogPre, err)
	}

	// if no file exists in --fastaout then all seqs need to be extracted
	if foundIndices == 0 {
		logrus.Infof("All paralog sequences of all target regions need to be serialized to index.")
	}

	for target, paralogs := range targetParalogs {
		paralogsToSerialize[target] = make(map[string]struct{})
		for paralog := range paralogs {
			_, exists := existingIndices[paralog]
			if !exists {
				logrus.Infof("Will have to serialize paralog region '%s' of target region '%s' since .gtai non-existent in --indexout '%s'.", paralog, target, *indexDirParalogPre)
				paralogsToSerialize[target][paralog] = struct{}{}
			}
		}
	}
	return paralogsToSerialize
}

// OLD PARALOG LOGIC (NOT NEEDED ANYMORE)
// func BuildAndSerializeAll(paralogsToSerialize map[string]map[string]struct{}, indexDirParalogPre *string, fastaDirParalogPre *string) {
// 	for target, paralogs := range paralogsToSerialize {
// 		for paralog := range paralogs {
// 			// create and format index file
// 			paralogIndexName := fmt.Sprintf("%s.gtai", paralog)
// 			logrus.Infof("Creating index for paralog %s of target region %s in: %s ", paralog, target, paralogIndexName)
// 			paralogIndexPath := filepath.Join(*indexDirParalogPre, paralogIndexName)
// 			paralogIndexFile, err := os.Create(paralogIndexPath)
// 			if err != nil {
// 				logrus.Fatalf("Error creating paralog index file: %s: %s", paralogIndexPath, err)
// 			}
//
// 			// open fa file
// 			paralogFastaName := fmt.Sprintf("%s.fa", paralog)
// 			paralogFastaPath := filepath.Join(*fastaDirParalogPre, paralogFastaName)
// 			paralogFastaFile, err := os.Open(paralogFastaPath)
// 			if err != nil {
// 				fmt.Println(err)
// 				logrus.Fatalf("Error reading paralog fa file: %s", paralogFastaPath)
// 			}
//
// 			BuildAndSerializeGenomeIndex(paralogFastaFile, paralogIndexFile)
// 		}
// 	}
// }
