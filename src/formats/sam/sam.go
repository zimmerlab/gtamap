package sam

import (
	"bufio"
	"errors"
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"io"
	"os"
	"strconv"
)

const Version string = "1.6"

type Header struct {
	Version       string         // version of the SAM format (e.g. 1.6)
	SequenceInfos []SequenceInfo // information about the contained sequences
	ToolVersion   string
}

type SequenceInfo struct {
	Name   string
	Length int
}

type TranscriptInfo struct {
	Id                  string
	TranscriptEnsemblId string
	TranscriptLength    int
}

func (header *Header) String() string {
	headerString := fmt.Sprintf("@HD\tVN:%s\n", header.Version)

	for _, seqInfo := range header.SequenceInfos {
		headerString += fmt.Sprintf("@SQ\tSN:%s\tLN:%d\n", seqInfo.Name, seqInfo.Length)
	}

	headerString += fmt.Sprintf("@PG\tID:%s\tPN:%s\tVN:%s\n", "GTAMap", "GTAMap", header.ToolVersion)

	return headerString
}

type RecordPair struct {
	First  *Record
	Second *Record
}

type Record struct {
	Qname         string // id of the read
	Flag          Flag   // bitwise flag
	Rname         string // id of the reference
	Pos           int    // 1-based leftmost mapping position
	Mapq          int    // mapping quality
	Cigar         string // CIGAR string (compact representation of alignment)
	Rnext         string // id of the mate reference
	Pnext         int    // position of the mate reference
	Tlen          int    // observed template length
	Seq           string // read sequence
	Qual          string // read quality
	TranscriptIds []int  // ids of the transcripts that this read was aligned to
}

func (entry *Record) String() string {

	if entry == nil {
		return ""
	}

	alignmentLineString := fmt.Sprintf("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s", entry.Qname, entry.Flag.String(), entry.Rname, entry.Pos, entry.Mapq, entry.Cigar, entry.Rnext, entry.Pnext, entry.Tlen, entry.Seq, entry.Qual)
	//alignmentLineString += fmt.Sprintf("\tXT:Z:T%d", entry.TranscriptId)

	for i, transcriptId := range entry.TranscriptIds {
		if i == 0 {
			alignmentLineString += fmt.Sprintf("\tXT:Z:T%d", transcriptId)
		} else {
			alignmentLineString += fmt.Sprintf(",T%d", transcriptId)
		}
	}

	return alignmentLineString
}

func WriteSamRecordLine(qname string, flag Flag, rname string, pos int, mapq int, cigar string, rnext string, pnext int, tlen int, seq string, qual string, transcriptIds []int) string {
	alignmentLineString := fmt.Sprintf("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s", qname, flag.String(), rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual)

	for i, transcriptId := range transcriptIds {
		if i == 0 {
			alignmentLineString += fmt.Sprintf("\tXT:Z:T%d", transcriptId)
		} else {
			alignmentLineString += fmt.Sprintf(",T%d", transcriptId)
		}
	}
	return alignmentLineString
}

// Flag is a struct that represents the bitwise flag of a SAM entry
type Flag struct {
	Value int
}

// String returns the string representation of the flag (the integer value)
func (flag *Flag) String() string {
	return fmt.Sprintf("%d", flag.Value)
}

// SetPaired sets the read as paired which means that the read is part of a pair and has a mate
func (flag *Flag) SetPaired() {
	flag.Value |= 1 << 0
}
func (flag *Flag) IsPaired() bool {
	return (flag.Value & (1 << 0)) != 0
}

// SetProperlyPaired sets the read as properly paired
// A read is properly paired if:
// - it is paired
// - it is mapped and its mate is mapped
// - the orientation and the distance between the mates are consistent with the sequencing protocol
func (flag *Flag) SetProperlyPaired() {
	flag.Value |= 1 << 1
}
func (flag *Flag) IsProperlyPaired() bool {
	return (flag.Value & (1 << 1)) != 0
}

// SetUnmapped sets the read as unmapped which means that the read could not be aligned to the reference
func (flag *Flag) SetUnmapped() {
	flag.Value |= 1 << 2
}
func (flag *Flag) IsUnmapped() bool {
	return (flag.Value & (1 << 2)) != 0
}

// SetMateUnmapped sets the mate as unmapped which means that the mate of the read could not be aligned to the reference
func (flag *Flag) SetMateUnmapped() {
	flag.Value |= 1 << 3
}
func (flag *Flag) IsMateUnmapped() bool {
	return (flag.Value & (1 << 3)) != 0
}

// SetReverseStrand sets the read as mapped to the reverse strand of the reference
func (flag *Flag) SetReverseStrand() {
	flag.Value |= 1 << 4
}
func (flag *Flag) IsReverseStrand() bool {
	return (flag.Value & (1 << 4)) != 0
}

// SetMateReverseStrand sets the mate as mapped to the reverse strand of the reference
func (flag *Flag) SetMateReverseStrand() {
	flag.Value |= 1 << 5
}
func (flag *Flag) IsMateReverseStrand() bool {
	return (flag.Value & (1 << 5)) != 0
}

// SetFirstInPair sets the read as the first in the pair which means that
func (flag *Flag) SetFirstInPair() {
	flag.Value |= 1 << 6
}
func (flag *Flag) IsFirstInPair() bool {
	return (flag.Value & (1 << 6)) != 0
}

// SetSecondInPair sets the read as the second in the pair which means that
func (flag *Flag) SetSecondInPair() {
	flag.Value |= 1 << 7
}
func (flag *Flag) IsSecondInPair() bool {
	return (flag.Value & (1 << 7)) != 0
}

// SetNotPrimaryAlignment sets the read as not primary alignment which means that this alignment is not
// the primary alignment for the read and another alignment is probably better.
// This flag is used when a read has multiple alignments and its mapping is ambiguous.
func (flag *Flag) SetNotPrimaryAlignment() {
	flag.Value |= 1 << 8
}
func (flag *Flag) IsNotPrimaryAlignment() bool {
	return (flag.Value & (1 << 8)) != 0
}

// SetFailsQualityCheck sets the read as failing quality check
func (flag *Flag) SetFailsQualityCheck() {
	flag.Value |= 1 << 9
}
func (flag *Flag) FailsQualityCheck() bool {
	return (flag.Value & (1 << 9)) != 0
}

// SetDuplicate sets the read as a duplicate which means that the read is a PCR or optical duplicate
func (flag *Flag) SetDuplicate() {
	flag.Value |= 1 << 10
}
func (flag *Flag) IsDuplicate() bool {
	return (flag.Value & (1 << 10)) != 0
}

// SetSupplementary sets the read as supplementary which means that the read is part of a chimeric alignment
func (flag *Flag) SetSupplementary() {
	flag.Value |= 1 << 11
}
func (flag *Flag) IsSupplementary() bool {
	return (flag.Value & (1 << 11)) != 0
}

type FlagBuilder struct {
	flag Flag
}

func NewFlagBuilder() FlagBuilder {
	return FlagBuilder{flag: Flag{}}
}
func (b *FlagBuilder) SetPaired() *FlagBuilder {
	b.flag.SetPaired()
	return b
}
func (b *FlagBuilder) SetProperlyPaired() *FlagBuilder {
	b.flag.SetProperlyPaired()
	return b
}
func (b *FlagBuilder) SetUnmapped() *FlagBuilder {
	b.flag.SetUnmapped()
	return b
}
func (b *FlagBuilder) SetMateUnmapped() *FlagBuilder {
	b.flag.SetMateUnmapped()
	return b
}
func (b *FlagBuilder) SetReverseStrand() *FlagBuilder {
	b.flag.SetReverseStrand()
	return b
}
func (b *FlagBuilder) SetMateReverseStrand() *FlagBuilder {
	b.flag.SetMateReverseStrand()
	return b
}
func (b *FlagBuilder) SetFirstInPair() *FlagBuilder {
	b.flag.SetFirstInPair()
	return b
}
func (b *FlagBuilder) SetSecondInPair() *FlagBuilder {
	b.flag.SetSecondInPair()
	return b
}
func (b *FlagBuilder) SetNotPrimaryAlignment() *FlagBuilder {
	b.flag.SetNotPrimaryAlignment()
	return b
}
func (b *FlagBuilder) SetFailsQualityCheck() *FlagBuilder {
	b.flag.SetFailsQualityCheck()
	return b
}
func (b *FlagBuilder) SetDuplicate() *FlagBuilder {
	b.flag.SetDuplicate()
	return b
}
func (b *FlagBuilder) SetSupplementary() *FlagBuilder {
	b.flag.SetSupplementary()
	return b
}

func (b *FlagBuilder) Build() Flag {
	return b.flag
}

// NOTE: I only added this for testing
type CigarOp struct {
	Length int
	Op     byte
}

func ParseCIGAR(cigar []byte) ([]CigarOp, error) {

	var ops []CigarOp

	start := 0
	for i := 0; i < len(cigar); i++ {
		// as soon as we hit a non-numeric byte, we can create a new cigarOp
		if cigar[i] < '0' || cigar[i] > '9' {
			cigraPrefix, err := strconv.Atoi(string(cigar[start:i]))
			if err != nil {
				return nil, err
			}
			ops = append(ops, CigarOp{Op: cigar[i], Length: cigraPrefix})
			start = i + 1
		}
	}

	return ops, nil
}

type SAMIndexEntry struct {
	Name           []byte
	Offset         int64
	FirstOfPair    bool
	DelimPositions [11]uint16 // this is hardcoded since we know how many fields each sam entry has. Using this array, we can quickly jump to a desired field
}

type SAMReadPair struct {
	First  *SAMIndexEntry
	Second *SAMIndexEntry
}

func ParseDelims(line []byte) ([]byte, [11]uint16, bool) {
	tabCount := 0
	delimPositions := [11]uint16{}
	var name []byte
	var flagStart int
	var isFist bool
	for i := 0; i < len(line); i++ {
		if line[i] == '\t' {

			if tabCount == 0 {
				name = make([]byte, i)
				copy(name, line[0:i])
				flagStart = i + 1
			}
			if tabCount == 1 {
				flagBytes := line[flagStart:i]
				flagInt, err := strconv.Atoi(string(flagBytes))
				if err != nil {
					panic("Error parsing flag of read")
				}
				isFist = flagInt&0x40 != 0
			}
			delimPositions[tabCount] = uint16(i)
			tabCount++
		}

	}
	return name, delimPositions, isFist
}

// ParseSparseSAMEntry takes a file and a sparseSAMEntry and parses the corresponding sam entry
func ParseSparseSAMEntry(file *os.File, sparseSAMEntry *SAMIndexEntry) (*mapperutils.ReadMatchResult, []byte, []byte, string, error) {
	if sparseSAMEntry == nil || (sparseSAMEntry.Offset < 0) {
		return nil, nil, nil, "", errors.New("invalid read: no offset available")
	}
	offset := sparseSAMEntry.Offset

	// jump to read pos in sam
	_, err := file.Seek(offset, io.SeekStart)
	if err != nil {
		return nil, nil, nil, "", err
	}

	reader := bufio.NewReader(file)
	line, err := reader.ReadBytes('\n')
	if err != nil && err != io.EOF {
		return nil, nil, nil, "", err
	}

	// remove trailing newline if present
	if len(line) > 0 && line[len(line)-1] == '\n' {
		line = line[:len(line)-1]
	}
	startBytes := line[sparseSAMEntry.DelimPositions[2]+1 : sparseSAMEntry.DelimPositions[3]]
	qnameBytes := line[0:sparseSAMEntry.DelimPositions[0]]
	qualityBytes := line[sparseSAMEntry.DelimPositions[3]+1 : sparseSAMEntry.DelimPositions[4]]
	cigarBytes := line[sparseSAMEntry.DelimPositions[4]+1 : sparseSAMEntry.DelimPositions[5]]
	sequenceBytes := line[sparseSAMEntry.DelimPositions[8]+1 : sparseSAMEntry.DelimPositions[9]]

	startPos, err := strconv.Atoi(string(startBytes))
	if err != nil {
		return nil, nil, nil, "", err
	}

	var sId int
	if sparseSAMEntry.FirstOfPair {
		sId = 0
	} else {
		sId = 1
	}

	result := &mapperutils.ReadMatchResult{
		SequenceIndex: sId,
		FourthPass:    false, // this is always false since the read is listed in the output sam
	}

	result.MatchedRead = regionvector.NewRegionVector()
	result.MatchedGenome = regionvector.NewRegionVector()

	// Process CIGAR string to build region vectors and identify mismatches
	err = processCigarForRegions(cigarBytes, startPos, sequenceBytes, result)
	if err != nil {
		return nil, nil, nil, "", err
	}

	return result, sequenceBytes, qualityBytes, string(qnameBytes), nil
}

// processCigarForRegions analyzes the CIGAR string to build read and genome region vectors
func processCigarForRegions(cigar []byte, refPos int, readSeq []byte, result *mapperutils.ReadMatchResult) error {
	readPos := 0

	// get cigar ops
	cigarOps, err := ParseCIGAR(cigar)
	if err != nil {
		return err
	}

	for _, op := range cigarOps {
		switch op.Op {
		case 'M', '=', 'X':
			result.MatchedRead.AddRegionNonOverlappingPanic(readPos, readPos+op.Length)
			result.MatchedGenome.AddRegionNonOverlappingPanic(refPos, refPos+op.Length)

			if op.Op == 'X' {
				for i := 0; i < op.Length; i++ {
					result.MismatchesRead = append(result.MismatchesRead, readPos+i)
				}
			}

			readPos += op.Length
			refPos += op.Length

		case 'I':
			readPos += op.Length

		case 'D', 'N':
			refPos += op.Length

		case 'S':
			readPos += op.Length

		case 'H':
			// nothing is consumed

		case 'P':
			// nothing is consumed
		}
	}

	return nil
}

func CreateSAMIndex(samfile *os.File) map[string]*SAMReadPair {

	scanner := bufio.NewScanner(samfile)
	const maxCapacity = 1024 * 1024 // 1MB
	buf := make([]byte, maxCapacity)
	scanner.Buffer(buf, maxCapacity)

	samIndex := make(map[string]*SAMReadPair)

	offset := int64(0)
	for scanner.Scan() {
		line := scanner.Bytes()
		if line[0] == '@' {
			offset += int64(len(line) + 1)
			continue
		}

		name, delimPositions, firstOfPair := parseDelims(line)
		indexEntry := &SAMIndexEntry{
			Name:           name,
			Offset:         offset,
			FirstOfPair:    firstOfPair,
			DelimPositions: delimPositions,
		}
		strName := string(name)
		_, exists := samIndex[strName]
		if !exists {
			samIndex[strName] = &SAMReadPair{}
		}

		if firstOfPair {
			samIndex[strName].First = indexEntry
		} else {
			samIndex[strName].Second = indexEntry
		}

		offset += int64(len(line)) + 1 // +1 for newline
	}

	return samIndex
}
func parseDelims(line []byte) ([]byte, [11]uint16, bool) {
	tabCount := 0
	delimPositions := [11]uint16{}
	var name []byte
	var flagStart int
	var isFist bool
	for i := 0; i < len(line); i++ {
		if line[i] == '\t' {

			if tabCount == 0 {
				name = make([]byte, i)
				copy(name, line[0:i])
				flagStart = i + 1
			}
			if tabCount == 1 {
				flagBytes := line[flagStart:i]
				flagInt, err := strconv.Atoi(string(flagBytes))
				if err != nil {
					panic("Error parsing flag of read")
				}
				isFist = flagInt&0x40 != 0
			}
			delimPositions[tabCount] = uint16(i)
			tabCount++
		}

	}
	return name, delimPositions, isFist
}
