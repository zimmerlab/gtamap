package sam

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
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
	Qname         string   // id of the read
	Flag          Flag     // bitwise flag
	Rname         string   // id of the reference
	Pos           int      // 1-based leftmost mapping position
	Mapq          int      // mapping quality
	Cigar         string   // CIGAR string (compact representation of alignment)
	Rnext         string   // id of the mate reference
	Pnext         int      // position of the mate reference
	Tlen          int      // observed template length
	Seq           string   // read sequence
	Qual          string   // read quality
	Tags          []string // optional tags (e.g. NM, XT, etc.)
	TranscriptIds []int    // ids of the transcripts that this read was aligned to
}

func (entry *Record) UniformCigar() string {
	// replace all occurrences of "X" and "=" in the CIGAR string with "M"
	cigar := strings.ReplaceAll(entry.Cigar, "X", "M")
	cigar = strings.ReplaceAll(cigar, "=", "M")
	return cigar
}

func (entry *Record) String() string {
	if entry == nil {
		return ""
	}

	alignmentLineString := fmt.Sprintf("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s", entry.Qname, entry.Flag.String(), entry.Rname, entry.Pos, entry.Mapq, entry.Cigar, entry.Rnext, entry.Pnext, entry.Tlen, entry.Seq, entry.Qual)
	// alignmentLineString += fmt.Sprintf("\tXT:Z:T%d", entry.TranscriptId)

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

type ParsedHeader struct {
	Type string
	Tags map[string]string
}

type ParsedFile struct {
	Headers []ParsedHeader
	Records []Record
}

func ParseSAMFile(filename string) (*ParsedFile, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("error opening file: %v", err)
	}
	defer file.Close()

	samFile := &ParsedFile{
		Headers: make([]ParsedHeader, 0),
		Records: make([]Record, 0),
	}

	scanner := bufio.NewScanner(file)
	lineNum := 0

	for scanner.Scan() {
		lineNum++
		line := strings.TrimSpace(scanner.Text())

		// Skip empty lines
		if len(line) == 0 {
			continue
		}

		// Parse header lines (start with @)
		if strings.HasPrefix(line, "@") {
			header, err := parseHeader(line)
			if err != nil {
				return nil, fmt.Errorf("error parsing header at line %d: %v", lineNum, err)
			}
			samFile.Headers = append(samFile.Headers, header)
		} else {
			// Parse alignment record
			record, err := parseRecord(line)
			if err != nil {
				return nil, fmt.Errorf("error parsing record at line %d: %v", lineNum, err)
			}
			samFile.Records = append(samFile.Records, record)
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading file: %v", err)
	}

	return samFile, nil
}

// parseHeader parses a SAM header line
func parseHeader(line string) (ParsedHeader, error) {
	fields := strings.Split(line, "\t")
	if len(fields) < 1 {
		return ParsedHeader{}, fmt.Errorf("invalid header format")
	}

	headerType := fields[0][1:] // Remove the @ prefix
	tags := make(map[string]string)

	for i := 1; i < len(fields); i++ {
		parts := strings.SplitN(fields[i], ":", 2)
		if len(parts) == 2 {
			tags[parts[0]] = parts[1]
		}
	}

	return ParsedHeader{
		Type: headerType,
		Tags: tags,
	}, nil
}

// parseRecord parses a SAM alignment record
func parseRecord(line string) (Record, error) {
	fields := strings.Split(line, "\t")
	if len(fields) < 11 {
		return Record{}, fmt.Errorf("SAM record must have at least 11 fields, got %d", len(fields))
	}

	// Parse FLAG
	flag, err := strconv.Atoi(fields[1])
	if err != nil {
		return Record{}, fmt.Errorf("invalid FLAG field: %v", err)
	}

	// Parse POS
	pos, err := strconv.Atoi(fields[3])
	if err != nil {
		return Record{}, fmt.Errorf("invalid POS field: %v", err)
	}

	// Parse MAPQ
	mapq, err := strconv.Atoi(fields[4])
	if err != nil {
		return Record{}, fmt.Errorf("invalid MAPQ field: %v", err)
	}

	// Parse PNEXT
	pnext, err := strconv.Atoi(fields[7])
	if err != nil {
		return Record{}, fmt.Errorf("invalid PNEXT field: %v", err)
	}

	// Parse TLEN
	tlen, err := strconv.Atoi(fields[8])
	if err != nil {
		return Record{}, fmt.Errorf("invalid TLEN field: %v", err)
	}

	// Collect optional tags
	var tags []string
	if len(fields) > 11 {
		tags = fields[11:]
	}

	return Record{
		Qname: fields[0],
		Flag:  Flag{flag},
		Rname: fields[2],
		Pos:   pos,
		Mapq:  mapq,
		Cigar: fields[5],
		Rnext: fields[6],
		Pnext: pnext,
		Tlen:  tlen,
		Seq:   fields[9],
		Qual:  fields[10],
		Tags:  tags,
	}, nil
}
