package sam

import (
	"fmt"
	"math"
	"strconv"
	"strings"

	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
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

func FormatMappedReadPairToSAM(resFw mapperutils.ReadMatchResult, resRv mapperutils.ReadMatchResult, readPair *fastq.ReadPair, genomeIndex *index.GenomeIndex) (string, bool) {
	flagFw := Flag{}
	flagRv := Flag{}

	flagFw.SetPaired()
	flagFw.SetProperlyPaired()
	flagFw.SetFirstInPair()
	if !genomeIndex.IsSequenceForward(resFw.SequenceIndex) {
		flagFw.SetReverseStrand()
		flagRv.SetMateReverseStrand()
	}

	flagRv.SetPaired()
	flagRv.SetProperlyPaired()
	flagRv.SetSecondInPair()
	if !genomeIndex.IsSequenceForward(resRv.SequenceIndex) {
		flagRv.SetReverseStrand()
		flagFw.SetMateReverseStrand()
	}

	geneStartRelativeToContig := int(genomeIndex.GetSequenceInfo(resFw.SequenceIndex).StartGenomic)
	geneEndRelativeToContig := int(genomeIndex.GetSequenceInfo(resRv.SequenceIndex).EndGenomic)

	// QNAME
	// use only the string before any whitespace as the header
	headerFw := strings.Split(readPair.ReadR1.Header, " ")[0]
	headerRv := strings.Split(readPair.ReadR2.Header, " ")[0]

	// FLAG
	flagFwStr := strconv.Itoa(flagFw.Value)
	flagRvStr := strconv.Itoa(flagRv.Value)

	// RNAME
	rnameFw := genomeIndex.GetSequenceInfo(resFw.SequenceIndex).Contig
	rnameRv := genomeIndex.GetSequenceInfo(resRv.SequenceIndex).Contig

	// POS
	// the offset is the start of the first region in the matched genome which corresponds to the first
	// mapped position of the read
	offsetFw := resFw.MatchedGenome.GetFirstRegion().Start
	// if the read is mapped to the reverse strand, we need to calculate the offset because the target sequence
	// was reverse complemented. therefore the start position in 5'-3' direction of the original sequence is the
	// end position of the reverse complemented sequence
	if flagFw.IsReverseStrand() {
		geneLength := geneEndRelativeToContig - geneStartRelativeToContig
		offsetFw = geneLength - resFw.MatchedGenome.GetLastRegion().End
	}
	// +1 because the sam format is 1-based
	startGenomeFw := geneStartRelativeToContig + offsetFw + 1

	offsetRv := resRv.MatchedGenome.GetFirstRegion().Start
	if flagRv.IsReverseStrand() {
		geneLength := geneEndRelativeToContig - geneStartRelativeToContig
		offsetRv = geneLength - resRv.MatchedGenome.GetLastRegion().End
	}
	// +1 because the sam format is 1-based
	startGenomeRv := geneStartRelativeToContig + offsetRv + 1

	// MAPQ
	// https://samtools.github.io/hts-specs/SAMv1.pdf
	// No alignments should be assigned mapping quality 255
	mapqFw := 254
	mapqRv := 254

	// CIGAR
	// https://samtools.github.io/hts-specs/SAMv1.pdf
	// Adjacent CIGAR operations should be different
	cigarFw := resFw.GetCigar()
	cigarRv := resRv.GetCigar()

	// RNEXT
	rnextFw := rnameRv
	rnextRv := rnameFw

	// PNEXT
	pnextFw := startGenomeRv
	pnextRv := startGenomeFw

	// TLEN
	// signed observed Template LENgth. For primary reads where the primary alignments of all reads
	// in the template are mapped to the same reference sequence, the absolute value of TLEN equals the
	// distance between the mapped end of the template and the mapped start of the template, inclusively
	// (i.e., end − start + 1).15 Note that mapped base is defined to be one that aligns to the reference as
	// described by CIGAR, hence excludes soft-clipped bases. The TLEN field is positive for the leftmost
	// segment of the template, negative for the rightmost, and the sign for any middle segment is undefined.
	// If segments cover the same coordinates then the choice of which is leftmost and rightmost is arbitrary
	// but the two ends must still have differing signs. It is set as 0 for a single-segment template or when
	// the information is unavailable (e.g., when the first or last segment of a multi-segment template is
	// unmapped or when the two are mapped to different reference sequences).
	// The intention of this field is to indicate where the other end of the template has been aligned without
	// needing to read the remainder of the SAM file. Unfortunately there has been no clear consensus on
	// the definitions of the template mapped start and end. Thus the exact definitions are implementation-defined.
	tlen := int(math.Abs(float64(startGenomeFw-startGenomeRv))) + len(*readPair.ReadR1.Sequence)
	tlenFw := tlen
	if startGenomeFw > startGenomeRv {
		tlenFw = -1 * tlen
	}
	tlenRv := tlen
	if startGenomeRv >= startGenomeFw {
		tlenRv = -1 * tlen
	}

	// SEQ
	seqFw := string(*readPair.ReadR1.Sequence)
	if flagFw.IsReverseStrand() {
		revCompSeq, revCompSeqErr := utils.ReverseComplementDnaBytes(*readPair.ReadR1.Sequence)
		if revCompSeqErr != nil {
			logrus.WithFields(logrus.Fields{
				"read": readPair.ReadR1.Header,
			}).Error("Error reversing complementing sequence", revCompSeqErr)
			return "", false
		}
		seqFw = string(revCompSeq)
	}

	seqRv := string(*readPair.ReadR2.Sequence)
	if flagRv.IsReverseStrand() {
		revCompSeq, revCompSeqErr := utils.ReverseComplementDnaBytes(*readPair.ReadR2.Sequence)
		if revCompSeqErr != nil {
			logrus.WithFields(logrus.Fields{
				"read": readPair.ReadR2.Header,
			}).Error("Error reversing complementing sequence", revCompSeqErr)
			return "", false
		}
		seqRv = string(revCompSeq)
	}

	// QUAL
	qualFw := string(*readPair.ReadR1.Quality)
	if flagFw.IsReverseStrand() {
		qualFw = string(utils.ReverseBytes(*readPair.ReadR1.Quality))
	}

	qualRv := string(*readPair.ReadR2.Quality)
	if flagRv.IsReverseStrand() {
		qualRv = string(utils.ReverseBytes(*readPair.ReadR2.Quality))
	}

	// ATTRIBUTES

	var builder strings.Builder

	// QNAME
	builder.WriteString(headerFw)
	builder.WriteString("\t")
	// FLAG
	builder.WriteString(flagFwStr)
	builder.WriteString("\t")
	// RNAME
	builder.WriteString(rnameFw)
	builder.WriteString("\t")
	// POS
	builder.WriteString(strconv.Itoa(startGenomeFw))
	builder.WriteString("\t")
	// MAPQ
	builder.WriteString(strconv.Itoa(mapqFw))
	builder.WriteString("\t")
	// CIGAR
	builder.WriteString(cigarFw)
	builder.WriteString("\t")
	// RNEXT
	builder.WriteString(rnextFw)
	builder.WriteString("\t")
	// PNEXT
	builder.WriteString(strconv.Itoa(pnextFw))
	builder.WriteString("\t")
	// TLEN
	builder.WriteString(strconv.Itoa(tlenFw))
	builder.WriteString("\t")
	// SEQ
	builder.WriteString(seqFw)
	builder.WriteString("\t")
	// QUAL
	builder.WriteString(qualFw)
	builder.WriteString("\n")

	// R2 READ

	// QNAME
	builder.WriteString(headerRv)
	builder.WriteString("\t")
	// FLAG
	builder.WriteString(flagRvStr)
	builder.WriteString("\t")
	// RNAME
	builder.WriteString(rnameRv)
	builder.WriteString("\t")
	// POS
	builder.WriteString(strconv.Itoa(startGenomeRv))
	builder.WriteString("\t")
	// MAPQ
	builder.WriteString(strconv.Itoa(mapqRv))
	builder.WriteString("\t")
	// CIGAR
	builder.WriteString(cigarRv)
	builder.WriteString("\t")
	// RNEXT
	builder.WriteString(rnextRv)
	builder.WriteString("\t")
	// PNEXT
	builder.WriteString(strconv.Itoa(pnextRv))
	builder.WriteString("\t")
	// TLEN
	builder.WriteString(strconv.Itoa(tlenRv))
	builder.WriteString("\t")
	// SEQ
	builder.WriteString(seqRv)
	builder.WriteString("\t")
	// QUAL
	builder.WriteString(qualRv)
	builder.WriteString("\n")

	return builder.String(), true
}
