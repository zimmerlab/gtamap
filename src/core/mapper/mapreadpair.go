package mapper

import (
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/core/timer"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
	"strconv"
	"strings"
)

func MapReadPair(readPair *fastq.ReadPair, genomeIndex *index.GenomeIndex,
	secondPassChan *mapperutils.SecondPassChannel,
	timerChannel chan<- *timer.Timer) (string, bool) {

	keepFw := Filter(readPair.ReadR1.Sequence, genomeIndex)
	keepRw := Filter(readPair.ReadR2.Sequence, genomeIndex)

	logrus.WithFields(logrus.Fields{
		"keepFw": keepFw,
		"keepRv": keepRw,
	}).Debug("Filter results")

	if !keepFw || !keepRw {
		return "", false
	}

	resultFw, isMappableFw := MapRead(readPair.ReadR1, genomeIndex)
	resultRv, isMappableRv := MapRead(readPair.ReadR2, genomeIndex)

	if !isMappableFw || !isMappableRv || len(resultFw) == 0 || len(resultRv) == 0 {
		return "", false
	}

	needSecondPass := false

	for i, resFw := range resultFw {
		if resFw.SecondPass {
			logrus.Debug("Second pass required for forward read: ", i)
			needSecondPass = true
		}
	}
	for i, resRv := range resultRv {
		if resRv.SecondPass {
			logrus.Debug("Second pass required for reverse read: ", i)
			needSecondPass = true
		}
	}

	if needSecondPass {
		secondPassChan.Send(&mapperutils.SecondPassTask{
			ReadPair: readPair,
			ResultFw: &resultFw,
			ResultRv: &resultRv,
		})
		return "", false
	}

	if len(resultFw) > 1 || len(resultRv) > 1 {
		// TODO: handle multimapping reads

		var builder strings.Builder

		builder.WriteString(readPair.ReadR1.Header)
		builder.WriteString("_fw\tMULTIMAPPING\t")
		for i, resFw := range resultFw {
			builder.WriteString(strconv.Itoa(i))
			builder.WriteString(":\t")
			builder.WriteString(resFw.GetCigar())
			builder.WriteString("\t")
			builder.WriteString(strconv.Itoa(resFw.MatchedRead.Length()))
			builder.WriteString("\t")
		}
		builder.WriteString("\n")

		builder.WriteString(readPair.ReadR2.Header)
		builder.WriteString("_rw\tMULTIMAPPING\t")
		for i, resRv := range resultRv {
			builder.WriteString(strconv.Itoa(i))
			builder.WriteString(":\t")
			builder.WriteString(resRv.GetCigar())
			builder.WriteString("\t")
			builder.WriteString(strconv.Itoa(resRv.MatchedRead.Length()))
			builder.WriteString("\t")
		}
		builder.WriteString("\n")

		return builder.String(), true
	}

	// multimapping was handled before so length is always 1
	resFw := resultFw[0]
	resRv := resultRv[0]

	postprocessReadMatch(genomeIndex, readPair.ReadR1, &resFw)
	postprocessReadMatch(genomeIndex, readPair.ReadR2, &resRv)

	flagFw := sam.Flag{}
	flagRv := sam.Flag{}

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

	var builder strings.Builder

	// QNAME
	// use only the string before any whitespace as the header
	headerFw := strings.Split(readPair.ReadR1.Header, " ")[0]
	builder.WriteString(headerFw)
	builder.WriteString("\t")
	// FLAG
	builder.WriteString(strconv.Itoa(flagFw.Value))
	builder.WriteString("\t")
	// RNAME
	builder.WriteString(genomeIndex.GetSequenceInfo(resFw.SequenceIndex).Contig)
	builder.WriteString("\t")
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

	startGenomeFw := geneStartRelativeToContig + offsetFw

	// +1 because the sam format is 1-based
	builder.WriteString(strconv.Itoa(startGenomeFw + 1))
	builder.WriteString("\t")
	// MAPQ
	builder.WriteString(strconv.Itoa(255))
	builder.WriteString("\t")
	// CIGAR
	builder.WriteString(resFw.GetCigar())
	builder.WriteString("\t")
	// RNEXT
	builder.WriteString("*")
	builder.WriteString("\t")
	// PNEXT
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// TLEN
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// SEQ
	if !flagFw.IsReverseStrand() {
		builder.WriteString(string(*readPair.ReadR1.Sequence))
	} else {
		builder.WriteString(string(utils.ReverseComplementDnaBytes(*readPair.ReadR1.Sequence)))
	}
	builder.WriteString("\t")
	// QUAL
	builder.WriteString(string(*readPair.ReadR1.Quality))
	builder.WriteString("\n")

	// R2 READ

	// QNAME
	headerRv := strings.Split(readPair.ReadR2.Header, " ")[0]
	builder.WriteString(headerRv)
	builder.WriteString("\t")
	// FLAG
	builder.WriteString(strconv.Itoa(flagRv.Value))
	builder.WriteString("\t")
	// RNAME
	builder.WriteString(genomeIndex.GetSequenceInfo(resRv.SequenceIndex).Contig)
	builder.WriteString("\t")
	// POS
	offsetRv := resRv.MatchedGenome.GetFirstRegion().Start

	if flagRv.IsReverseStrand() {
		geneLength := geneEndRelativeToContig - geneStartRelativeToContig
		offsetRv = geneLength - resRv.MatchedGenome.GetLastRegion().End
	}

	startGenomeRv := geneStartRelativeToContig + offsetRv

	// +1 because the sam format is 1-based
	builder.WriteString(strconv.Itoa(startGenomeRv + 1))
	builder.WriteString("\t")
	// MAPQ
	builder.WriteString(strconv.Itoa(255))
	builder.WriteString("\t")
	// CIGAR
	builder.WriteString(resRv.GetCigar())
	builder.WriteString("\t")
	// RNEXT
	builder.WriteString("*")
	builder.WriteString("\t")
	// PNEXT
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// TLEN
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// SEQ
	if !flagRv.IsReverseStrand() {
		builder.WriteString(string(*readPair.ReadR2.Sequence))
	} else {
		builder.WriteString(string(utils.ReverseComplementDnaBytes(*readPair.ReadR2.Sequence)))
	}
	builder.WriteString("\t")
	// QUAL
	builder.WriteString(string(*readPair.ReadR2.Quality))
	builder.WriteString("\n")

	return builder.String(), true
}

func postprocessReadMatch(genomeIndex *index.GenomeIndex, read *fastq.Read, result *mapperutils.ReadMatchResult) {
	applyLeftNormalization(genomeIndex, read, result)
}

// applyLeftNormalization shifts the mapping in genome to the left if there are gaps in the genome
// that have ambiguous bases (N) before and after the gap. If this is the case, the mapping can not
// be determined, and it is up to the mapper to decide where to place the read.
// This function shifts the gap to the leftmost position.
func applyLeftNormalization(genomeIndex *index.GenomeIndex, read *fastq.Read, result *mapperutils.ReadMatchResult) {

	if len(result.MatchedGenome.Regions) < 2 {
		logrus.WithFields(logrus.Fields{
			"read":      read.Header,
			"isForward": genomeIndex.IsSequenceForward(result.SequenceIndex),
		}).Debug("No left normalization required")
		return
	}

	logrus.WithFields(logrus.Fields{
		"read":      read.Header,
		"isForward": genomeIndex.IsSequenceForward(result.SequenceIndex),
	}).Debug("apply left normalization")

	logrus.WithFields(logrus.Fields{
		"genome": result.MatchedGenome,
	}).Debug("match before left normalization")

	if genomeIndex.IsSequenceForward(result.SequenceIndex) {
		determineLeftNormalizationShiftFw(genomeIndex, read, result)
	} else {
		determineLeftNormalizationShiftRv(genomeIndex, read, result)
	}

	logrus.WithFields(logrus.Fields{
		"genome": result.MatchedGenome,
	}).Debug("match after left normalization")
}

// determineLeftNormalizationShiftFw determines the left normalization shift for forward reads.
// It finds gaps in the genome and shifts the mapping in genome to the left based on the
// value of shift.
// The shift is the number of equal bases between the suffix of the read sequence before the gap
// and the suffix of the genome sequence within the gap.
// The maximum value of the shift is the size if the gap.
func determineLeftNormalizationShiftFw(genomeIndex *index.GenomeIndex, read *fastq.Read, result *mapperutils.ReadMatchResult) {

	logrus.Debug("determine left normalization for forward seq")

	// find gaps in genome (bases in reference that are not present in read)
	for i := 0; i < len(result.MatchedGenome.Regions)-1; i++ {

		gap := result.MatchedGenome.GetGap(i)

		// the rightmost position of the subsequence to be tested in the read
		// (last position in the read before the gap in forward direction)
		rRead, sizeErr := result.MatchedGenome.GetSizeLeftIncluding(i)
		if sizeErr != nil {
			logrus.Fatal("Error getting size left including", sizeErr)
		}

		// the rightmost position of the subsequence to be tested in the genome
		// (end of the gap in the genome)
		rGenome := gap.End

		// DEBUG
		//readSeqBeforeGap := (*read.Sequence)[rRead-gap.Length() : rRead]
		//fmt.Println("readSeqBeforeGap\t", string(readSeqBeforeGap))
		//genomeSeqInGap := (*genomeIndex.Sequences[result.SequenceIndex])[rGenome-gap.Length() : rGenome]
		//fmt.Println("genomeSeqInGap\t\t", string(genomeSeqInGap))

		shift := 0

		for i := 0; i < gap.Length(); i++ {

			charRead := (*read.Sequence)[rRead-1-i]

			charGenome := (*genomeIndex.Sequences[result.SequenceIndex])[rGenome-1-i]

			if charRead != charGenome {
				break
			}

			shift++
		}

		logrus.WithFields(logrus.Fields{
			"shift": shift,
		}).Debug("determined shift")

		// shift the mapping in genome to left based on value of shift
		result.MatchedGenome.Regions[i].End -= shift
		result.MatchedGenome.Regions[i+1].Start -= shift
	}
}

// determineLeftNormalizationShiftRv determines the left normalization shift for reverse reads.
// The procedure is equivalent to the one for forward reads, but the direction of the
// comparison is reversed.
// The shift is now determined by the number of equal bases between the prefix of the read sequence
// after the gap and the prefix of the genome sequence within the gap.
func determineLeftNormalizationShiftRv(genomeIndex *index.GenomeIndex, read *fastq.Read, result *mapperutils.ReadMatchResult) {

	logrus.Debug("determine left normalization for reverse seq")

	// find gaps in genome (bases in reference that are not present in read)
	for i := 0; i < len(result.MatchedGenome.Regions)-1; i++ {

		gap := result.MatchedGenome.GetGap(i)

		// the rightmost position of the subsequence to be tested in the read
		// (last position in the read before the gap in forward direction)
		// actually this is the first position in the read after the gap (because rv)
		rRead, sizeErr := result.MatchedGenome.GetSizeLeftIncluding(i)
		if sizeErr != nil {
			logrus.Fatal("Error getting size left including", sizeErr)
		}

		// the rightmost position of the subsequence to be tested in the genome
		// (end of the gap in the genome in forward direction)
		// actually this is the first position in the gap (because rv)
		rGenome := gap.Start

		// DEBUG
		//readSeqAfterGap := (*read.Sequence)[rRead : rRead+gap.Length()]
		//fmt.Println("readSeqAfterGap\t", string(readSeqAfterGap))
		//genomeSeqInGap := (*genomeIndex.Sequences[result.SequenceIndex])[rGenome : rGenome+gap.Length()]
		//fmt.Println("genomeSeqInGap\t", string(genomeSeqInGap))

		shift := 0

		for i := 0; i < gap.Length(); i++ {

			charRead := (*read.Sequence)[rRead+i]

			charGenome := (*genomeIndex.Sequences[result.SequenceIndex])[rGenome+i]

			if charRead != charGenome {
				break
			}

			shift++
		}

		logrus.WithFields(logrus.Fields{
			"shift": shift,
		}).Debug("determined shift")

		// shift the mapping in genome to right based on value of shift
		result.MatchedGenome.Regions[i].End += shift
		result.MatchedGenome.Regions[i+1].Start += shift
	}
}
