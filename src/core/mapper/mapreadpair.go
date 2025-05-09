package mapper

import (
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mappedreadpair"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/core/timer"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
	"math"
	"strconv"
	"strings"
)

func MapReadPair(readPair *fastq.ReadPair, genomeIndex *index.GenomeIndex,
	fourthPassChan *mapperutils.FourthPassChannel,
	timerChannel chan<- *timer.Timer) ([]mappedreadpair.ReadPairMatchResult, bool) {

	keepFw := Filter(readPair.ReadR1.Sequence, genomeIndex)
	keepRw := Filter(readPair.ReadR2.Sequence, genomeIndex)

	logrus.WithFields(logrus.Fields{
		"keepFw": keepFw,
		"keepRv": keepRw,
	}).Debug("Filter results")

	if !keepFw || !keepRw {
		return nil, false
	}

	resultFw, isMappableFw := MapRead(readPair.ReadR1, genomeIndex)
	resultRv, isMappableRv := MapRead(readPair.ReadR2, genomeIndex)

	// TODO: REMOVE DEBUG
	//if isMappableFw {
	//	debugout.GenerateAlignmentView(genomeIndex, resultFw[0], readPair.ReadR1)
	//}
	//if isMappableRv {
	//	debugout.GenerateAlignmentView(genomeIndex, resultRv[0], readPair.ReadR2)
	//}

	if !isMappableFw || !isMappableRv || len(resultFw) == 0 || len(resultRv) == 0 {
		logrus.WithFields(logrus.Fields{
			"isMappableFw":  isMappableFw,
			"isMappableRv":  isMappableRv,
			"num resultsFw": len(resultFw),
			"num resultsRv": len(resultRv),
		}).Debug("readpair not mappable")
		return nil, false
	}

	needFourthPass := false

	for i, resFw := range resultFw {
		if resFw.FourthPass {
			logrus.Debug("Fourth pass required for forward read: ", i)
			needFourthPass = true
		}
	}
	for i, resRv := range resultRv {
		if resRv.FourthPass {
			logrus.Debug("Fourth pass required for reverse read: ", i)
			needFourthPass = true
		}
	}

	if needFourthPass {
		fourthPassChan.Send(&mapperutils.FourthPassTask{
			ReadPair: readPair,
			ResultFw: &resultFw,
			ResultRv: &resultRv,
		})
		return nil, false
	}

	if len(resultFw) > 1 || len(resultRv) > 1 {
		// TODO: handle multimapping reads

		logrus.WithFields(logrus.Fields{
			"numResultsFw": len(resultFw),
			"numResultsRv": len(resultRv),
			"readFw":       readPair.ReadR1.Header,
			"readRv":       readPair.ReadR2.Header,
		}).Warn("multimapping reads not handled yet")

		//var builder strings.Builder
		//
		//builder.WriteString(readPair.ReadR1.Header)
		//builder.WriteString("_fw\tMULTIMAPPING\t")
		//for i, resFw := range resultFw {
		//	builder.WriteString(strconv.Itoa(i))
		//	builder.WriteString(":\t")
		//	builder.WriteString(resFw.GetCigar())
		//	builder.WriteString("\t")
		//	builder.WriteString(strconv.Itoa(resFw.MatchedRead.Length()))
		//	builder.WriteString("\t")
		//}
		//builder.WriteString("\n")
		//
		//builder.WriteString(readPair.ReadR2.Header)
		//builder.WriteString("_rw\tMULTIMAPPING\t")
		//for i, resRv := range resultRv {
		//	builder.WriteString(strconv.Itoa(i))
		//	builder.WriteString(":\t")
		//	builder.WriteString(resRv.GetCigar())
		//	builder.WriteString("\t")
		//	builder.WriteString(strconv.Itoa(resRv.MatchedRead.Length()))
		//	builder.WriteString("\t")
		//}
		//builder.WriteString("\n")
		//
		//return builder.String(), true

		return nil, false
	}

	// multimapping was handled before so length is always 1
	//postprocessReadMatch(genomeIndex, readPair.ReadR1, &resFw)
	//postprocessReadMatch(genomeIndex, readPair.ReadR2, &resRv)

	//return FormatMappedReadPairToSAM(resFw, resRv, readPair, genomeIndex)
	mappedReadPairs := make([]mappedreadpair.ReadPairMatchResult, 0)
	for i, _ := range resultFw {
		mappedReadPairs = append(mappedReadPairs, mappedreadpair.ReadPairMatchResult{
			ReadPair: readPair,
			Fw:       &resultFw[i],
			Rv:       &resultRv[i],
			Index:    genomeIndex,
		})
	}
	return mappedReadPairs, true
}

func FormatMappedReadPairToSAM(resFw mapperutils.ReadMatchResult, resRv mapperutils.ReadMatchResult, readPair *fastq.ReadPair, genomeIndex *index.GenomeIndex) (string, bool) {

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
func determineLeftNormalizationShiftFw(
	genomeIndex *index.GenomeIndex,
	read *fastq.Read,
	result *mapperutils.ReadMatchResult) {

	logrus.Debug("determine left normalization for forward seq")

	regionIndexBeforeGap := result.MatchedGenome.GetGapIndexAfterPos(0)

	// find gaps in genome (bases in reference that are not present in read)
	for regionIndexBeforeGap > -1 {

		gapGenome := result.MatchedGenome.GetGapAfterRegionIndex(regionIndexBeforeGap)

		// only handle gaps in genome that have no mutual gap in read
		// (only introns or deletions)
		gapRead := result.MatchedRead.GetGapAfterRegionIndex(regionIndexBeforeGap)
		if gapRead != nil {
			logrus.WithFields(logrus.Fields{
				"gapRead":   gapRead,
				"gapGenome": gapGenome,
			}).Debug("skip gap in genome that has gap in read")
			regionIndexBeforeGap = result.MatchedGenome.GetGapIndexAfterPos(gapGenome.End + 1)
			continue
		}

		// the rightmost position of the subsequence to be tested in the read
		// (last position in the read before the gap in forward direction)
		rRead, sizeErr := result.MatchedGenome.GetSizeLeftIncluding(regionIndexBeforeGap)
		if sizeErr != nil {
			logrus.Fatal("Error getting size left including", sizeErr)
		}

		// the rightmost position of the subsequence to be tested in the genome
		// (end of the gap in the genome)
		rGenome := gapGenome.End

		// the amount of bases that the gap must be shifted to be left normalized
		shift := 0
		// the maximum amount that the gap can be shifted is until the beginning of the read
		maxShift := rRead

		for i := 0; i < maxShift; i++ {

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
		result.MatchedGenome.Regions[regionIndexBeforeGap].End -= shift
		result.MatchedGenome.Regions[regionIndexBeforeGap+1].Start -= shift

		regionIndexBeforeGap = result.MatchedGenome.GetGapIndexAfterPos(gapGenome.End + 1)
	}
}

// determineLeftNormalizationShiftRv determines the left normalization shift for reverse reads.
// The procedure is equivalent to the one for forward reads, but the direction of the
// comparison is reversed.
// The shift is now determined by the number of equal bases between the prefix of the read sequence
// after the gap and the prefix of the genome sequence within the gap.
func determineLeftNormalizationShiftRv(
	genomeIndex *index.GenomeIndex,
	read *fastq.Read,
	result *mapperutils.ReadMatchResult) {

	logrus.Debug("determine left normalization for reverse seq")

	regionIndexBeforeGap := result.MatchedGenome.GetGapIndexAfterPos(0)

	// find gaps in genome (bases in reference that are not present in read)
	for regionIndexBeforeGap > -1 {

		gapGenome := result.MatchedGenome.GetGapAfterRegionIndex(regionIndexBeforeGap)

		// only handle gaps in genome that have no mutual gap in read
		// (only introns or deletions)
		gapRead := result.MatchedRead.GetGapAfterRegionIndex(regionIndexBeforeGap)
		if gapRead != nil {
			logrus.WithFields(logrus.Fields{
				"gapRead":   gapRead,
				"gapGenome": gapGenome,
			}).Debug("skip gap in genome that has gap in read")
			regionIndexBeforeGap = result.MatchedGenome.GetGapIndexAfterPos(gapGenome.End + 1)
			continue
		}

		// the rightmost position of the subsequence to be tested in the read
		// (last position in the read before the gap in forward direction)
		// actually this is the first position in the read after the gap (because rv)
		rRead, sizeErr := result.MatchedGenome.GetSizeLeftIncluding(regionIndexBeforeGap)
		if sizeErr != nil {
			logrus.Fatal("Error getting size left including", sizeErr)
		}

		// the rightmost position of the subsequence to be tested in the genome
		// (end of the gap in the genome in forward direction)
		// actually this is the first position in the gap (because rv)
		rGenome := gapGenome.Start

		// the amount of bases that the gap must be shifted to be left normalized
		shift := 0
		// the maximum amount that the gap can be shifted is until the end of the read
		maxShift := result.MatchedGenome.Length() - rRead

		for i := 0; i < maxShift; i++ {

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
		//result.MatchedGenome.Regions[i].End += shift
		//result.MatchedGenome.Regions[i+1].Start += shift

		result.MatchedGenome.Regions[regionIndexBeforeGap].End += shift
		result.MatchedGenome.Regions[regionIndexBeforeGap+1].Start += shift

		regionIndexBeforeGap = result.MatchedGenome.GetGapIndexAfterPos(gapGenome.End + 1)
	}
}
