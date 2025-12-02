package thirdpass

import (
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
	"sync"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
)

func ThirdPassWorker(
	thirdPassChan *ThirdPassChannel,
	wgThirdPass *sync.WaitGroup,
	outputChan chan<- string,
	index *index.GenomeIndex,
) {
	defer wgThirdPass.Done()
	total := 0
	mmTotal := 0

	f, err := os.Create(config.Mapper.GetMappingStatsFilePath())
	if err != nil {
		panic(err)
	}
	defer f.Close()

	if err := WriteHeader(f); err != nil {
		panic(err)
	}

	for {

		task, ok := thirdPassChan.Receive()
		if !ok {
			logrus.Info("Done with second pass")
			break
		}
		// logrus.Debugf("Third pass: %s", task.ReadPairId)

		var builder strings.Builder

		if config.Mapper.Mapping.Output.IncludeAllPairings {
			for i := 0; i < len(task.TargetInfo.Fw); i++ {

				if task.TargetInfo.Fw[i].IncompleteMap {
					continue
				}

				for j := 0; j < len(task.TargetInfo.Rv); j++ {

					if task.TargetInfo.Rv[j].IncompleteMap {
						continue
					}

					total += len(*task.TargetInfo.ReadPair.ReadR1.Sequence)
					mmTotal += len(task.TargetInfo.Fw[i].MismatchesRead)
					total += len(*task.TargetInfo.ReadPair.ReadR1.Sequence)
					mmTotal += len(task.TargetInfo.Rv[j].MismatchesRead)

					s, err := readPairResultToSamString(
						index,
						task.TargetInfo.ReadPair,
						task.TargetInfo.Fw[i],
						task.TargetInfo.Rv[j],
					)
					if err != nil {
						logrus.Error("Error converting read pair result to SAM string: ", err)
						continue
					}
					builder.WriteString(s)
				}
			}
		} else {

			// here we do not pair any reads, we just output single sam entries
			// FW

			numRecordsR1 := 0

			for i := 0; i < len(task.TargetInfo.Fw); i++ {

				if task.TargetInfo.Fw[i].IncompleteMap {
					continue
				}

				total += len(*task.TargetInfo.ReadPair.ReadR1.Sequence)
				mmTotal += len(task.TargetInfo.Fw[i].MismatchesRead)

				s, err := readPairResultToSamString(
					index,
					task.TargetInfo.ReadPair,
					task.TargetInfo.Fw[i],
					nil,
				)

				task.TargetInfo.Fw[i].WriteTSV(
					f,
					index.GetSequenceInfo(
						task.TargetInfo.Fw[i].SequenceIndex,
					).GeneId,
					1,
					i,
					task.TargetInfo.ReadPair.ReadR1.Header,
				)

				if err != nil {
					logrus.Error("Error converting read pair result to SAM "+
						"string: ", err)
					continue
				}

				builder.WriteString(s)

				numRecordsR1++
			}

			// RV
			numRecordsR2 := 0

			if task.TargetInfo.ReadPair.ReadR2 != nil {
				for j := 0; j < len(task.TargetInfo.Rv); j++ {

					if task.TargetInfo.Rv[j].IncompleteMap {
						continue
					}

					total += len(*task.TargetInfo.ReadPair.ReadR2.Sequence)
					mmTotal += len(task.TargetInfo.Rv[j].MismatchesRead)

					s, err := readPairResultToSamString(
						index,
						task.TargetInfo.ReadPair,
						nil,
						task.TargetInfo.Rv[j],
					)

					task.TargetInfo.Rv[j].WriteTSV(
						f,
						index.GetSequenceInfo(
							task.TargetInfo.Rv[j].SequenceIndex,
						).GeneId,
						1,
						j,
						task.TargetInfo.ReadPair.ReadR2.Header,
					)

					if err != nil {
						logrus.Error("Error converting read pair result to SAM string: ", err)
						continue
					}

					builder.WriteString(s)

					numRecordsR2++
				}

				if (numRecordsR1 == 0 && numRecordsR2 != 0) ||
					(numRecordsR1 != 0 && numRecordsR2 == 0) {
					// one mate is unmapped
					// write read mate as unmapped to sam
					builder.WriteString(
						unmappedReadMateToSamString(
							task.TargetInfo.ReadPair,
							numRecordsR1 == 0,
						),
					)
				}
			}
		}

		outputChan <- builder.String()
	}

	logrus.Info("Done with output")
	logrus.WithFields(logrus.Fields{
		"Average MM":        float64(mmTotal) / float64(total),
		"Aligned Positions": total,
		"Mismatches":        mmTotal,
	}).Info("Stats")
}

func unmappedReadMateToSamString(
	readPair *fastq.ReadPair,
	isR1 bool,
) string {
	flag := sam.Flag{}

	flag.SetUnmapped()
	flag.SetPaired()

	var header string
	var seq string
	var qual string

	if isR1 {
		flag.SetFirstInPair()
		header = strings.Split(readPair.ReadR1.Header, " ")[0]
		seq = string(*readPair.ReadR1.Sequence)
		qual = string(*readPair.ReadR1.Quality)
	} else {
		flag.SetSecondInPair()
		header = strings.Split(readPair.ReadR2.Header, " ")[0]
		seq = string(*readPair.ReadR2.Sequence)
		qual = string(*readPair.ReadR2.Quality)
	}

	var builder strings.Builder

	// QNAME
	builder.WriteString(header)
	builder.WriteString("\t")
	// FLAG
	builder.WriteString(flag.String())
	builder.WriteString("\t")
	// RNAME
	builder.WriteString("*")
	builder.WriteString("\t")
	// POS
	builder.WriteString("0")
	builder.WriteString("\t")
	// MAPQ
	builder.WriteString("0")
	builder.WriteString("\t")
	// CIGAR
	builder.WriteString("*")
	builder.WriteString("\t")
	// RNEXT
	builder.WriteString("*")
	builder.WriteString("\t")
	// PNEXT
	builder.WriteString("0")
	builder.WriteString("\t")
	// TLEN
	builder.WriteString("0")
	builder.WriteString("\t")
	// SEQ
	builder.WriteString(seq)
	builder.WriteString("\t")
	// QUAL
	builder.WriteString(qual)
	builder.WriteString("\n")

	return builder.String()
}

func readSingleResultToSamString(
	genomeIndex *index.GenomeIndex,
	read *fastq.Read,
	resFw *mapperutils.ReadMatchResult,
) (
	string,
	error,
) {
	flag := sam.Flag{}

	if resFw == nil {
		// TODO: handle parameter for writing unmapped results to output
		return "", fmt.Errorf("No mapping result provided for single read: %s", read.Header)
	}

	if !genomeIndex.IsSequenceForward(resFw.SequenceIndex) {
		flag.SetReverseStrand()
	}

	geneStartRelativeToContig := int(genomeIndex.GetSequenceInfo(resFw.SequenceIndex).StartGenomic)
	geneEndRelativeToContig := int(genomeIndex.GetSequenceInfo(resFw.SequenceIndex).EndGenomic)

	headerFw := strings.Split(read.Header, " ")[0]

	// FLAG
	flagStr := strconv.Itoa(flag.Value)

	// RNAME
	// the name of the reference sequence (contig) to which the read is aligned
	// it is supposed that both pairs must map to the same contig
	rname := genomeIndex.GetSequenceInfo(resFw.SequenceIndex).Contig

	// POS
	// the offset is the start of the first region in the matched genome which corresponds to the first
	// mapped position of the read
	startGenomeFw := 0
	firstRegionGenome, _ := resFw.MatchedGenome.GetFirstRegion()
	offsetFw := firstRegionGenome.Start
	// if the read is mapped to the reverse strand, we need to calculate the offset because the target sequence
	// was reverse complemented. therefore the start position in 5'-3' direction of the original sequence is the
	// end position of the reverse complemented sequence
	if flag.IsReverseStrand() {
		geneLength := geneEndRelativeToContig - geneStartRelativeToContig
		lastRegionGenome, _ := resFw.MatchedGenome.GetLastRegion()
		offsetFw = geneLength - lastRegionGenome.End
	}
	// +1 because the sam format is 1-based
	startGenomeFw = geneStartRelativeToContig + offsetFw + 1

	// MAPQ
	// https://samtools.github.io/hts-specs/SAMv1.pdf
	// No alignments should be assigned mapping quality 255
	// TODO: implement mapping quality
	mapqFw := 254

	// CIGAR
	// https://samtools.github.io/hts-specs/SAMv1.pdf
	// Adjacent CIGAR operations should be different

	cigar, errCigar := resFw.GetCigar()
	if errCigar != nil {
		logrus.WithFields(logrus.Fields{
			"read":              read.Header,
			"length of mapping": resFw.MatchedGenome.Length(),
			"expected length":   len(*read.Sequence),
		}).Error("Error getting CIGAR string of FW read: ", errCigar)

		// TODO: handle error
		// return "", false
		cigar = "*"
	}

	// SEQ
	seqFw := string(*read.Sequence)
	if flag.IsReverseStrand() {
		revCompSeq, revCompSeqErr := utils.ReverseComplementDnaBytes(*read.Sequence)
		if revCompSeqErr != nil {
			logrus.WithFields(logrus.Fields{
				"read": read.Header,
			}).Error("Error reversing complementing sequence", revCompSeqErr)
			return "", revCompSeqErr
		}
		seqFw = string(revCompSeq)
	}

	// QUAL
	qual := string(*read.Quality)
	if flag.IsReverseStrand() {
		qual = string(utils.ReverseBytes(*read.Quality))
	}

	// ATTRIBUTES

	var builder strings.Builder

	// R1 READ
	// QNAME
	builder.WriteString(headerFw)
	builder.WriteString("\t")
	// FLAG
	builder.WriteString(flagStr)
	builder.WriteString("\t")
	// RNAME
	builder.WriteString(rname)
	builder.WriteString("\t")
	// POS
	builder.WriteString(strconv.Itoa(startGenomeFw))
	builder.WriteString("\t")
	// MAPQ
	builder.WriteString(strconv.Itoa(mapqFw))
	builder.WriteString("\t")
	// CIGAR
	builder.WriteString(cigar)
	builder.WriteString("\t")
	// RNEXT
	builder.WriteString("*")
	builder.WriteString("\t")
	// PNEXT
	builder.WriteString("0")
	builder.WriteString("\t")
	// TLEN
	builder.WriteString("0")
	builder.WriteString("\t")
	// SEQ
	builder.WriteString(seqFw)
	builder.WriteString("\t")
	// QUAL
	builder.WriteString(qual)
	builder.WriteString("\n")

	return builder.String(), nil
}

func readPairResultToSamString(
	genomeIndex *index.GenomeIndex,
	readPair *fastq.ReadPair,
	resFw *mapperutils.ReadMatchResult,
	resRv *mapperutils.ReadMatchResult,
) (
	string,
	error,
) {
	if readPair.ReadR2 == nil {
		return readSingleResultToSamString(
			genomeIndex,
			readPair.ReadR1,
			resFw,
		)
	}

	flagFw := sam.Flag{}
	flagRv := sam.Flag{}

	flagFw.SetPaired()
	flagFw.SetFirstInPair()
	if resFw != nil && resRv != nil {
		flagFw.SetProperlyPaired()
	}
	if resFw != nil {
		if !genomeIndex.IsSequenceForward(resFw.SequenceIndex) {
			flagFw.SetReverseStrand()
			flagRv.SetMateReverseStrand()
		}
	} else {
		flagFw.SetUnmapped()
		flagRv.SetMateUnmapped()
	}

	flagRv.SetPaired()
	flagRv.SetSecondInPair()
	if resFw != nil && resRv != nil {
		flagRv.SetProperlyPaired()
	}
	if resRv != nil {
		if !genomeIndex.IsSequenceForward(resRv.SequenceIndex) {
			flagRv.SetReverseStrand()
			flagFw.SetMateReverseStrand()
		}
	} else {
		flagRv.SetUnmapped()
		flagFw.SetMateUnmapped()
	}

	geneStartRelativeToContig := -1
	geneEndRelativeToContig := -1
	if resFw != nil {
		geneStartRelativeToContig = int(genomeIndex.GetSequenceInfo(resFw.SequenceIndex).StartGenomic)
		geneEndRelativeToContig = int(genomeIndex.GetSequenceInfo(resFw.SequenceIndex).EndGenomic)
	} else if resRv != nil {
		geneStartRelativeToContig = int(genomeIndex.GetSequenceInfo(resRv.SequenceIndex).StartGenomic)
		geneEndRelativeToContig = int(genomeIndex.GetSequenceInfo(resRv.SequenceIndex).EndGenomic)
	}

	// QNAME
	// use only the string before any whitespace as the header
	headerFw := strings.Split(readPair.ReadR1.Header, " ")[0]
	headerRv := strings.Split(readPair.ReadR2.Header, " ")[0]

	// FLAG
	flagFwStr := strconv.Itoa(flagFw.Value)
	flagRvStr := strconv.Itoa(flagRv.Value)

	// RNAME
	// the name of the reference sequence (contig) to which the read is aligned
	// it is supposed that both pairs must map to the same contig
	rname := ""
	if resFw != nil {
		rname = genomeIndex.GetSequenceInfo(resFw.SequenceIndex).Contig
	}
	if resRv != nil {
		rname = genomeIndex.GetSequenceInfo(resRv.SequenceIndex).Contig
	}

	// POS
	// the offset is the start of the first region in the matched genome which corresponds to the first
	// mapped position of the read
	startGenomeFw := 0
	if resFw != nil {
		firstRegionGenome, _ := resFw.MatchedGenome.GetFirstRegion()
		offsetFw := firstRegionGenome.Start
		// if the read is mapped to the reverse strand, we need to calculate the offset because the target sequence
		// was reverse complemented. therefore the start position in 5'-3' direction of the original sequence is the
		// end position of the reverse complemented sequence
		if flagFw.IsReverseStrand() {
			geneLength := geneEndRelativeToContig - geneStartRelativeToContig
			lastRegionGenome, _ := resFw.MatchedGenome.GetLastRegion()
			offsetFw = geneLength - lastRegionGenome.End
		}
		// +1 because the sam format is 1-based
		startGenomeFw = geneStartRelativeToContig + offsetFw + 1
	}

	startGenomeRv := 0
	if resRv != nil {
		firstRegionGenome, _ := resRv.MatchedGenome.GetFirstRegion()
		offsetRv := firstRegionGenome.Start
		if flagRv.IsReverseStrand() {
			geneLength := geneEndRelativeToContig - geneStartRelativeToContig
			lastRegionGenome, _ := resRv.MatchedGenome.GetLastRegion()
			offsetRv = geneLength - lastRegionGenome.End
		}
		// +1 because the sam format is 1-based
		startGenomeRv = geneStartRelativeToContig + offsetRv + 1
	}

	// MAPQ
	// https://samtools.github.io/hts-specs/SAMv1.pdf
	// No alignments should be assigned mapping quality 255
	// TODO: implement mapping quality
	mapqFw := 0
	if resFw != nil {
		mapqFw = 254
	}
	mapqRv := 0
	if resRv != nil {
		mapqRv = 254
	}

	// CIGAR
	// https://samtools.github.io/hts-specs/SAMv1.pdf
	// Adjacent CIGAR operations should be different
	cigarFw := "*"
	if resFw != nil {
		var errCigarFw error

		cigarFw, errCigarFw = resFw.GetCigar()
		if errCigarFw != nil {
			logrus.WithFields(logrus.Fields{
				"read":              readPair.ReadR1.Header,
				"length of mapping": resFw.MatchedGenome.Length(),
				"expected length":   len(*readPair.ReadR1.Sequence),
			}).Error("Error getting CIGAR string of FW read: ", errCigarFw)

			// TODO: handle error
			// return "", false
			cigarFw = "*"
		}
	}
	cigarRv := "*"
	if resRv != nil {
		var errCigarRv error

		cigarRv, errCigarRv = resRv.GetCigar()

		if errCigarRv != nil {
			logrus.WithFields(logrus.Fields{
				"read":              readPair.ReadR2.Header,
				"length of mapping": resRv.MatchedGenome.Length(),
				"expected length":   len(*readPair.ReadR2.Sequence),
			}).Error("Error getting CIGAR string of RV read: ", errCigarRv)

			// TODO: handle error
			// return "", false
			cigarRv = "*"
		}
	}

	// RNEXT
	rnextFw := rname
	rnextRv := rname
	if resRv == nil {
		rnextFw = "*"
	}
	if resFw == nil {
		rnextRv = "*"
	}

	// PNEXT
	pnextFw := startGenomeRv
	pnextRv := startGenomeFw
	if resFw == nil {
		pnextRv = 0
	}
	if resRv == nil {
		pnextFw = 0
	}

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
	tlenFw := 0
	tlenRv := 0
	if resFw != nil && resRv != nil {
		tlen := int(math.Abs(float64(startGenomeFw-startGenomeRv))) + len(*readPair.ReadR1.Sequence)
		tlenFw = tlen
		if startGenomeFw > startGenomeRv {
			tlenFw = -1 * tlen
		}
		tlenRv = tlen
		if startGenomeRv >= startGenomeFw {
			tlenRv = -1 * tlen
		}
	}

	// SEQ
	seqFw := string(*readPair.ReadR1.Sequence)
	if flagFw.IsReverseStrand() {
		revCompSeq, revCompSeqErr := utils.ReverseComplementDnaBytes(*readPair.ReadR1.Sequence)
		if revCompSeqErr != nil {
			logrus.WithFields(logrus.Fields{
				"read": readPair.ReadR1.Header,
			}).Error("Error reversing complementing sequence", revCompSeqErr)
			return "", revCompSeqErr
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
			return "", revCompSeqErr
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

	// R1 READ
	if resFw != nil {
		// QNAME
		builder.WriteString(headerFw)
		builder.WriteString("\t")
		// FLAG
		builder.WriteString(flagFwStr)
		builder.WriteString("\t")
		// RNAME
		builder.WriteString(rname)
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

	}

	// R2 READ
	if resRv != nil {
		// QNAME
		builder.WriteString(headerRv)
		builder.WriteString("\t")
		// FLAG
		builder.WriteString(flagRvStr)
		builder.WriteString("\t")
		// RNAME
		builder.WriteString(rname)
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
	}

	return builder.String(), nil
}

func WriteHeader(f *os.File) error {
	header := []string{
		"Gene", "ReadId", "SId", "IsFw", "AltId", "MM",
		"IsFixPoint", "MainAnchorLength", "MainAnchorMM",
		"TotalLeftOptions", "TotalRightOptions",
		"ValidLeftOptions", "ValidRightOptions",
		"LeftFixpointLength", "RightFixpointLength",
		"IsOverhangCorrected",
		"IsGapFill", "IsGapFillOverflow", "GapsFilled", "GapsFilledOverflow",
		"IsSymInErr", "SymInErrLen",
	}
	line := strings.Join(header, "\t") + "\n"
	_, err := f.WriteString(line)
	return err
}
