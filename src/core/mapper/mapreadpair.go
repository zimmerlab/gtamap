package mapper

import (
	"fmt"
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
			fmt.Println("Second pass required for forward read: ", i)
			needSecondPass = true
		}
	}
	for i, resRv := range resultRv {
		if resRv.SecondPass {
			fmt.Println("Second pass required for reverse read: ", i)
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
	builder.WriteString(readPair.ReadR1.Header)
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
	builder.WriteString(readPair.ReadR2.Header)
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
