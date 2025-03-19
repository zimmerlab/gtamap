package mapper

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/core/timer"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/KleinSamuel/gtamap/src/utils"
	"strconv"
	"strings"
)

func MapReadPair(readPair *fastq.ReadPair, genomeIndex *index.GenomeIndex,
	secondPassChan *mapperutils.SecondPassChannel,
	timerChannel chan<- *timer.Timer) (string, bool) {

	keepFw := Filter(readPair.ReadR1.Sequence, genomeIndex)
	keepRw := Filter(readPair.ReadR2.Sequence, genomeIndex)

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

	resFw := resultFw[0]
	resRv := resultRv[0]

	var builder strings.Builder

	// QNAME
	builder.WriteString(readPair.ReadR1.Header)
	builder.WriteString("\t")
	// FLAG
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// RNAME
	builder.WriteString(genomeIndex.GetSequenceHeader(resFw.SequenceIndex))
	builder.WriteString("\t")
	// POS
	startRelativeFw := resFw.MatchedGenome.GetFirstRegion().Start
	// TODO: use actual start position of sequence from index
	startGenomeFw := 45884425 + startRelativeFw

	builder.WriteString(strconv.Itoa(startGenomeFw))
	builder.WriteString("\t")
	// MAPQ
	builder.WriteString(strconv.Itoa(255))
	builder.WriteString("\t")
	// CIGAR
	builder.WriteString(resFw.GetCigar())
	builder.WriteString("\t")
	// PNEXT
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// RNEXT
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// TLEN
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// SEQ
	if genomeIndex.IsSequenceForward(resFw.SequenceIndex) {
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
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// RNAME
	builder.WriteString(genomeIndex.GetSequenceHeader(resRv.SequenceIndex))
	builder.WriteString("\t")
	// POS
	startRelativeRw := resRv.MatchedGenome.GetLastRegion().End
	startGenomeRw := 45903174 - startRelativeRw + 1
	builder.WriteString(strconv.Itoa(startGenomeRw))
	builder.WriteString("\t")
	// MAPQ
	builder.WriteString(strconv.Itoa(255))
	builder.WriteString("\t")
	// CIGAR
	builder.WriteString(resRv.GetCigar())
	builder.WriteString("\t")
	// PNEXT
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// RNEXT
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// TLEN
	builder.WriteString(strconv.Itoa(0))
	builder.WriteString("\t")
	// SEQ
	if genomeIndex.IsSequenceForward(resRv.SequenceIndex) {
		builder.WriteString(string(*readPair.ReadR2.Sequence))
	} else {
		builder.WriteString(string(utils.ReverseComplementDnaBytes(*readPair.ReadR2.Sequence)))
	}
	builder.WriteString("\t")
	// QUAL
	builder.WriteString(string(*readPair.ReadR2.Quality))
	builder.WriteString("\n")

	//resString := readPair.ReadR1.Header + "_fw\t"
	//for _, resFw := range resultFw {
	//	startRelative := resFw.MatchedGenome.GetFirstRegion().Start
	//	startGenome := 45884425 + startRelative
	//	resString += strconv.Itoa(resFw.SequenceIndex) + ":" + strconv.Itoa(startRelative) + "|" + strconv.Itoa(startGenome) + ":" + strconv.Itoa(resFw.MatchedRead.Length()) + " "
	//}
	//resString = resString + "\n"
	//
	//resString += readPair.ReadR2.Header + "_rw\t"
	//for _, resRw := range resultRw {
	//	startRelative := resRw.MatchedGenome.GetLastRegion().End
	//	startGenome := 45903174 - startRelative + 1
	//	resString += strconv.Itoa(resRw.SequenceIndex) + ":" + strconv.Itoa(startRelative) + "|" + strconv.Itoa(startGenome) + ":" + strconv.Itoa(resRw.MatchedRead.Length()) + " "
	//}
	//resString = resString + "\n"

	return builder.String(), true
}
