package thirdpass

import (
	"math"
	"strconv"
	"strings"
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
)

func ThirdPassWorker(thirdPassChan *ThirdPassChannel, wgThirdPass *sync.WaitGroup, outputChan chan<- string, index *index.GenomeIndex) {
	defer wgThirdPass.Done()
	logrus.Info("Started third pass")

	for {

		task, ok := thirdPassChan.Receive()
		if !ok {
			break
		}
		// logrus.Infof("Third pass: %s", task.ReadPairId)

		var builder strings.Builder

		for i := 0; i < len(task.TargetInfo.Fw); i++ {
			for j := 0; j < len(task.TargetInfo.Rv); j++ {
				s, isOk := readPairResultToSamString(index, task.TargetInfo.ReadPair, task.TargetInfo.Fw[i], task.TargetInfo.Rv[j])
				if !isOk {
					continue
				}
				builder.WriteString(s)
			}
		}
		outputChan <- builder.String()
	}
	logrus.Info("Done with third pass")
}

func readPairResultToSamString(genomeIndex *index.GenomeIndex, readPair *fastq.ReadPair,
	resFw *mapperutils.ReadMatchResult, resRv *mapperutils.ReadMatchResult,
) (string, bool) {
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
	if resFw != nil {
		geneStartRelativeToContig = int(genomeIndex.GetSequenceInfo(resFw.SequenceIndex).StartGenomic)
	}
	geneEndRelativeToContig := -1
	if resRv != nil {
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
		offsetFw := resFw.MatchedGenome.GetFirstRegion().Start
		// if the read is mapped to the reverse strand, we need to calculate the offset because the target sequence
		// was reverse complemented. therefore the start position in 5'-3' direction of the original sequence is the
		// end position of the reverse complemented sequence
		if flagFw.IsReverseStrand() {
			geneLength := geneEndRelativeToContig - geneStartRelativeToContig
			offsetFw = geneLength - resFw.MatchedGenome.GetLastRegion().End
		}
		// +1 because the sam format is 1-based
		startGenomeFw = geneStartRelativeToContig + offsetFw + 1
	}

	startGenomeRv := 0
	if resRv != nil {
		offsetRv := resRv.MatchedGenome.GetFirstRegion().Start
		if flagRv.IsReverseStrand() {
			geneLength := geneEndRelativeToContig - geneStartRelativeToContig
			offsetRv = geneLength - resRv.MatchedGenome.GetLastRegion().End
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
				"read": readPair.ReadR1.Header,
			}).Error("Error getting CIGAR string", errCigarFw)

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
				"read": readPair.ReadR2.Header,
			}).Error("Error getting CIGAR string", errCigarRv)

			// TODO: handle error
			// return "", false
			cigarRv = "*"
		}
	}

	// RNEXT
	rnextFw := rname
	rnextRv := rname

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

	// R2 READ

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

	return builder.String(), true
}

func getCoverageSlice(mappedReadPairs []*mapperutils.ReadPairMatchResults, geneLength int, index *index.GenomeIndex) []int {
	coverage := make([]int, geneLength)
	for _, mappedReadPair := range mappedReadPairs {
		// forward reads
		for _, region := range mappedReadPair.Fw[0].MatchedGenome.Regions {
			if index.IsSequenceForward(mappedReadPair.Fw[0].SequenceIndex) {
				for pos := region.Start; pos < region.End; pos++ {
					coverage[pos]++
				}
			} else {
				// reverse strand reads need complementary position
				// as done in FormatMappedReadPairToSAM
				for pos := region.Start; pos < region.End; pos++ {
					complementaryPos := geneLength - pos - 1
					coverage[complementaryPos]++
				}
			}
		}

		// reverse reads
		for _, region := range mappedReadPair.Rv[0].MatchedGenome.Regions {
			if index.IsSequenceForward(mappedReadPair.Rv[0].SequenceIndex) {
				for pos := region.Start; pos < region.End; pos++ {
					coverage[pos]++
				}
			} else {
				// get compl pos
				for pos := region.Start; pos < region.End; pos++ {
					complementaryPos := geneLength - pos - 1
					coverage[complementaryPos]++
				}
			}
		}
	}
	return coverage
}

func getHighCoverageRegons(coverageSlice []int) [][2]int {
	highCov := make([][2]int, 0)
	var deltaOne int
	var deltaTwo int
	foundStart := false
	var start int
	var covAtStart int
	for i := 0; i < len(coverageSlice)-8; i++ {
		deltaOne = coverageSlice[i+1] - coverageSlice[i]
		deltaTwo = coverageSlice[i+8] - coverageSlice[i]

		if !foundStart && deltaOne > 30 && deltaTwo > 200 {
			foundStart = true
			start = i
			covAtStart = coverageSlice[i]
		}

		if foundStart && coverageSlice[i] < covAtStart {
			foundStart = false
			highCov = append(highCov, [2]int{start, i})
		}

	}
	return highCov
}

func getMultimappingReads(mappedReadPairs []*mapperutils.ReadPairMatchResults, fwIntervals [][2]int, rvIntervals [][2]int) []*mapperutils.ReadPairMatchResults {
	multimappings := make([]*mapperutils.ReadPairMatchResults, 0)
main:
	for _, mappedReadPair := range mappedReadPairs {
		// check fw
		for _, interval := range fwIntervals {
			// append multimapp only if #mm > 0
			if spansInterval(mappedReadPair.Fw[0].MatchedGenome, interval) && len(mappedReadPair.Fw[0].MismatchesRead) != 0 {
				multimappings = append(multimappings, mappedReadPair)
				continue main
			}
		}

		// check rv
		for _, interval := range rvIntervals {
			// append multimapp only if #mm > 0
			if spansInterval(mappedReadPair.Rv[0].MatchedGenome, interval) && len(mappedReadPair.Rv[0].MismatchesRead) != 0 {
				multimappings = append(multimappings, mappedReadPair)
				continue main
			}
		}
	}
	return multimappings
}

func spansInterval(regions *regionvector.RegionVector, interval [2]int) bool {
	for _, region := range regions.Regions {
		if (region.Start <= interval[0] && region.End > interval[0]) || // region starts before interval and extends into it
			(region.Start < interval[1] && region.End >= interval[1]) || // region ends after interval and starts within it
			(region.Start >= interval[0] && region.End <= interval[1]) { // region is completely contained within interval
			return true
		}
	}
	return false
}

func getCoverageSlices(mappedReadPairs []*mapperutils.ReadPairMatchResults, geneLength int, index *index.GenomeIndex) ([]int, []int) {
	coverageFw := make([]int, geneLength)
	coverageRv := make([]int, geneLength)

	for _, mappedReadPair := range mappedReadPairs {
		// fw reads
		for _, region := range mappedReadPair.Fw[0].MatchedGenome.Regions {
			if index.IsSequenceForward(mappedReadPair.Fw[0].SequenceIndex) {
				for pos := region.Start; pos < region.End; pos++ {
					coverageFw[pos]++
				}
			} else {
				for pos := region.Start; pos < region.End; pos++ {
					coverageRv[pos]++
				}
			}
		}

		// rv reads
		for _, region := range mappedReadPair.Rv[0].MatchedGenome.Regions {
			if index.IsSequenceForward(mappedReadPair.Rv[0].SequenceIndex) {
				for pos := region.Start; pos < region.End; pos++ {
					coverageFw[pos]++
				}
			} else {
				for pos := region.Start; pos < region.End; pos++ {
					coverageRv[pos]++
				}
			}
		}
	}

	return coverageFw, coverageRv
}

func getTotalCoverage(coverageFw, coverageRv []int) []int {
	// this can be used to get the total coverage of fw and rv by mirrowing rv coords
	totalCoverage := make([]int, len(coverageFw))
	for i := 0; i < len(coverageFw); i++ {
		totalCoverage[i] = coverageFw[i] + coverageRv[len(coverageFw)-1-i]
	}
	return totalCoverage
}
