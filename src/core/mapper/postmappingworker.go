package mapper

import (
	"fmt"
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

func PostMappingWorker(mappedReadPairChan <-chan *ReadPairMatchResults, wg *sync.WaitGroup, outputChan chan<- string) {
	logrus.Info("Started accumulating mapped readpairs")
	defer wg.Done()

	// I. Extract and sort mappring readpairs into map: mainSeqIndex -> [mappings]
	resultsPerSeqIndex := make(map[int][]*ReadPairMatchResults)
	for mappedReadPair := range mappedReadPairChan {
		// TODO: Handle multimappings
		if len(mappedReadPair.Fw) > 1 || len(mappedReadPair.Rv) > 1 {

			if len(mappedReadPair.Fw) != len(mappedReadPair.Rv) {
			}

			logrus.Fatalf("Multimapping readpair found. Currently not handled: Read ID: %s", mappedReadPair.ReadPair.ReadR1.Header)
		}

		if mappedReadPair.Fw[0].SequenceIndex/2 != mappedReadPair.Rv[0].SequenceIndex/2 {
			logrus.Infof("Mates of Readpair don't agree on seq index %s", mappedReadPair.ReadPair.ReadR1.Header)
		}

		// For now I simply take the first Fw and Rv
		// projects 0,1 -> 0; 2,3 -> 1; 4,5 -> 2 to map to main taret seq id
		parentSeqIndex := mappedReadPair.Fw[0].SequenceIndex / 2

		resultsPerSeqIndex[parentSeqIndex] = append(resultsPerSeqIndex[parentSeqIndex], mappedReadPair)
	}
	logrus.Info("Finished accumulating mapped readpairs")

	// II. extract coverage
	// Here we assume that there are already only one fw and rv per mappedREadPair
	// I will have to create a new struct which only holds one fw/rv but other than that is the exact same
	// as ReadPairMatchResults
	multiMapsPerSeqIndex := make(map[string][]*ReadPairMatchResults)
	for seqIndex, mappedReadPairs := range resultsPerSeqIndex {
		coverageSliceFw, coverageSliceRv := getCoverageSlices(mappedReadPairs, len((*mappedReadPairs[0].Index.Sequences[seqIndex])))
		highCovRegionsfw := getHighCoverageRegons(coverageSliceFw)
		highCovRegionsrv := getHighCoverageRegons(coverageSliceRv)
		multimaps := getMultimappingReads(mappedReadPairs, highCovRegionsfw, highCovRegionsrv)
		targetRegionId := mappedReadPairs[0].Index.GetSequenceInfo(seqIndex).GeneId
		multiMapsPerSeqIndex[targetRegionId] = multimaps
	}
	for seqIndex, ambiguousReadPairs := range multiMapsPerSeqIndex {
		for _, ambiguousReadPair := range ambiguousReadPairs {
			alternativeMappingFw, _ := MapRead(ambiguousReadPair.ReadPair.ReadR1, ambiguousReadPair.Index.ParalogRegions[seqIndex])
			alternativeMappingRv, _ := MapRead(ambiguousReadPair.ReadPair.ReadR2, ambiguousReadPair.Index.ParalogRegions[seqIndex])
			fmt.Println(alternativeMappingFw)
			fmt.Println(alternativeMappingRv)

		}
	}

	// WRITE TO SAM
	for _, mappedReadPairs := range resultsPerSeqIndex {
		for _, mappedReadPair := range mappedReadPairs {
			fwRes := mappedReadPair.Fw[0]
			rvRes := mappedReadPair.Rv[0]
			output, isMapping := FormatMappedReadPairToSAM(fwRes, rvRes, mappedReadPair.ReadPair, mappedReadPair.Index)
			if isMapping {
				outputChan <- output
			}
		}
	}
}

// holds all relevant mapping information of potentially several hits per readpair
// bundled into one type
type ReadPairMatchResults struct {
	ReadPair *fastq.ReadPair
	Fw       []mapperutils.ReadMatchResult
	Rv       []mapperutils.ReadMatchResult
	Index    *index.GenomeIndex
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

func getCoverageSlice(mappedReadPairs []*ReadPairMatchResults, geneLength int) []int {
	coverage := make([]int, geneLength)
	for _, mappedReadPair := range mappedReadPairs {
		// forward reads
		for _, region := range mappedReadPair.Fw[0].MatchedGenome.Regions {
			if mappedReadPair.Index.IsSequenceForward(mappedReadPair.Fw[0].SequenceIndex) {
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
			if mappedReadPair.Index.IsSequenceForward(mappedReadPair.Rv[0].SequenceIndex) {
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

func getMultimappingReads(mappedReadPairs []*ReadPairMatchResults, fwIntervals [][2]int, rvIntervals [][2]int) []*ReadPairMatchResults {
	multimappings := make([]*ReadPairMatchResults, 0)
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

func getCoverageSlices(mappedReadPairs []*ReadPairMatchResults, geneLength int) ([]int, []int) {
	coverageFw := make([]int, geneLength)
	coverageRv := make([]int, geneLength)

	for _, mappedReadPair := range mappedReadPairs {
		// fw reads
		for _, region := range mappedReadPair.Fw[0].MatchedGenome.Regions {
			if mappedReadPair.Index.IsSequenceForward(mappedReadPair.Fw[0].SequenceIndex) {
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
			if mappedReadPair.Index.IsSequenceForward(mappedReadPair.Rv[0].SequenceIndex) {
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
