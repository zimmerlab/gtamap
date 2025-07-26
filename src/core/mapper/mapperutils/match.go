package mapperutils

import (
	"fmt"
	"strconv"
	"strings"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/sirupsen/logrus"
)

type Match struct {
	SequenceIndex int // the index of the sequence in the genome
	FromGenome    int // the start position of the match in the genome
	ToGenome      int // the end position of the match in the genome
	FromRead      int // the start position of the match in the read
	ToRead        int // the end position of the match in the read
	StartGenome   int // the start position of the match in the genome (diagonal)
	Used          bool
}

func (m *Match) String() string {
	return fmt.Sprintf("[%d, %d, %d, %d, %t]", m.FromRead, m.ToRead, m.FromGenome, m.ToGenome, m.Used)
}

type GlobalMatchResult struct {
	MatchesPerSequence []*SequenceMatchResult
}

type SequenceMatchResult struct {
	MatchesPerDiagonal map[int][]*Match
}

type ReadMatchResult struct {
	SequenceIndex   int                        // the index of the sequence in the genome
	MatchedRead     *regionvector.RegionVector // region vector containing the matched positions in the read
	MatchedGenome   *regionvector.RegionVector // region vector containing the matched positions in the genome
	MismatchesRead  []int                      // the positions of the mismatches in the read
	diagonalHandler *DiagonalHandler
	IncompleteMap   bool
	SpliceSitesInfo []bool // corresponds to the number of junctions of match result. is true if junction follows known SpliceSites
}

func (r *ReadMatchResult) HasUnknownSpliceSites() bool {
	if r.SpliceSitesInfo == nil {
		return false
	}
	for _, site := range r.SpliceSitesInfo {
		if !site {
			return true
		}
	}
	return false
}

// GetLargestNachor returns the largest anchor region and also accounts for mms
// WARN: Before calling this, the func NormalizeRegions needs to be called (once)
func (r ReadMatchResult) GetLargestAnchor(introns *regionvector.RegionSet) (regionvector.Region, int, int) {
	if !r.MatchedGenome.HasGaps() {
		return r.MatchedGenome.Regions[0], 0, 0
	}

	usedIndices := make(map[int]float32)
	maxLength := -1
	mainAnchorIndex := 0
	var largestAnchor regionvector.Region
	for len(usedIndices) < len(r.MatchedGenome.Regions) {
		for i, region := range r.MatchedGenome.Regions {

			// have we used this anchor before?
			_, seen := usedIndices[i]
			if seen {
				continue // skip if yes
			}

			if region.Length() > maxLength {
				maxLength = region.Length()
				largestAnchor = region
				mainAnchorIndex = i
			}
		}

		readMainAnchor := r.MatchedRead.Regions[mainAnchorIndex]
		mmsInAnchor := float32(0.0)

		// check how many mm are in main anchor
		for _, mm := range r.MismatchesRead {
			if mm >= readMainAnchor.Start && mm <= readMainAnchor.End {
				mmsInAnchor++
			}
		}

		ratio := mmsInAnchor / float32(largestAnchor.Length())
		if ratio > 0.2 { // main anchor is not supposed to have more than 2% mm
			usedIndices[mainAnchorIndex] = ratio
			maxLength = -1
		} else {
			break
		}
	}

	// all regions had a bad ratio, choose region with smallest ratio
	if len(usedIndices) == len(r.MatchedGenome.Regions) {
		minRatio := float32(2)
		for anchorIndex, ratio := range usedIndices {
			if ratio < minRatio {
				largestAnchor = r.MatchedGenome.Regions[anchorIndex]
				minRatio = ratio
				mainAnchorIndex = anchorIndex
			}
		}
	}

	prevIntron := introns.GetPrevIntron(largestAnchor.Start)
	if prevIntron == nil {
		return largestAnchor, 0, mainAnchorIndex
	}

	return largestAnchor, prevIntron.Rank, mainAnchorIndex // anchor rank == intron rank of prev intron
}

// original which is wrong
func (r *ReadMatchResult) NormalizeRegions() {
	// I. In some cases the regions differ in length but in different positions in genome and read:
	// READ [10,30] [30,40]
	// GENO [10,20] [20,40] -> split into regions of kmer length before normalizing
	kmerLength := int(config.KmerLength())
	for _, region := range r.MatchedGenome.Regions {
		if region.Length() != kmerLength && region.Length()%kmerLength == 0 {
			start := region.Start
			// remove old region
			r.MatchedGenome.RemoveRegion(region.Start, region.End)
			for i := 0; i < region.Length(); i += kmerLength {
				r.MatchedGenome.AddRegionNonOverlappingPanic(start+i, start+i+kmerLength)
			}
		}
	}
	for _, region := range r.MatchedRead.Regions {
		if region.Length() != kmerLength && region.Length()%kmerLength == 0 {
			start := region.Start
			// remove old region
			r.MatchedRead.RemoveRegion(region.Start, region.End)
			for i := 0; i < region.Length(); i += kmerLength {
				r.MatchedRead.AddRegionNonOverlappingPanic(start+i, start+i+kmerLength)
			}
		}
	}

	genomeBlocks := make([]regionvector.Region, 0)
	readBlocks := make([]regionvector.Region, 0)

	startIndex := 0

	hasGap := func(a, b regionvector.Region) bool {
		return a.End < b.Start
	}

	for i := 0; i < len(r.MatchedGenome.Regions)-1; i++ {
		gapInGenome := hasGap(r.MatchedGenome.Regions[i], r.MatchedGenome.Regions[i+1])
		gapInRead := hasGap(r.MatchedRead.Regions[i], r.MatchedRead.Regions[i+1])

		if gapInGenome || gapInRead {
			// account for this case
			// Read [0,10] [10,20] [20,30] [30,40]
			// GEN  [0,10] [10,20] [20,30] [30,36] -> due to either leftNorm or bcause of best split
			//if r.MatchedGenome.Regions[i].Length() < r.MatchedRead.Regions[i].Length() {
			//	padding := r.MatchedRead.Regions[i].Length() - r.MatchedGenome.Regions[i].Length()
			//	r.MatchedRead.Regions[i].End -= padding
			//	r.MatchedRead.Regions[i+1].Start -= padding
			//} else if r.MatchedGenome.Regions[i].Length() > r.MatchedRead.Regions[i].Length() {
			//	padding := r.MatchedGenome.Regions[i].Length() - r.MatchedRead.Regions[i].Length()
			//	r.MatchedRead.Regions[i].End += padding
			//	r.MatchedRead.Regions[i+1].Start += padding
			//}
			genomeBlocks = append(genomeBlocks, regionvector.Region{
				Start: r.MatchedGenome.Regions[startIndex].Start,
				End:   r.MatchedGenome.Regions[i].End,
			})
			readBlocks = append(readBlocks, regionvector.Region{
				Start: r.MatchedRead.Regions[startIndex].Start,
				End:   r.MatchedRead.Regions[i].End,
			})
			startIndex = i + 1
		}
	}

	if startIndex < len(r.MatchedGenome.Regions) {
		genomeBlocks = append(genomeBlocks, regionvector.Region{
			Start: r.MatchedGenome.Regions[startIndex].Start,
			End:   r.MatchedGenome.Regions[len(r.MatchedGenome.Regions)-1].End,
		})
		readBlocks = append(readBlocks, regionvector.Region{
			Start: r.MatchedRead.Regions[startIndex].Start,
			End:   r.MatchedRead.Regions[len(r.MatchedRead.Regions)-1].End,
		})
	}

	r.MatchedGenome.Regions = genomeBlocks
	r.MatchedRead.Regions = readBlocks
}

func (r *ReadMatchResult) Copy() *ReadMatchResult {
	var dhCopy *DiagonalHandler
	if r.diagonalHandler != nil {
		dhCopy = r.diagonalHandler.Copy()
	}

	return &ReadMatchResult{
		SequenceIndex:   r.SequenceIndex,
		MatchedRead:     r.MatchedRead.Copy(),
		MatchedGenome:   r.MatchedGenome.Copy(),
		MismatchesRead:  append([]int{}, r.MismatchesRead...),
		diagonalHandler: dhCopy,
	}
}

func AssignReadMatchResults(fwMappings []*ReadMatchResult, rvMappings []*ReadMatchResult) (map[int][]*ReadMatchResult, map[int][]*ReadMatchResult, map[int]struct{}) {
	fwMapPerSeqIndex := make(map[int][]*ReadMatchResult)
	rvMapPerSeqIndex := make(map[int][]*ReadMatchResult)
	mappedRegionIds := make(map[int]struct{})

	// accumulate mapping results
	for _, mapping := range fwMappings {
		regionIndex := mapping.SequenceIndex / 2 // maps back to the main sequence index
		fwMapPerSeqIndex[regionIndex] = append(fwMapPerSeqIndex[regionIndex], mapping)
		mappedRegionIds[regionIndex] = struct{}{}
	}

	for _, mapping := range rvMappings {
		regionIndex := mapping.SequenceIndex / 2
		rvMapPerSeqIndex[regionIndex] = append(rvMapPerSeqIndex[regionIndex], mapping)
		mappedRegionIds[regionIndex] = struct{}{}
	}

	return fwMapPerSeqIndex, rvMapPerSeqIndex, mappedRegionIds
}

// receives fw and rv matches of one seqID. Returns possible combinations of fw/rv.
// E.g. fw -> 25 and rv -> 25 doesnt work since they need to map to separate strands etc
// Currently returns the best combination with less or equal amount of mm the provided maxMismatches param
func GetBestPossibleMappingCombination(fwMatches []*ReadMatchResult, rvMatches []*ReadMatchResult) *ValidReadPairCombination {
	var bestCombination *ValidReadPairCombination
	maxMismatches := config.MaxConfMm

	for i := 0; i < len(fwMatches); i++ {
		fwMatch := fwMatches[i]
		for j := 0; j < len(rvMatches); j++ {
			rvMatch := rvMatches[j]
			if fwMatch.SequenceIndex-1 == rvMatch.SequenceIndex || fwMatch.SequenceIndex == rvMatch.SequenceIndex-1 {
				// dont allow IncompleteMaps
				if fwMatch.IncompleteMap || rvMatch.IncompleteMap {
					continue
				}
				currMM := len(fwMatch.MismatchesRead) + len(rvMatch.MismatchesRead)
				if currMM < maxMismatches {
					// skip possible combinations if they have short diagonals
					// (like a diag of length 10 in the end; these cases are likely wrong and should not be classified as confident map)
					if !(hasLongDiagonals(rvMatch) && hasLongDiagonals(fwMatch)) {
						continue
					}
					maxMismatches = currMM
					if bestCombination == nil {
						bestCombination = &ValidReadPairCombination{
							Fw:            fwMatch,
							Rv:            rvMatch,
							NumMismatches: currMM,
						}
					} else {
						maxMismatches = currMM
						bestCombination.Fw = fwMatch
						bestCombination.Rv = rvMatch
						bestCombination.NumMismatches = currMM
					}
				}
			}
		}
	}

	return bestCombination
}

// returns false as soon as one region is smaller than 2*kmerlength
// iterates over gaps and checks lengths of region before and after gaps
func hasLongDiagonals(mapping *ReadMatchResult) bool {
	mapping.NormalizeRegions()
	for _, region := range mapping.MatchedGenome.Regions {
		if region.Length() < config.MinConfAnchorLength {
			return false
		}
	}
	return true
	// // used to keep track of the read position for the next gap
	// readGapPos := 0
	// // returns the index of the first region after which a gap occurs (-1 if no gap)
	// indexRegionBeforeGap := mapping.MatchedGenome.GetGapIndexAfterPos(readGapPos)
	//
	// startIndex := 0
	//
	// // loop through all gaps in the read (-1 means there is no more gap)
	// for indexRegionBeforeGap > -1 {
	// 	if mapping.MatchedGenome.Regions[indexRegionBeforeGap].End-mapping.MatchedGenome.Regions[startIndex].Start < config.MinConfAnchorLength {
	// 		return false
	// 	}
	// 	startIndex = indexRegionBeforeGap + 1
	// 	readGapPos = mapping.MatchedGenome.Regions[indexRegionBeforeGap].End + 1
	// 	indexRegionBeforeGap = mapping.MatchedGenome.GetGapIndexAfterPos(readGapPos)
	// }
	//
	// if indexRegionBeforeGap == -1 && startIndex != 0 {
	// 	if mapping.MatchedGenome.GetLastRegion().End-mapping.MatchedGenome.Regions[startIndex].Start < config.MinConfAnchorLength {
	// 		return false
	// 	}
	// }
	// return true
}

func (r *ReadMatchResult) GetCigar() (string, error) {
	var builder strings.Builder

	isForwardStrand := r.SequenceIndex == 0

	if isForwardStrand {

		// number of cumulative aligned positions (matches or mismatches) between the read and the genome
		numMatchesSum := 0

		if len(r.MatchedGenome.Regions) != len(r.MatchedRead.Regions) {
			r.SyncRegions()
		}

		// the regions in MatchedGenome and MatchedRead have the same dimensions
		for i := 0; i < len(r.MatchedGenome.Regions); i++ {

			// each pair of regions represent a match
			numMatchesSum += r.MatchedGenome.Regions[i].Length()

			// if this is the last match, add the number of matches
			if i == len(r.MatchedGenome.Regions)-1 {
				builder.WriteString(strconv.Itoa(numMatchesSum))
				builder.WriteString("M")
				break
			}

			gapInGenome := r.MatchedGenome.Regions[i].End < r.MatchedGenome.Regions[i+1].Start
			gapInRead := r.MatchedRead.Regions[i].End < r.MatchedRead.Regions[i+1].Start

			if gapInGenome && gapInRead {
				logrus.WithFields(logrus.Fields{
					"read":   r.MatchedRead,
					"genome": r.MatchedGenome,
				}).Warn("Gap in genome and read at the same time")

				return "", fmt.Errorf("gap in genome and read at the same time")
			}

			if !gapInGenome && !gapInRead {
				continue
			}

			builder.WriteString(strconv.Itoa(numMatchesSum))
			builder.WriteString("M")
			numMatchesSum = 0

			// intron or deletion
			if gapInGenome {
				// number of skipped bases in the reference
				numSkipped := r.MatchedGenome.Regions[i+1].Start - r.MatchedGenome.Regions[i].End

				builder.WriteString(strconv.Itoa(numSkipped))

				// determine whether the gap is an intron or a deletion
				if numSkipped < config.IntronLengthMin() {
					builder.WriteString("D")
				} else {
					builder.WriteString("N")
				}
			}

			// insertion
			if gapInRead {
				// number of skipped bases in the read
				numSkipped := r.MatchedRead.Regions[i+1].Start - r.MatchedRead.Regions[i].End

				builder.WriteString(strconv.Itoa(numSkipped))
				builder.WriteString("I")
			}
		}

	} else {

		// reverse the order of the regions in MatchedGenome
		// this is because the read is on the reverse strand

		numMatchesSum := 0

		if len(r.MatchedGenome.Regions) != len(r.MatchedRead.Regions) {
			r.SyncRegions()
		}

		// the regions in MatchedGenome and MatchedRead have the same dimensions
		for i := len(r.MatchedGenome.Regions) - 1; i >= 0; i-- {

			numMatchesSum += r.MatchedGenome.Regions[i].Length()

			// if this is the last match, add the number of matches
			if i == 0 {
				builder.WriteString(strconv.Itoa(numMatchesSum))
				builder.WriteString("M")
				break
			}

			gapInGenome := r.MatchedGenome.Regions[i].Start > r.MatchedGenome.Regions[i-1].End
			gapInRead := r.MatchedRead.Regions[i].Start > r.MatchedRead.Regions[i-1].End

			if gapInGenome && gapInRead {
				logrus.WithFields(logrus.Fields{
					"read":   r.MatchedRead,
					"genome": r.MatchedGenome,
				}).Warn("Gap in genome and read at the same time")

				return "", fmt.Errorf("gap in genome and read at the same time")
			}

			if !gapInGenome && !gapInRead {
				continue
			}

			builder.WriteString(strconv.Itoa(numMatchesSum))
			builder.WriteString("M")
			numMatchesSum = 0

			// intron or deletion
			if gapInGenome {
				// number of skipped bases in the reference
				numSkipped := r.MatchedGenome.Regions[i].Start - r.MatchedGenome.Regions[i-1].End

				builder.WriteString(strconv.Itoa(numSkipped))

				// determine whether the gap is an intron or a deletion
				if numSkipped < config.IntronLengthMin() {
					builder.WriteString("D")
				} else {
					builder.WriteString("N")
				}
			}

			// insertion
			if gapInRead {
				// number of skipped bases in the read
				numSkipped := r.MatchedRead.Regions[i].Start - r.MatchedRead.Regions[i-1].End

				builder.WriteString(strconv.Itoa(numSkipped))
				builder.WriteString("I")
			}

		}
	}

	return builder.String(), nil
}

// SyncRegions : if a ReadMatchResult was previously incomplete and had a readRegions rv like so
// [[50,60], [60, 70] , ..., [140,150]], after the remap happened, it looks like this
// [[0,50], [50,60], [60, 70] , ..., [140,150]], after the remap happened, it looks like this
// now the dimensions don't always add up -> we have to normalize the readRegions to match the genomeRegions
// E.g if genomeRegions contains three regions (with a total length of 150), the read regions also should only
// contain 3 regions
func (r *ReadMatchResult) SyncRegions() {
	adaptedReadRegions := make([]regionvector.Region, len(r.MatchedGenome.Regions))
	originalStart, ok := r.MatchedRead.GetFirstRegion()
	if ok {
		for i, region := range r.MatchedGenome.Regions {
			startRead, err := regionvector.GenomicCoordToReadCoord(originalStart.Start, region.Start, r.MatchedGenome.Regions)
			if err != nil {
				logrus.Errorf("Error while converting genomic coord to read coord")
				logrus.Fatal(err)
			}
			stopRead, err := regionvector.GenomicCoordToReadCoord(originalStart.Start, region.End, r.MatchedGenome.Regions)
			if err != nil {
				logrus.Errorf("Error while converting genomic coord to read coord")
				logrus.Fatal(err)
			}
			adaptedReadRegions[i] = regionvector.Region{Start: startRead, End: stopRead}
		}
		r.MatchedRead.Regions = adaptedReadRegions
	}
}

// holds all relevant mapping information of potentially several hits per readpair
// bundled into one type
type ReadPairMatchResults struct {
	ReadPair *fastq.ReadPair
	Fw       []*ReadMatchResult
	Rv       []*ReadMatchResult
}

type ValidReadPairCombination struct {
	Fw            *ReadMatchResult
	Rv            *ReadMatchResult
	NumMismatches int
}

// hold information for main seq id
type TargetAnnotation struct {
	PreferedStrand int                             // 0 -> + (fw+ and rv-); 1 -> - (fw- and rv+)
	Confidence     float32                         // percentage of reads contributing to PreferedStrand
	Introns        map[int]*regionvector.RegionSet // maps to sub sequence index, meaning each gene has two slices of introns 0 -> plusOrintation 1 -> minusOrientation
	// zero bases coords, start incl, stop excl
}

func (t TargetAnnotation) String() string {
	var sb strings.Builder
	strand := "+"
	if t.PreferedStrand == 1 {
		strand = "-"
	}

	sb.WriteString(fmt.Sprintf("Preferred Strand: [%s] ", strand))
	sb.WriteString(fmt.Sprintf("(Confidence: %.2f%%)", t.Confidence*100))

	for orientation, introns := range t.Introns {
		orientationLabel := "+"
		if orientation == 1 {
			orientationLabel = "-"
		}
		sb.WriteString(fmt.Sprintf("\n(%s):", orientationLabel))
		for _, intron := range introns.Regions {
			sb.WriteString(intron.String())
		}
	}

	return sb.String()
}

func (t TargetAnnotation) LogInfo() {
	strand := "+"
	if t.PreferedStrand == 1 {
		strand = "-"
	}

	logrus.WithFields(logrus.Fields{
		"preferred strand": strand,
		"confidence score": t.Confidence * 100,
	}).Info("Done with Annotation")

	for orientation, introns := range t.Introns {
		orientationLabel := "[+]"
		if orientation == 1 {
			orientationLabel = "[-]"
		}
		logrus.WithFields(logrus.Fields{
			"orientation": orientationLabel,
		}).Info("Inferred Introns of")

		for _, intron := range introns.Regions {
			logrus.WithFields(logrus.Fields{
				"Intron": intron.String(),
			}).Info()
		}
	}
}

func (r ReadPairMatchResults) String() string {
	var builder strings.Builder
	builder.Write([]byte("ReadPairR1 Header: "))
	builder.Write([]byte(r.ReadPair.ReadR1.Header))
	builder.Write([]byte("\n"))
	builder.Write([]byte("  <== FW MAPPINGS ==>"))
	builder.Write([]byte("\n"))
	for _, mapping := range r.Fw {
		builder.Write([]byte("\t SeqIndex: "))
		seqIndex := strconv.Itoa(mapping.SequenceIndex)
		builder.WriteString(seqIndex)
		builder.WriteString("\n")
		builder.WriteString("\t GENOME -> ")
		builder.WriteString(mapping.MatchedGenome.String())
		builder.WriteString("\n")
		builder.WriteString("\t READ   -> ")
		builder.WriteString(mapping.MatchedRead.String())
		builder.WriteString("\n")
		builder.WriteString("\t MISMAT -> ")
		ints := mapping.MismatchesRead
		strs := make([]string, len(ints))
		for i, v := range ints {
			strs[i] = strconv.Itoa(v)
		}
		builder.WriteString(strings.Join(strs, ","))
		builder.WriteString("\n")
	}
	builder.Write([]byte("  <== RV MAPPINGS ==>"))
	builder.Write([]byte("\n"))
	for _, mapping := range r.Rv {
		builder.Write([]byte("\t SeqIndex: "))
		seqIndex := strconv.Itoa(mapping.SequenceIndex)
		builder.WriteString(seqIndex)
		builder.WriteString("\n")
		builder.WriteString("\t GENOME -> ")
		builder.WriteString(mapping.MatchedGenome.String())
		builder.WriteString("\n")
		builder.WriteString("\t READ   -> ")
		builder.WriteString(mapping.MatchedRead.String())
		builder.WriteString("\n")
		builder.WriteString("\t MISMAT -> ")
		ints := mapping.MismatchesRead
		strs := make([]string, len(ints))
		for i, v := range ints {
			strs[i] = strconv.Itoa(v)
		}
		builder.WriteString(strings.Join(strs, ","))
		builder.WriteString("\n")
	}
	return builder.String()
}
