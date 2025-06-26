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

func (r ReadMatchResult) HasUnknownSpliceSites() bool {
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

func (m ReadMatchResult) Copy() *ReadMatchResult {
	var dhCopy *DiagonalHandler
	if m.diagonalHandler != nil {
		dhCopy = m.diagonalHandler.Copy()
	}

	return &ReadMatchResult{
		SequenceIndex:   m.SequenceIndex,
		MatchedRead:     m.MatchedRead.Copy(),
		MatchedGenome:   m.MatchedGenome.Copy(),
		MismatchesRead:  append([]int{}, m.MismatchesRead...),
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
func GetBestPossibleMappingCombination(fwMatches []*ReadMatchResult, rvMatches []*ReadMatchResult, maxMismatches int) *ValidReadPairCombination {
	var bestCombination *ValidReadPairCombination

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
	// used to keep track of the read position for the next gap
	readGapPos := 0
	// returns the index of the first region after which a gap occurs (-1 if no gap)
	indexRegionBeforeGap := mapping.MatchedGenome.GetGapIndexAfterPos(readGapPos)

	startIndex := 0

	// loop through all gaps in the read (-1 means there is no more gap)
	for indexRegionBeforeGap > -1 {
		if mapping.MatchedGenome.Regions[indexRegionBeforeGap].End-mapping.MatchedGenome.Regions[startIndex].Start < int(config.KmerLength())*2 {
			return false
		}
		startIndex = indexRegionBeforeGap + 1
		readGapPos = mapping.MatchedGenome.Regions[indexRegionBeforeGap].End + 1
		indexRegionBeforeGap = mapping.MatchedGenome.GetGapIndexAfterPos(readGapPos)
	}

	if indexRegionBeforeGap == -1 && startIndex != 0 {
		if mapping.MatchedGenome.GetLastRegion().End-mapping.MatchedGenome.Regions[startIndex].Start < int(config.KmerLength())*2 {
			return false
		}
	}
	return true
}

func (m ReadMatchResult) GetCigar() (string, error) {
	var builder strings.Builder
	if m.MatchedGenome.Length() != m.MatchedRead.Length() {
		return "", fmt.Errorf("length of matched genome unequal to length matched read: Genome: %d vs Read: %d", m.MatchedGenome.Length(), m.MatchedRead.Length())
	}

	isForwardStrand := m.SequenceIndex == 0

	if isForwardStrand {

		// number of cumulative aligned positions (matches or mismatches) between the read and the genome
		numMatchesSum := 0

		if len(m.MatchedGenome.Regions) != len(m.MatchedRead.Regions) {
			m.normalizeRegions()
		}

		// the regions in MatchedGenome and MatchedRead have the same dimensions
		for i := 0; i < len(m.MatchedGenome.Regions); i++ {

			// each pair of regions represent a match
			numMatchesSum += m.MatchedGenome.Regions[i].Length()

			// if this is the last match, add the number of matches
			if i == len(m.MatchedGenome.Regions)-1 {
				builder.WriteString(strconv.Itoa(numMatchesSum))
				builder.WriteString("M")
				break
			}
			if m.MatchedRead == nil {
				println()
			}

			gapInGenome := m.MatchedGenome.Regions[i].End < m.MatchedGenome.Regions[i+1].Start
			gapInRead := m.MatchedRead.Regions[i].End < m.MatchedRead.Regions[i+1].Start

			if gapInGenome && gapInRead {
				logrus.WithFields(logrus.Fields{
					"read":   m.MatchedRead,
					"genome": m.MatchedGenome,
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
				numSkipped := m.MatchedGenome.Regions[i+1].Start - m.MatchedGenome.Regions[i].End

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
				numSkipped := m.MatchedRead.Regions[i+1].Start - m.MatchedRead.Regions[i].End

				builder.WriteString(strconv.Itoa(numSkipped))
				builder.WriteString("I")
			}
		}

	} else {

		// reverse the order of the regions in MatchedGenome
		// this is because the read is on the reverse strand

		numMatchesSum := 0

		if len(m.MatchedGenome.Regions) != len(m.MatchedRead.Regions) {
			m.normalizeRegions()
		}

		// the regions in MatchedGenome and MatchedRead have the same dimensions
		for i := len(m.MatchedGenome.Regions) - 1; i >= 0; i-- {

			numMatchesSum += m.MatchedGenome.Regions[i].Length()

			// if this is the last match, add the number of matches
			if i == 0 {
				builder.WriteString(strconv.Itoa(numMatchesSum))
				builder.WriteString("M")
				break
			}

			gapInGenome := m.MatchedGenome.Regions[i].Start > m.MatchedGenome.Regions[i-1].End
			gapInRead := m.MatchedRead.Regions[i].Start > m.MatchedRead.Regions[i-1].End

			if gapInGenome && gapInRead {
				logrus.WithFields(logrus.Fields{
					"read":   m.MatchedRead,
					"genome": m.MatchedGenome,
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
				numSkipped := m.MatchedGenome.Regions[i].Start - m.MatchedGenome.Regions[i-1].End

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
				numSkipped := m.MatchedRead.Regions[i].Start - m.MatchedRead.Regions[i-1].End

				builder.WriteString(strconv.Itoa(numSkipped))
				builder.WriteString("I")
			}

		}
	}

	return builder.String(), nil
}

func (m ReadMatchResult) normalizeRegions() {
	// normalize regions: if a ReadMatchResult was previously incomplete and had a readRegions rv like so
	// [[50,60], [60, 70] , ..., [140,150]], after the remap happened, it looks like this
	// [[0,50], [50,60], [60, 70] , ..., [140,150]], after the remap happened, it looks like this
	// now the dimensions don't always add up -> we have to normalize the readRegions to match the genomeRegions
	// E.g if genomeRegions contains three regions (with a total length of 150), the read regions also should only
	// contain 3 regions

	adaptedReadRegions := make([]*regionvector.Region, len(m.MatchedGenome.Regions))
	originalStart := m.MatchedRead.GetFirstRegion().Start
	var startRead int
	var stopRead int
	for i, region := range m.MatchedGenome.Regions {
		startRead = regionvector.GenomicCoordToReadCoord(originalStart, region.Start, m.MatchedGenome.Regions)
		stopRead = regionvector.GenomicCoordToReadCoord(originalStart, region.End, m.MatchedGenome.Regions)
		adaptedReadRegions[i] = &regionvector.Region{Start: startRead, End: stopRead}
	}
	m.MatchedRead.Regions = adaptedReadRegions
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
