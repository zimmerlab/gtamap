package mapperutils

import (
	"fmt"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/datastructure"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/interval"
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
	SequenceIndex            int                        // the index of the sequence in the genome
	MatchedRead              *regionvector.RegionVector // region vector containing the matched positions in the read
	MatchedGenome            *regionvector.RegionVector // region vector containing the matched positions in the genome
	MismatchesRead           []int                      // the positions of the mismatches in the read
	MismatchCounts           map[string]int
	MismatchConstraintGlobal int
	diagonalHandler          *DiagonalHandler
	IncompleteMap            bool
	SpliceSitesInfo          []int // corresponds to the number of junctions of match result. canonical splice site -> 2 = canonical, 1 = non-can, 0 = no splice site

	Id int // fileds for logging

	IsFixPoint          bool
	MainAnchorLength    int
	MainAnchorMM        int
	TotalLeftOptions    int
	TotalRightOptions   int
	ValidLeftOptions    int
	ValidRightOptions   int
	LeftFixpointLength  int
	RightFixpointLength int

	IsOverhangCorrected bool

	IsGapFill          bool
	IsGapFillOverflow  bool
	GapsFilled         []int
	GapsFilledOverflow []int

	IsSymInErr  bool
	SymInErrLen []int
}

func (r *ReadMatchResult) String() string {
	sb := strings.Builder{}

	sb.Write(fmt.Appendf(nil, "ReadMatchResult on seq index %d (incomplete: %t)\n", r.SequenceIndex, r.IncompleteMap))

	maxLenReadStart := 0
	maxLenReadEnd := 0
	maxLenGenomeStart := 0
	maxLenGenomeEnd := 0

	for i := 0; i < len(r.MatchedRead.Regions); i++ {

		read := r.MatchedRead.Regions[i]
		genome := r.MatchedGenome.Regions[i]

		maxLenReadStart = max(maxLenReadStart, len(strconv.Itoa(read.Start)))
		maxLenReadEnd = max(maxLenReadEnd, len(strconv.Itoa(read.End)))
		maxLenGenomeStart = max(maxLenGenomeStart, len(strconv.Itoa(genome.Start)))
		maxLenGenomeEnd = max(maxLenGenomeEnd, len(strconv.Itoa(genome.End)))
	}

	for i := 0; i < len(r.MatchedRead.Regions); i++ {

		read := r.MatchedRead.Regions[i]
		genome := r.MatchedGenome.Regions[i]

		if i > 0 && read.Start > r.MatchedRead.Regions[i-1].End {
			sb.Write(fmt.Appendf(nil, "%3d: %*d - %*d <-> \n", i, maxLenReadStart, r.MatchedRead.Regions[i-1].End, maxLenReadEnd, read.Start))
		}

		sb.Write(fmt.Appendf(nil, "%3d: %*d - %*d <-> %*d - %*d\n", i, maxLenReadStart, read.Start, maxLenReadEnd, read.End, maxLenGenomeStart, genome.Start, maxLenGenomeEnd, genome.End))
	}

	sb.Write(fmt.Appendf(nil, "Mismatches in read: %v\n", r.MismatchesRead))

	return sb.String()
}

func (r *ReadMatchResult) HasUnknownSpliceSites() bool {
	if r.SpliceSitesInfo == nil {
		return false
	}
	for _, site := range r.SpliceSitesInfo {
		if site > 0 {
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

func (r *ReadMatchResult) MergeRegions() {
	mergedRead := []regionvector.Region{r.MatchedRead.Regions[0]}
	mergedGenome := []regionvector.Region{r.MatchedGenome.Regions[0]}

	for i := 1; i < len(r.MatchedRead.Regions); i++ {
		last := &mergedRead[len(mergedRead)-1]
		current := r.MatchedRead.Regions[i]

		if current.Start <= last.End { // Overlapping or adjacent
			if current.End > last.End {
				last.End = current.End
			}
		} else {
			mergedRead = append(mergedRead, current)
		}
	}

	for i := 1; i < len(r.MatchedGenome.Regions); i++ {
		last := &mergedGenome[len(mergedGenome)-1]
		current := r.MatchedGenome.Regions[i]

		if current.Start <= last.End { // Overlapping or adjacent
			if current.End > last.End {
				last.End = current.End
			}
		} else {
			mergedGenome = append(mergedGenome, current)
		}
	}
	r.MatchedGenome.Regions = mergedGenome
	r.MatchedRead.Regions = mergedRead
}

// MergeRegions merges the read and genome regions if both are consecutive.
// This preserves the property that the number of read and genome regions are
// the same and that each read region is associated with the genome region at
// the same index.
func (r *ReadMatchResult) MergeRegionsIfBothConsecutive() {

	if len(r.MatchedRead.Regions) != len(r.MatchedGenome.Regions) {
		logrus.WithFields(logrus.Fields{
			"len read regions":   len(r.MatchedRead.Regions),
			"len genome regions": len(r.MatchedGenome.Regions),
		}).Error("Cant merge regions of ReadMatchResult: different number of read and genome regions")
		return
	}

	mergedRead := make([]regionvector.Region, 1)
	mergedRead[0] = regionvector.Region{
		Start: r.MatchedRead.Regions[0].Start,
		End:   r.MatchedRead.Regions[0].End,
	}

	mergedGenome := make([]regionvector.Region, 1)
	mergedGenome[0] = regionvector.Region{
		Start: r.MatchedGenome.Regions[0].Start,
		End:   r.MatchedGenome.Regions[0].End,
	}

	iMerged := 0

	for i := 1; i < len(r.MatchedRead.Regions); i++ {

		isConsecutiveRead := r.MatchedRead.Regions[i].Start == mergedRead[iMerged].End
		isConsecutiveGenome := r.MatchedGenome.Regions[i].Start == mergedGenome[iMerged].End

		if isConsecutiveRead && isConsecutiveGenome {
			// extend the last region
			mergedRead[iMerged].End = r.MatchedRead.Regions[i].End
			mergedGenome[iMerged].End = r.MatchedGenome.Regions[i].End
		} else {
			// add a new region
			mergedRead = append(mergedRead, regionvector.Region{
				Start: r.MatchedRead.Regions[i].Start,
				End:   r.MatchedRead.Regions[i].End,
			})
			mergedGenome = append(mergedGenome, regionvector.Region{
				Start: r.MatchedGenome.Regions[i].Start,
				End:   r.MatchedGenome.Regions[i].End,
			})
			iMerged++
		}
	}

	r.MatchedRead.Regions = mergedRead
	r.MatchedGenome.Regions = mergedGenome
}

func (r *ReadMatchResult) IsValid() bool {

	numGapsGenome := 0

	for i := 0; i < len(r.MatchedGenome.Regions)-1; i++ {

		lenGapGenome := r.MatchedGenome.Regions[i+1].Start - r.MatchedGenome.Regions[i].End

		if lenGapGenome == 0 {
			continue
		}

		numGapsGenome++

		if config.Mapper.Mapping.IsReadOriginRna {

		} else {
			if lenGapGenome > config.Mapper.Mapping.DnaMode.MaxGapLength {
				return false
			}
			if numGapsGenome > config.Mapper.Mapping.DnaMode.MaxGapCount {
				return false
			}
		}
	}

	return true
}

func (r *ReadMatchResult) NormalizeRegions() {

	r.MergeRegions()
	readRegions := r.MatchedRead.Regions
	genomeRegions := r.MatchedGenome.Regions

	var normalizedRead []regionvector.Region
	var normalizedGenome []regionvector.Region

	readIndex, genomeIndex := 0, 0
	readOffset, genomeOffset := 0, 0

	for readIndex < len(readRegions) && genomeIndex < len(genomeRegions) {
		readRegion := readRegions[readIndex]
		genomeRegion := genomeRegions[genomeIndex]

		readLen := readRegion.End - readRegion.Start - readOffset
		genomeLen := genomeRegion.End - genomeRegion.Start - genomeOffset

		blockLen := min(readLen, genomeLen)

		newRead := regionvector.Region{
			Start: readRegion.Start + readOffset,
			End:   readRegion.Start + readOffset + blockLen,
		}
		newGenome := regionvector.Region{
			Start: genomeRegion.Start + genomeOffset,
			End:   genomeRegion.Start + genomeOffset + blockLen,
		}

		normalizedRead = append(normalizedRead, newRead)
		normalizedGenome = append(normalizedGenome, newGenome)

		readOffset += blockLen
		genomeOffset += blockLen

		if readOffset >= readRegion.Length() {
			readIndex++
			readOffset = 0
		}
		if genomeOffset >= genomeRegion.Length() {
			genomeIndex++
			genomeOffset = 0
		}
	}

	r.MatchedRead.Regions = normalizedRead
	r.MatchedGenome.Regions = normalizedGenome
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// func (r *ReadMatchResult) _NormalizeRegions() { // old version DID NOT WORK
// 	// I. In some cases the regions differ in length but in different positions in genome and read:
// 	// READ [10,30] [30,40]
// 	// GENO [10,20] [20,40] -> split into regions of kmer length before normalizing
// 	kmerLength := int(config.KmerLength())
// 	for _, region := range r.MatchedGenome.Regions {
// 		if region.Length() != kmerLength && region.Length()%kmerLength == 0 {
// 			start := region.Start
// 			// remove old region
// 			r.MatchedGenome.RemoveRegion(region.Start, region.End)
// 			for i := 0; i < region.Length(); i += kmerLength {
// 				r.MatchedGenome.AddRegionNonOverlappingPanic(start+i, start+i+kmerLength)
// 			}
// 		}
// 	}
// 	for _, region := range r.MatchedRead.Regions {
// 		if region.Length() != kmerLength && region.Length()%kmerLength == 0 {
// 			start := region.Start
// 			// remove old region
// 			r.MatchedRead.RemoveRegion(region.Start, region.End)
// 			for i := 0; i < region.Length(); i += kmerLength {
// 				r.MatchedRead.AddRegionNonOverlappingPanic(start+i, start+i+kmerLength)
// 			}
// 		}
// 	}
//
// 	genomeBlocks := make([]regionvector.Region, 0)
// 	readBlocks := make([]regionvector.Region, 0)
//
// 	startIndex := 0
//
// 	hasGap := func(a, b regionvector.Region) bool {
// 		return a.End < b.Start
// 	}
//
// 	for i := 0; i < len(r.MatchedGenome.Regions)-1; i++ {
// 		gapInGenome := hasGap(r.MatchedGenome.Regions[i], r.MatchedGenome.Regions[i+1])
// 		gapInRead := hasGap(r.MatchedRead.Regions[i], r.MatchedRead.Regions[i+1])
//
// 		if gapInGenome || gapInRead {
// 			// account for this case
// 			// Read [0,10] [10,20] [20,30] [30,40]
// 			// GEN  [0,10] [10,20] [20,30] [30,36] -> due to either leftNorm or bcause of best split
// 			//if r.MatchedGenome.Regions[i].Length() < r.MatchedRead.Regions[i].Length() {
// 			//	padding := r.MatchedRead.Regions[i].Length() - r.MatchedGenome.Regions[i].Length()
// 			//	r.MatchedRead.Regions[i].End -= padding
// 			//	r.MatchedRead.Regions[i+1].Start -= padding
// 			//} else if r.MatchedGenome.Regions[i].Length() > r.MatchedRead.Regions[i].Length() {
// 			//	padding := r.MatchedGenome.Regions[i].Length() - r.MatchedRead.Regions[i].Length()
// 			//	r.MatchedRead.Regions[i].End += padding
// 			//	r.MatchedRead.Regions[i+1].Start += padding
// 			//}
// 			genomeBlocks = append(genomeBlocks, regionvector.Region{
// 				Start: r.MatchedGenome.Regions[startIndex].Start,
// 				End:   r.MatchedGenome.Regions[i].End,
// 			})
// 			readBlocks = append(readBlocks, regionvector.Region{
// 				Start: r.MatchedRead.Regions[startIndex].Start,
// 				End:   r.MatchedRead.Regions[i].End,
// 			})
// 			startIndex = i + 1
// 		}
// 	}
//
// 	if startIndex < len(r.MatchedGenome.Regions) {
// 		genomeBlocks = append(genomeBlocks, regionvector.Region{
// 			Start: r.MatchedGenome.Regions[startIndex].Start,
// 			End:   r.MatchedGenome.Regions[len(r.MatchedGenome.Regions)-1].End,
// 		})
// 		readBlocks = append(readBlocks, regionvector.Region{
// 			Start: r.MatchedRead.Regions[startIndex].Start,
// 			End:   r.MatchedRead.Regions[len(r.MatchedRead.Regions)-1].End,
// 		})
// 	}
//
// 	r.MatchedGenome.Regions = genomeBlocks
// 	r.MatchedRead.Regions = readBlocks
// }

func (r *ReadMatchResult) Copy() *ReadMatchResult {
	var dhCopy *DiagonalHandler
	if r.diagonalHandler != nil {
		dhCopy = r.diagonalHandler.Copy()
	}

	mismatchCounts := make(map[string]int)
	for key, val := range r.MismatchCounts {
		mismatchCounts[key] = val
	}

	cSymInErrLen := make([]int, 0)
	if r.SymInErrLen != nil {
		cSymInErrLen := make([]int, len(r.SymInErrLen))
		copy(cSymInErrLen, r.SymInErrLen)
	}

	cGapsFilled := make([]int, 0)
	if r.GapsFilled != nil {
		cGapsFilled = make([]int, len(r.GapsFilled))
		copy(cGapsFilled, r.GapsFilled)
	}
	cGapsFilledOverflow := make([]int, 0)
	if r.GapsFilledOverflow != nil {
		cGapsFilledOverflow = make([]int, len(r.GapsFilledOverflow))
		copy(cGapsFilledOverflow, r.GapsFilledOverflow)
	}

	return &ReadMatchResult{
		SequenceIndex:            r.SequenceIndex,
		MatchedRead:              r.MatchedRead.Copy(),
		MatchedGenome:            r.MatchedGenome.Copy(),
		MismatchesRead:           append([]int{}, r.MismatchesRead...),
		MismatchCounts:           mismatchCounts,
		MismatchConstraintGlobal: r.MismatchConstraintGlobal,
		diagonalHandler:          dhCopy,

		IsFixPoint:          r.IsFixPoint,
		MainAnchorLength:    r.MainAnchorLength,
		MainAnchorMM:        r.MainAnchorMM,
		TotalLeftOptions:    r.TotalLeftOptions,
		TotalRightOptions:   r.TotalRightOptions,
		ValidLeftOptions:    r.ValidLeftOptions,
		ValidRightOptions:   r.ValidRightOptions,
		LeftFixpointLength:  r.LeftFixpointLength,
		RightFixpointLength: r.RightFixpointLength,

		IsOverhangCorrected: r.IsOverhangCorrected,

		IsGapFill:          r.IsGapFill,
		GapsFilled:         cGapsFilled,
		IsGapFillOverflow:  r.IsGapFillOverflow,
		GapsFilledOverflow: cGapsFilledOverflow,

		IsSymInErr:  r.IsSymInErr,
		SymInErrLen: cSymInErrLen,
	}
}

func AssignReadMatchResults(
	fwMappings []*ReadMatchResult,
	rvMappings []*ReadMatchResult,
) (
	map[int][]*ReadMatchResult,
	map[int][]*ReadMatchResult,
	map[int]struct{},
) {

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
	// INFO: DNA RNA MODE
	// For DNA, l and r regions of junction should be longer
	for _, region := range mapping.MatchedGenome.Regions {
		if config.IsOriginRNA {
			if region.Length() < config.MinConfAnchorLengthRNA {
				return false
			}
		} else {
			if region.Length() < config.MinConfAnchorLengthDNA {
				return false
			}
		}
	}
	return true
}

// AddMismatch tries to add a mismatch at the given position in the read.
// It checks whether adding the mismatch would exceed the global or region-specific
// mismatch constraints. If adding the mismatch is valid, it updates the mismatch
// counts and returns true. If adding the mismatch would violate any constraints,
// it returns false and does not modify the state.
// posInRead is the position of the mismatch in the read (0-based).
// posInTargetRegion is the position of the mismatch in the target region (0-based).
// It throws a fatal error if a duplicate mismatch position is added.
func (r *ReadMatchResult) AddMismatch(
	regionMask *interval.PriorityList,
	posInRead int,
	posInTargetRegion int,
) bool {

	for _, mm := range r.MismatchesRead {
		if mm == posInRead {
			logrus.Fatal("Trying to add duplicate mismatch position to "+
				"ReadMatchResult:", r.MismatchesRead, "posInRead:", posInRead)
		}
	}

	// mismatches exceed global mismatch constraint
	if len(r.MismatchesRead)+1 > r.MismatchConstraintGlobal {
		return false
	}

	found, name, _, maxMismatches := regionMask.GetItemAtPosition(posInTargetRegion)
	if !found {
		name = "unmasked"
	}

	if _, foundName := r.MismatchCounts[name]; !foundName {
		r.MismatchCounts[name] = 0
	}

	// mismatches exceed region specific mismatch constraint
	if found && r.MismatchCounts[name]+1 >= maxMismatches {
		return false
	}

	// add the mismatch count to the region specific mismatch counts
	r.MismatchCounts[name] += 1
	// add the read position of the mismatch to the result
	r.MismatchesRead = append(r.MismatchesRead, posInRead)

	// the read is still valid even after adding the mismatch
	return true
}

func (r *ReadMatchResult) GetCigar() (string, error) {

	var builder strings.Builder

	// slices.Sort(r.MismatchesRead)
	sort.Ints(r.MismatchesRead)

	// for i := 1; i < len(r.MismatchesRead); i++ {
	// 	if r.MismatchesRead[i] == r.MismatchesRead[i-1] {
	// 		logrus.Fatal("Duplicate mismatch positions in read:", r.MismatchesRead)
	// 	}
	// }

	r.NormalizeRegions()

	// TODO: remove this sanity check
	if len(r.MatchedRead.Regions) != len(r.MatchedGenome.Regions) {
		logrus.Fatal("len(r.MatchedRead.Regions) != len(r.MatchedGenome.Regions) even after NormalizeRegions()")
	}
	if r.MatchedRead.Length() != r.MatchedGenome.Length() {
		logrus.Fatal("r.MatchedRead.Length() != r.MatchedGenome.Length() even after NormalizeRegions()")
	}

	isForwardStrand := r.SequenceIndex == 0

	if isForwardStrand {

		// number of cumulative aligned positions (matches or mismatches) between the read and the genome
		numMatchesSum := 0

		// TODO: is this needed when NormalizeRegions is called above?
		if len(r.MatchedGenome.Regions) != len(r.MatchedRead.Regions) {
			// r.SyncRegions()
			r.NormalizeRegions()
		}

		// the regions in MatchedGenome and MatchedRead have the same dimensions
		for i := 0; i < len(r.MatchedGenome.Regions); i++ {

			// each pair of regions represent a match
			numMatchesSum += r.MatchedGenome.Regions[i].Length()

			// if this is the last match, add the number of matches
			if i == len(r.MatchedGenome.Regions)-1 {
				if config.IncludeMMinSAM {
					builder.WriteString(ParseMatchedRegion(r.MatchedRead.Regions[i], r.MismatchesRead, false))
				} else {
					builder.WriteString(strconv.Itoa(numMatchesSum))
					builder.WriteString("M")
				}
				break
			}

			gapInGenome := r.MatchedGenome.Regions[i].End < r.MatchedGenome.Regions[i+1].Start
			gapInRead := r.MatchedRead.Regions[i].End < r.MatchedRead.Regions[i+1].Start

			if gapInGenome && gapInRead {
				logrus.WithFields(logrus.Fields{
					"read":   r.MatchedRead,
					"genome": r.MatchedGenome,
				}).Warn("Gap in genome and read at the same time (fw)")

				return "", fmt.Errorf("gap in genome and read at the same time (fw)")
			}

			if !gapInGenome && !gapInRead {
				continue
			}

			if config.IncludeMMinSAM {
				builder.WriteString(ParseMatchedRegion(r.MatchedRead.Regions[i], r.MismatchesRead, false))
			} else {
				builder.WriteString(strconv.Itoa(numMatchesSum))
				builder.WriteString("M")
			}

			numMatchesSum = 0

			// intron or deletion
			if gapInGenome {
				// number of skipped bases in the reference
				numSkipped := r.MatchedGenome.Regions[i+1].Start - r.MatchedGenome.Regions[i].End

				builder.WriteString(strconv.Itoa(numSkipped))

				if !config.IsOriginRNA || numSkipped < config.IntronLengthMin() {
					// there are no introns in DNA mode
					// configured threshold is used to determine whether a gap is a deletion or an intron
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
			// r.SyncRegions()
			r.NormalizeRegions()
		}

		// the regions in MatchedGenome and MatchedRead have the same dimensions
		for i := len(r.MatchedGenome.Regions) - 1; i >= 0; i-- {

			numMatchesSum += r.MatchedGenome.Regions[i].Length()

			// if this is the last match, add the number of matches
			if i == 0 {
				if config.IncludeMMinSAM {
					builder.WriteString(ParseMatchedRegion(r.MatchedRead.Regions[i], r.MismatchesRead, true))
				} else {
					builder.WriteString(strconv.Itoa(numMatchesSum))
					builder.WriteString("M")
				}
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

			if config.IncludeMMinSAM {
				builder.WriteString(ParseMatchedRegion(r.MatchedRead.Regions[i], r.MismatchesRead, true))
			} else {
				builder.WriteString(strconv.Itoa(numMatchesSum))
				builder.WriteString("M")
			}

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
// func (r *ReadMatchResult) SyncRegions() {
// 	adaptedReadRegions := make([]regionvector.Region, len(r.MatchedGenome.Regions))
// 	originalStart, ok := r.MatchedRead.GetFirstRegion()
// 	if ok {
// 		for i, region := range r.MatchedGenome.Regions {
// 			startRead, err := regionvector.GenomicCoordToReadCoord(originalStart.Start, region.Start, r.MatchedGenome.Regions)
// 			if err != nil {
// 				logrus.Errorf("Error while converting genomic coord to read coord")
// 				logrus.Fatal(err)
// 			}
// 			stopRead, err := regionvector.GenomicCoordToReadCoord(originalStart.Start, region.End, r.MatchedGenome.Regions)
// 			if err != nil {
// 				logrus.Errorf("Error while converting genomic coord to read coord")
// 				logrus.Fatal(err)
// 			}
// 			adaptedReadRegions[i] = regionvector.Region{Start: startRead, End: stopRead}
// 		}
// 		r.MatchedRead.Regions = adaptedReadRegions
// 	}
// }

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
	IntronTrees    map[int]*datastructure.INTree   // interval trees for intronsmaps
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
		"pct of reads pairs with preffered orientation": t.Confidence * 100,
	}).Debug("Done with Annotation")

	for orientation, introns := range t.Introns {
		orientationLabel := "[+]"
		if orientation == 1 {
			orientationLabel = "[-]"
		}
		logrus.WithFields(logrus.Fields{
			"orientation": orientationLabel,
		}).Debug("Inferred Introns of")

		for _, intron := range introns.Regions {
			logrus.WithFields(logrus.Fields{
				"Intron": intron.String(),
			}).Debug()
		}
	}
}

func ParseMatchedRegion(
	region regionvector.Region,
	mismatchesInRead []int,
	isRev bool,
) string {

	var builder strings.Builder

	if len(mismatchesInRead) == 0 {
		return fmt.Sprintf("%d=", region.Length())
	}

	lastStart := 0      // relative to region start
	var lastMM int = -2 // to track consecutive mismatches
	count := 0          // count of consecutive mismatches

	flushMismatch := func() {
		if count > 0 {
			builder.WriteString(fmt.Sprintf("%dX", count))
			count = 0
		}
	}

	if isRev {
		for i := len(mismatchesInRead) - 1; i >= 0; i-- {
			mm := mismatchesInRead[i]
			if mm < region.Start || mm >= region.End {
				continue
			}
			relativePos := region.End - mm - 1

			if relativePos > lastStart {
				flushMismatch()
				builder.WriteString(fmt.Sprintf("%d=", relativePos-lastStart))
			}

			if lastMM == relativePos-1 {
				count++
			} else {
				flushMismatch()
				count = 1
			}
			lastMM = relativePos
			lastStart = relativePos + 1
		}
	} else {

		for _, mm := range mismatchesInRead {
			if mm < region.Start || mm >= region.End {
				continue
			}
			relativePos := mm - region.Start

			if relativePos > lastStart {
				flushMismatch()
				builder.WriteString(fmt.Sprintf("%d=", relativePos-lastStart))
			}

			if lastMM == relativePos-1 {
				count++
			} else {
				flushMismatch()
				count = 1
			}
			lastMM = relativePos
			lastStart = relativePos + 1
		}
	}

	flushMismatch()

	if lastStart < region.Length() {
		builder.WriteString(fmt.Sprintf("%d=", region.Length()-lastStart))
	}

	return builder.String()
}

func sliceToString(s []int) string {
	if len(s) == 0 {
		return "NA"
	}
	strs := make([]string, len(s))
	for i, v := range s {
		strs[i] = strconv.Itoa(v)
	}
	return strings.Join(strs, ",")
}

func (r *ReadMatchResult) WriteTSV(f *os.File, gene string, isFw int, altId int, readId string) error {
	fields := []string{
		gene,
		readId,
		strconv.Itoa(r.SequenceIndex),
		strconv.Itoa(isFw),
		strconv.Itoa(altId),
		strconv.Itoa(len(r.MismatchesRead)),
		strconv.FormatBool(r.IsFixPoint),
		strconv.Itoa(r.MainAnchorLength),
		strconv.Itoa(r.MainAnchorMM),
		strconv.Itoa(r.TotalLeftOptions),
		strconv.Itoa(r.TotalRightOptions),
		strconv.Itoa(r.ValidLeftOptions),
		strconv.Itoa(r.ValidRightOptions),
		strconv.Itoa(r.LeftFixpointLength),
		strconv.Itoa(r.RightFixpointLength),
		strconv.FormatBool(r.IsOverhangCorrected),
		strconv.FormatBool(r.IsGapFill),
		strconv.FormatBool(r.IsGapFillOverflow),
		sliceToString(r.GapsFilled),
		sliceToString(r.GapsFilledOverflow),
		strconv.FormatBool(r.IsSymInErr),
		sliceToString(r.SymInErrLen),
	}
	line := strings.Join(fields, "\t") + "\n"
	_, err := f.WriteString(line)
	return err
}
