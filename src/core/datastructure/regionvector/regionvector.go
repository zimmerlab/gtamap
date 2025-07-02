package regionvector

import (
	"fmt"
	"slices"
	"sort"
	"strings"

	"github.com/sirupsen/logrus"
)

type Region struct {
	Start int // 0-based
	End   int // end-exlusive
}

type Gap struct {
	Start           int // 0-based
	End             int // end-exlusive
	KnownSpliceSite bool
}

type Intron struct {
	Start          int // 0-based
	End            int // end-exclusive
	Evidence       int
	Rank           int // rank of intron in seq
	TrueSpliceSite bool
}

func (r *Region) Length() int {
	return r.End - r.Start
}

type RegionVector struct {
	Regions []*Region
}

// GetLargestAnchor GetLargestNachor returns the largest anchor region
// WARN: Before calling this, the func MergeAlignmentBlocks needs to be called (once)
func (rv RegionVector) GetLargestAnchor(introns *RegionSet) (*Region, int, int) {
	if !rv.HasGaps() {
		return rv.Regions[0], 0, 0
	}

	maxLength := -1
	index := 0
	var largestAnchor *Region
	for i, region := range rv.Regions {
		if region.Length() > maxLength {
			maxLength = region.Length()
			largestAnchor = region
			index = i
		}
	}
	prevIntron := introns.GetPrevIntron(largestAnchor.Start)
	if prevIntron == nil {
		return largestAnchor, 0, index
	}

	return largestAnchor, prevIntron.Rank, index // anchor rank == intron rank of prev intron
}

func (r *Region) String() string {
	return fmt.Sprintf("[%d, %d]", r.Start, r.End)
}

func NewRegionVector() *RegionVector {
	return &RegionVector{
		Regions: make([]*Region, 0),
	}
}

// HasGaps checks if the region vector has any gaps between the regions.
// A gap is defined as a region where the end of one region does not equal the start of the next region.
// Only works if:
// - the regions are sorted in ascending order by their start position
// - the regions are non-overlapping
func (rv *RegionVector) HasGaps() bool {
	if len(rv.Regions) <= 1 {
		return false
	}

	for i := 0; i < len(rv.Regions)-1; i++ {
		if rv.Regions[i].End != rv.Regions[i+1].Start {
			return true
		}
	}

	return false
}

// AddRegionNonOverlappingPanic adds a region to the region vector.
// It does not allow overlapping regions and will panic if the region overlaps with any existing region.
// The regions are sorted in ascending order based on their start position.
func (rv *RegionVector) AddRegionNonOverlappingPanic(start int, end int) {
	err := rv.AddRegionNonOverlapping(start, end)
	if err != nil {
		logrus.WithFields(logrus.Fields{
			"start": start,
			"end":   end,
			"rv":    rv,
		}).Fatal("Error adding region to region vector", err)
	}
}

// AddRegionNonOverlapping adds a region to the region vector.
// It does not allow overlapping regions and will return an error if the region overlaps with any existing region.
// The regions are sorted in ascending order based on their start position.
func (rv *RegionVector) AddRegionNonOverlapping(start int, end int) error {
	return rv.AddRegionObjNonOverlapping(&Region{Start: start, End: end})
}

// AddRegionObjNonOverlapping adds a region to the region vector.
// It does not allow overlapping regions and will return an error if the region overlaps with any existing region.
// The regions are sorted in ascending order based on their start position.
func (rv *RegionVector) AddRegionObjNonOverlapping(region *Region) error {
	if region.Start == region.End {
		return fmt.Errorf("region start and end are equal")
	}
	if region.Start > region.End {
		return fmt.Errorf("region start is greater than end")
	}

	// just add the region if the vector is empty
	if len(rv.Regions) == 0 {
		rv.Regions = append(rv.Regions, region)
		return nil
	}

	// find the index where the region should be inserted based on its start position
	i := sort.Search(len(rv.Regions), func(i int) bool {
		return rv.Regions[i].Start >= region.Start
	})

	if i < len(rv.Regions) {
		if i > 0 && rv.Regions[i-1].End > region.Start {
			// the region overlaps with the previous region
			return fmt.Errorf("region overlaps with previous region")
		}
		if rv.Regions[i].Start < region.End {
			// the region is overlaps the next region
			return fmt.Errorf("region overlaps with next region")
		}
	}

	rv.Regions = slices.Insert(rv.Regions, i, region)

	return nil
}

// AddRegion adds a region to the region vector.
// The resulting region vector may contain overlapping regions and is not sorted.
// The given region will be appended to the existing regions without any checks.
func (rv *RegionVector) AddRegion(start int, end int) {
	rv.Regions = append(rv.Regions, &Region{Start: start, End: end})
}

// AddRegionAndMerge adds a region to the region vector and merges it with any overlapping regions.
// The regions are sorted in ascending order based on their start position.
func (rv *RegionVector) AddRegionAndMerge(start int, end int) {
	rv.AddRegionObjAndMerge(&Region{Start: start, End: end})
}

// AddRegionObjAndMerge adds a region to the region vector and merges it with any overlapping regions.
// The regions are sorted in ascending order based on their start position.
func (rv *RegionVector) AddRegionObjAndMerge(r *Region) {
	if len(rv.Regions) == 0 {
		rv.Regions = append(rv.Regions, r)
		return
	}

	i := sort.Search(len(rv.Regions), func(i int) bool {
		return rv.Regions[i].Start >= r.Start
	})

	rv.Regions = slices.Insert(rv.Regions, i, r)

	indexInserted := i

	if i > 0 {
		if rv.Regions[i-1].End >= r.Start {
			if rv.Regions[i-1].End < r.End {
				rv.Regions[i-1].End = r.End
			}
			rv.Regions = slices.Delete(rv.Regions, i, i+1)
			indexInserted = i - 1
		}
	}

	startIndex := indexInserted
	endIndex := indexInserted

	for i := startIndex + 1; i < len(rv.Regions); i++ {
		if rv.Regions[indexInserted].End >= rv.Regions[i].Start {
			endIndex = i
		}
	}

	if endIndex > startIndex {

		if rv.Regions[endIndex].End > rv.Regions[startIndex].End {
			rv.Regions[startIndex].End = rv.Regions[endIndex].End
		}
		rv.Regions = slices.Delete(rv.Regions, startIndex+1, endIndex+1)
	}
}

// String returns a string representation of the region vector
func (rv *RegionVector) String() string {
	result := "["
	for i, r := range rv.Regions {
		if i > 0 {
			result += ", "
		}
		result += fmt.Sprintf("[%d, %d]", r.Start, r.End)
	}
	result += "]"
	return result
}

func (rv *RegionVector) Length() int {
	length := 0
	for _, r := range rv.Regions {
		length += r.End - r.Start
	}
	return length
}

func (rv *RegionVector) GetFirstGap() *Region {
	return rv.GetGap(0)
}

// GetGapAfterRegionIndex returns the gap after the region with given index.
// It returns nil if the region index is out of bounds or if there is no gap after the region.
func (rv *RegionVector) GetGapAfterRegionIndex(regionIndex int) *Region {
	if regionIndex >= len(rv.Regions) {
		return nil
	}

	gap := Region{
		Start: rv.Regions[regionIndex].End,
		End:   rv.Regions[regionIndex+1].Start,
	}

	if gap.Start >= gap.End {
		return nil
	}

	return &gap
}

func (rv *RegionVector) GetGap(num int) *Region {
	if num < 0 || num >= len(rv.Regions)-1 {
		return nil
	}

	if len(rv.Regions) <= 1 {
		return nil
	}

	counter := 0

	for i := 0; i < len(rv.Regions)-1; i++ {
		if rv.Regions[i].End != rv.Regions[i+1].Start {
			if counter == num {
				return &Region{
					Start: rv.Regions[i].End,
					End:   rv.Regions[i+1].Start,
				}
			}
			counter++
		}
	}

	return nil
}

// GetGapIndexAfterPos returns the index of the first region after which a gap occurs,
// that starts after the given position.
// It returns -1 if there is no gap after the given position.
//
// Example:
// If the region vector contains the regions [0,10], [10,15], [20,30]
// When the position is 5, the function will return index 1 -> gap [15,20]
// When the position is 14, the function will return index 1 -> gap [15,20]
// When the position is 15, the function will return -1
// When the position is 20, the function will return -1
func (rv *RegionVector) GetGapIndexAfterPos(position int) int {
	for i := 0; i < len(rv.Regions)-1; i++ {
		if rv.Regions[i].End > position && rv.Regions[i].End < rv.Regions[i+1].Start {
			return i
		}
	}

	return -1
}

//func (rv *RegionVector) GetFirstGap() *Region {
//	if len(rv.Regions) <= 1 {
//		return nil
//	}
//
//	return &Region{
//		Start: rv.Regions[0].End,
//		End:   rv.Regions[1].Start,
//	}
//}

//func (rv *RegionVector) GetGap(gapIndex int) *Region {
//	if gapIndex >= len(rv.Regions) {
//		return nil
//	}
//
//	return &Region{
//		Start: rv.Regions[gapIndex].End,
//		End:   rv.Regions[gapIndex+1].Start,
//	}
//}

func (rv *RegionVector) GetFirstRegion() *Region {
	if len(rv.Regions) == 0 {
		return nil
	}
	return rv.Regions[0]
}

func (rv *RegionVector) GetLastRegion() *Region {
	if len(rv.Regions) == 0 {
		return nil
	}
	return rv.Regions[len(rv.Regions)-1]
}

// GetSizeLeftIncluding returns the size of all regions in the region vector that come before the
// region at the given index including the size of the region at that index.
// Returns an error if the given index is out of bounds.
func (rv *RegionVector) GetSizeLeftIncluding(regionIndex int) (int, error) {
	if regionIndex < 0 || regionIndex >= len(rv.Regions) {
		return -1, fmt.Errorf("region index out of bounds")
	}

	size := 0

	for i := 0; i <= regionIndex; i++ {
		size += rv.Regions[i].Length()
	}

	return size, nil
}

// GetSizeRightIncluding returns the size of all regions in the region vector that come after the
// region at the given index including the size of the region at that index.
// Returns an error if the given index is out of bounds.
func (rv *RegionVector) GetSizeRightIncluding(regionIndex int) (int, error) {
	if regionIndex < 0 || regionIndex >= len(rv.Regions) {
		return -1, fmt.Errorf("region index out of bounds")
	}

	size := 0

	for i := regionIndex; i < len(rv.Regions); i++ {
		size += rv.Regions[i].Length()
	}

	return size, nil
}

// GetRegionIndexContainingPosRelative returns the index of the region that contains
// the given relative position.
//
// Example:
// If the region vector contains the regions [0, 10], [20, 30], [40, 50] and the relative position
// is 5, the function will return the first region 0 [0, 10].
// If the relative position is 25, the function will return the second region 1 [20, 30].
// If the relative position is 45, the function will return the third region 2 [40, 50].
//
// If the relative position is out of bounds, the function will return an error.
func (rv *RegionVector) GetRegionIndexContainingPosRelative(relPos int) (int, error) {
	if relPos < 0 || relPos >= rv.Length() {
		return -1, fmt.Errorf("relative position out of bounds")
	}

	posDone := 0

	for i, r := range rv.Regions {
		if relPos >= posDone && relPos < posDone+r.Length() {
			return i, nil
		}
		posDone += r.Length()
	}

	return -1, fmt.Errorf("relative position not found in any region")
}

type RegionSet struct {
	Regions []*Intron
	Starts  []int
}

func (i Intron) String() string {
	var sb strings.Builder
	sb.WriteString(fmt.Sprintf("%d: [%d, %d) Confident SpliceSite: [%t] Evidence: [%d]", i.Rank, i.Start, i.End, i.TrueSpliceSite, i.Evidence))
	return sb.String()
}

func GenomicCoordToReadCoord(startInRead, genomeCoord int, genomeIntervals []*Region) int {
	pos := 0
	for _, genomicRegion := range genomeIntervals {
		if genomicRegion.End < genomeCoord {
			pos += genomicRegion.Length()
		} else {
			pos += genomeCoord - genomicRegion.Start
			break
		}
	}

	return pos + startInRead
}

func NewRegionSet(regions []*Intron) *RegionSet {
	sort.Slice(regions, func(i, j int) bool {
		return regions[i].Start < regions[j].Start
	})

	// rank introns
	pos := 0
	for _, intron := range regions {
		intron.Rank = pos
		pos++
	}

	starts := make([]int, len(regions))
	for i, r := range regions {
		starts[i] = r.Start
	}

	return &RegionSet{Regions: regions, Starts: starts}
}

// IntersectsIntrons check if readMatch is partly inside an intron
func (rs *RegionSet) IntersectsIntrons(B []*Region) bool {
	for _, b := range B {
		idx := sort.Search(len(rs.Starts), func(i int) bool {
			return rs.Starts[i] > b.End
		})

		if idx > 0 && overlaps(rs.Regions[idx-1], b) {
			return true
		}
		if idx < len(rs.Regions) && overlaps(rs.Regions[idx], b) {
			return true
		}
	}
	return false
}

func (rs *RegionSet) GetNextIntron(pos int) *Intron {
	i := 0
	for _, start := range rs.Starts {
		if start >= pos {
			return rs.Regions[i]
		}
		i++
	}
	return nil
}

func (rs *RegionSet) GetPrevIntron(pos int) *Intron {
	i := 0
	intronIndex := -1
	found := false
	for _, start := range rs.Starts {
		if start <= pos {
			found = true
			intronIndex = i
		}
		i++
	}
	if found {
		return rs.Regions[intronIndex]
	}
	return nil
}

// GetIntersectingIntron returns coords
func (rs *RegionSet) GetIntersectingIntron(b *Region) *Intron {
	idx := sort.Search(len(rs.Starts), func(i int) bool {
		return rs.Starts[i] > b.End
	})

	if idx > 0 && overlaps(rs.Regions[idx-1], b) {
		return rs.Regions[idx-1]
	}
	if idx < len(rs.Regions) && overlaps(rs.Regions[idx], b) {
		return rs.Regions[idx]
	}
	return nil
}

// GetIntersectingIntrons performs binary search to jump to first intersecting interval and
// then counts how many intervals from that point are overlapping
// stops search as soon as no overlap can be found
// returns a list of all overlapping introns
func (rs *RegionSet) GetIntersectingIntrons(b *Region) []*Intron {
	introns := make([]*Intron, 0)
	idx := sort.Search(len(rs.Starts), func(i int) bool {
		return rs.Starts[i] >= b.Start
	})

	for i := idx - 1; i >= 0 && rs.Regions[i].End > b.Start; i-- {
		if overlaps(rs.Regions[i], b) {
			introns = append(introns, rs.Regions[i])
		} else {
			break
		}
	}

	for i := idx; i < len(rs.Regions) && rs.Regions[i].Start < b.End; i++ {
		if overlaps(rs.Regions[i], b) {
			introns = append(introns, rs.Regions[i])
		} else {
			break
		}
	}
	return introns
}

// overlaps is needed to check if a region overlaps an intron and since intron is a different struct compared to region I made an extra func
func overlaps(a *Intron, b *Region) bool {
	return a.Start < b.End && b.Start < a.End
}

// OverlapsByRegion checks if the region vector overlaps with the given region.
func (rv *RegionVector) OverlapsByRegion(region *Region) bool {
	return rv.Overlaps(region.Start, region.End)
}

// Overlaps checks if the region vector overlaps with the given start and end positions.
func (rv *RegionVector) Overlaps(start int, end int) bool {
	for _, r := range rv.Regions {
		if r.Start < end && r.End > start {
			return true
		}
	}
	return false
}

// MergeAlignmentBlocks merges all blocks of the region vec and overwrites them
func (rv *RegionVector) MergeAlignmentBlocks() {
	blocks := make([]*Region, 0)
	// if there are no gaps return start and end of region vector
	if !rv.HasGaps() {
		blocks = append(blocks, &Region{
			Start: rv.GetFirstRegion().Start,
			End:   rv.GetLastRegion().End,
		})
	}

	// used to keep track of the read position for the next gap
	readGapPos := 0
	// returns the index of the first region after which a gap occurs (-1 if no gap)
	indexRegionBeforeGap := rv.GetGapIndexAfterPos(readGapPos)

	startIndex := 0

	// loop through all gaps in the read (-1 means there is no more gap)
	for indexRegionBeforeGap > -1 {
		blockStart := rv.Regions[startIndex].Start
		blockEnd := rv.Regions[indexRegionBeforeGap].End

		blocks = append(blocks, &Region{
			Start: blockStart,
			End:   blockEnd,
		})

		startIndex = indexRegionBeforeGap + 1
		readGapPos = rv.Regions[indexRegionBeforeGap].End + 1
		indexRegionBeforeGap = rv.GetGapIndexAfterPos(readGapPos)
	}

	if indexRegionBeforeGap == -1 && startIndex != 0 {
		blockStart := rv.Regions[startIndex].Start
		blockEnd := rv.GetLastRegion().End

		blocks = append(blocks, &Region{
			Start: blockStart,
			End:   blockEnd,
		})
	}
	rv.Regions = blocks
}

// GetFirstRegion returns start and ends of aligned diagonals / diagonal borders
func (rv *RegionVector) GetAlignmentBlocks() []*Region {
	blocks := make([]*Region, 0)
	// if there are no gaps return start and end of region vector
	if !rv.HasGaps() {
		blocks = append(blocks, &Region{
			Start: rv.GetFirstRegion().Start,
			End:   rv.GetLastRegion().End,
		})
	}

	// used to keep track of the read position for the next gap
	readGapPos := 0
	// returns the index of the first region after which a gap occurs (-1 if no gap)
	indexRegionBeforeGap := rv.GetGapIndexAfterPos(readGapPos)

	startIndex := 0

	// loop through all gaps in the read (-1 means there is no more gap)
	for indexRegionBeforeGap > -1 {
		blockStart := rv.Regions[startIndex].Start
		blockEnd := rv.Regions[indexRegionBeforeGap].End

		blocks = append(blocks, &Region{
			Start: blockStart,
			End:   blockEnd,
		})

		startIndex = indexRegionBeforeGap + 1
		readGapPos = rv.Regions[indexRegionBeforeGap].End + 1
		indexRegionBeforeGap = rv.GetGapIndexAfterPos(readGapPos)
	}

	if indexRegionBeforeGap == -1 && startIndex != 0 {
		blockStart := rv.Regions[startIndex].Start
		blockEnd := rv.GetLastRegion().End

		blocks = append(blocks, &Region{
			Start: blockStart,
			End:   blockEnd,
		})
	}
	return blocks
}

// Equals checks if the region vector is equal to another region vector.
// Both regions vectors are taken as is and are not sorted by this function.
// Even if the regions are the same but in different order, the function will return false.
func (rv *RegionVector) Equals(other *RegionVector) bool {
	if len(rv.Regions) == 0 && len(other.Regions) == 0 {
		return true
	}
	if len(rv.Regions) != len(other.Regions) {
		return false
	}

	for i := 0; i < len(rv.Regions); i++ {
		if rv.Regions[i].Start != other.Regions[i].Start || rv.Regions[i].End != other.Regions[i].End {
			return false
		}
	}

	return true
}

func (rv *RegionVector) Copy() *RegionVector {
	newRV := NewRegionVector()
	for _, r := range rv.Regions {
		newRV.Regions = append(newRV.Regions, &Region{
			Start: r.Start,
			End:   r.End,
		})
	}
	return newRV
}

// RemoveRegion removes a region from the region vector.
// The existing regions will be updated such that any position within the given region is not contained in
// any region anymore.
// The region vector must be sorted and non-overlapping.
// Example:
// If the region vector contains the regions [0, 10], [10, 20], [20, 30] and the given region is [5, 25]
// The resulting region vector will contain the regions [0, 5], [25, 30].
// If the region vector contains the regions [0, 30] and the given region is [5, 25],
// The resulting region vector will contain the regions [0, 5], [25, 30].
func (rv *RegionVector) _RemoveRegion(start int, end int) { // original version
	if len(rv.Regions) == 0 {
		// the region vector is empty
		return
	}

	if end <= start {
		// the region to remove is empty
		return
	}

	if end <= rv.Regions[0].Start {
		// the region to remove is before the first region
		return
	}
	if start >= rv.Regions[len(rv.Regions)-1].End {
		// the region to remove is after the last region
		return
	}
	if start <= rv.Regions[0].Start && end >= rv.Regions[len(rv.Regions)-1].End {
		// the region to remove contains all regions
		rv.Regions = make([]*Region, 0)
		return
	}

	regions := make([]*Region, 0)

	for _, r := range rv.Regions {
		if r.End <= start {
			// the region is before the region to remove
			regions = append(regions, &Region{r.Start, r.End})
			continue
		}
		if r.Start >= end {
			// the region is after the region to remove
			regions = append(regions, &Region{r.Start, r.End})
			continue
		}
		if r.Start >= start && r.End <= end {
			// the region is completely contained in the region to remove
			continue
		}
		if r.Start < start && r.End > end {
			// the region to remove is contained within the current region
			// r.Start == start && r.End == end should be handled in the previous case
			regions = append(regions, &Region{r.Start, start})
			regions = append(regions, &Region{end, r.End})
			continue
		}
		if r.Start >= start && r.End > end {
			regions = append(regions, &Region{end, r.End})
			continue
		}
		if r.Start < start && r.End <= end {
			regions = append(regions, &Region{r.Start, start})
			continue
		}
	}

	rv.Regions = regions
}

func (rv *RegionVector) UncoveredRegionsBySelf(min int, max int) *RegionVector {
	rvUncovered := NewRegionVector()

	// the total region for which the uncovered regions should be computed
	rvUncovered.AddRegion(min, max)

	// remove the covered regions from the region vector itself
	for _, r := range rv.Regions {
		rvUncovered.RemoveRegion(r.Start, r.End)
	}

	return rvUncovered
}

func (rv *RegionVector) UncoveredRegionsBySelfAndOther(rvCovered *RegionVector, min int, max int) *RegionVector {
	// the uncovered regions by the region vector itself
	rvUncovered := rv.UncoveredRegionsBySelf(min, max)

	// remove the covered regions from the given region vector
	for _, r := range rvCovered.Regions {
		rvUncovered.RemoveRegion(r.Start, r.End)
	}

	return rvUncovered
}

// RemoveRegion removes a region from the region vector.
// The existing regions will be updated such that any position within the given region is not contained in
// any region anymore.
// The region vector must be sorted and non-overlapping.
// Example:
// If the region vector contains the regions [0, 10], [10, 20], [20, 30] and the given region is [5, 25]
// The resulting region vector will contain the regions [0, 5], [25, 30].
// If the region vector contains the regions [0, 30] and the given region is [5, 25],
// The resulting region vector will contain the regions [0, 5], [25, 30].
func (rv *RegionVector) RemoveRegion(start int, end int) { // optimized
	if len(rv.Regions) == 0 {
		// the region vector is empty
		return
	}
	if end <= start {
		// the region to remove is empty
		return
	}
	if end <= rv.Regions[0].Start {
		// the region to remove is before the first region
		return
	}
	if start >= rv.Regions[len(rv.Regions)-1].End {
		// the region to remove is after the last region
		return
	}
	if start <= rv.Regions[0].Start && end >= rv.Regions[len(rv.Regions)-1].End {
		// the region to remove contains all regions
		rv.Regions = rv.Regions[:0] // reuse slice capacity
		return
	}

	// pre-allocate with worst-case capacity to avoid slice growth
	// worst case: one region splits into two, so max +1 additional region
	newRegions := make([]*Region, 0, len(rv.Regions)+1)

	for _, r := range rv.Regions {
		if r.End <= start {
			// the region is before the region to remove - reuse existing region
			newRegions = append(newRegions, r)
			continue
		}
		if r.Start >= end {
			// the region is after the region to remove - reuse existing region
			newRegions = append(newRegions, r)
			continue
		}
		if r.Start >= start && r.End <= end {
			// the region is completely contained in the region to remove
			continue
		}
		if r.Start < start && r.End > end {
			// the region to remove is contained within the current region
			// Split into two regions - only allocate new ones when needed
			newRegions = append(newRegions, &Region{r.Start, start})
			newRegions = append(newRegions, &Region{end, r.End})
			continue
		}
		if r.Start >= start && r.End > end {
			// Modify existing region in-place if possible, or create new one
			if r.Start == end {
				// Can reuse the existing region
				newRegions = append(newRegions, r)
			} else {
				newRegions = append(newRegions, &Region{end, r.End})
			}
			continue
		}
		if r.Start < start && r.End <= end {
			// Modify existing region in-place if possible, or create new one
			if r.End == start {
				// Can reuse the existing region
				newRegions = append(newRegions, r)
			} else {
				newRegions = append(newRegions, &Region{r.Start, start})
			}
			continue
		}
	}
	rv.Regions = newRegions
}
