package regionvector

import (
	"fmt"
	"slices"
	"sort"
)

type Region struct {
	Start int // 0-based
	End   int // end-exlusive
}

func (r *Region) Length() int {
	return r.End - r.Start
}

type RegionVector struct {
	Regions []*Region
}

func NewRegionVector() *RegionVector {
	return &RegionVector{
		Regions: make([]*Region, 0),
	}
}

func (rv *RegionVector) AddRegion(start int, end int) {
	rv.AddRegionObj(&Region{Start: start, End: end})
}

func (rv *RegionVector) AddRegionObj(r *Region) {

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
	if len(rv.Regions) <= 1 {
		return nil
	}

	return &Region{
		Start: rv.Regions[0].End,
		End:   rv.Regions[1].Start,
	}
}

func (rv *RegionVector) GetGap(gapIndex int) *Region {
	if gapIndex >= len(rv.Regions) {
		return nil
	}

	return &Region{
		Start: rv.Regions[gapIndex].End,
		End:   rv.Regions[gapIndex+1].Start,
	}
}

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
