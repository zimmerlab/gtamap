package mapperutils

import (
	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
)

// ComputeGapsInDiagonal TODO: write comment
// Gaps are dis-continuous regions in the diagonal which do not contain any region already mapped.
func ComputeGapsInDiagonal(diagonalRead *regionvector.RegionVector, diagonalGenome *regionvector.RegionVector,
	result *ReadMatchResult,
) (*regionvector.RegionVector, *regionvector.RegionVector) {
	gapsRead := regionvector.NewRegionVector()
	gapsGenome := regionvector.NewRegionVector()

	for i := 0; i < len(diagonalRead.Regions)-1; i++ {

		// no gap, skip region
		if diagonalRead.Regions[i].End == diagonalRead.Regions[i+1].Start {
			continue
		}
		// the gap contains regions that are already mapped
		if result.MatchedRead.Overlaps(diagonalRead.Regions[i].End, diagonalRead.Regions[i+1].Start) {
			continue
		}
		// the gap does not include any region that was already mapped and can be filled
		gapsRead.AddRegionNonOverlappingPanic(diagonalRead.Regions[i].End, diagonalRead.Regions[i+1].Start)
		gapsGenome.AddRegionNonOverlappingPanic(diagonalGenome.Regions[i].End, diagonalGenome.Regions[i+1].Start)
	}

	return gapsRead, gapsGenome
}
