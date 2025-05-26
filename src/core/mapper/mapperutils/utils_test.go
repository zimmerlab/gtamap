package mapperutils

import (
	"testing"

	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
)

func TestComputeGapsInDiagonal(t *testing.T) {
	resultMatchedRead := regionvector.NewRegionVector()
	resultMatchedRead.AddRegionNonOverlappingPanic(100, 110)
	resultMatchedRead.AddRegionNonOverlappingPanic(110, 120)
	resultMatchedRead.AddRegionNonOverlappingPanic(130, 140)
	resultMatchedRead.AddRegionNonOverlappingPanic(190, 200)

	resultMatchedGenome := regionvector.NewRegionVector()
	resultMatchedGenome.AddRegionNonOverlappingPanic(1000, 1010)
	resultMatchedGenome.AddRegionNonOverlappingPanic(1010, 1020)
	resultMatchedGenome.AddRegionNonOverlappingPanic(1030, 1040)
	resultMatchedGenome.AddRegionNonOverlappingPanic(1090, 1100)

	result := ReadMatchResult{
		SequenceIndex:  0,
		MatchedRead:    resultMatchedRead,
		MatchedGenome:  resultMatchedGenome,
		MismatchesRead: nil,
		NeedRemap:      false,
	}

	// testcase 1 where the gap in the diagonal is not overlapping with the matched regions
	diagonalReadCase1 := regionvector.NewRegionVector()
	diagonalReadCase1.AddRegionNonOverlappingPanic(0, 10)
	diagonalReadCase1.AddRegionNonOverlappingPanic(40, 50)
	diagonalGenomeCase1 := regionvector.NewRegionVector()
	diagonalGenomeCase1.AddRegionNonOverlappingPanic(900, 910)
	diagonalGenomeCase1.AddRegionNonOverlappingPanic(940, 950)
	gapsReadCase1 := regionvector.NewRegionVector()
	gapsReadCase1.AddRegionNonOverlappingPanic(10, 40)
	gapsGenomeCase1 := regionvector.NewRegionVector()
	gapsGenomeCase1.AddRegionNonOverlappingPanic(910, 940)

	// no gap in the diagonal should return nil
	diagonalReadCase2 := regionvector.NewRegionVector()
	diagonalReadCase2.AddRegionNonOverlappingPanic(0, 10)
	diagonalGenomeCase2 := regionvector.NewRegionVector()
	diagonalGenomeCase2.AddRegionNonOverlappingPanic(900, 910)

	// multiple gaps in the diagonal of which one overlaps a region already mapped
	diagonalReadCase3 := regionvector.NewRegionVector()
	diagonalReadCase3.AddRegionNonOverlappingPanic(0, 10)
	diagonalReadCase3.AddRegionNonOverlappingPanic(40, 50)
	diagonalReadCase3.AddRegionNonOverlappingPanic(80, 90)
	diagonalReadCase3.AddRegionNonOverlappingPanic(120, 130)
	diagonalGenomeCase3 := regionvector.NewRegionVector()
	diagonalGenomeCase3.AddRegionNonOverlappingPanic(900, 910)
	diagonalGenomeCase3.AddRegionNonOverlappingPanic(940, 950)
	diagonalGenomeCase3.AddRegionNonOverlappingPanic(980, 990)
	diagonalGenomeCase3.AddRegionNonOverlappingPanic(1020, 1030)
	gapsReadCase3 := regionvector.NewRegionVector()
	gapsReadCase3.AddRegionNonOverlappingPanic(10, 40)
	gapsReadCase3.AddRegionNonOverlappingPanic(50, 80)
	gapsGenomeCase3 := regionvector.NewRegionVector()
	gapsGenomeCase3.AddRegionNonOverlappingPanic(910, 940)
	gapsGenomeCase3.AddRegionNonOverlappingPanic(950, 980)

	testCases := []struct {
		name               string
		diagonalRead       *regionvector.RegionVector
		diagonalGenome     *regionvector.RegionVector
		expectedGapsRead   *regionvector.RegionVector
		expectedGapsGenome *regionvector.RegionVector
	}{
		{
			name:               "Gap in diagonal not overlapping with matched regions",
			diagonalRead:       diagonalReadCase1,
			diagonalGenome:     diagonalGenomeCase1,
			expectedGapsRead:   gapsReadCase1,
			expectedGapsGenome: gapsGenomeCase1,
		},
		{
			name:               "No gap in diagonal",
			diagonalRead:       diagonalReadCase2,
			diagonalGenome:     diagonalGenomeCase2,
			expectedGapsRead:   regionvector.NewRegionVector(),
			expectedGapsGenome: regionvector.NewRegionVector(),
		},
		{
			name:               "Multiple gaps in diagonal with one overlapping with matched regions",
			diagonalRead:       diagonalReadCase3,
			diagonalGenome:     diagonalGenomeCase3,
			expectedGapsRead:   gapsReadCase3,
			expectedGapsGenome: gapsGenomeCase3,
		},
	}

	for _, testCase := range testCases {
		t.Run(testCase.name, func(t *testing.T) {
			gotGapsRead, gotGapsGenome := ComputeGapsInDiagonal(testCase.diagonalRead, testCase.diagonalGenome, &result)

			if testCase.expectedGapsRead == nil {
				if gotGapsRead != nil {
					t.Errorf("expected gaps read nil, got %v", gotGapsRead)
				}
				if gotGapsGenome != nil {
					t.Errorf("expected gaps genome nil, got %v", gotGapsGenome)
				}
				return
			}

			if !testCase.expectedGapsRead.Equals(gotGapsRead) {
				t.Errorf("expected gaps read %v, got %v", testCase.expectedGapsRead, gotGapsRead)
			}
			if !testCase.expectedGapsGenome.Equals(gotGapsGenome) {
				t.Errorf("expected gaps genome %v, got %v", testCase.expectedGapsGenome, gotGapsGenome)
			}
		})
	}
}
