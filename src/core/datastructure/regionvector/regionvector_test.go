package regionvector

import "testing"

func TestRemoveRegion(t *testing.T) {

	// the region vector does not contain any regions
	case1InputRv := NewRegionVector()
	case1InputRegion := &Region{0, 10}
	case1ExpectedRv := NewRegionVector()

	// the deleted region is empty
	case2InputRv := NewRegionVector()
	case2InputRv.AddRegionNonOverlappingPanic(0, 10)
	case2InputRegion := &Region{0, 0}
	case2ExpectedRv := NewRegionVector()
	case2ExpectedRv.AddRegionNonOverlappingPanic(0, 10)

	// the deleted region is before the first region
	case3InputRv := NewRegionVector()
	case3InputRv.AddRegionNonOverlappingPanic(10, 15)
	case3InputRegion := &Region{0, 5}
	case3ExpectedRv := NewRegionVector()
	case3ExpectedRv.AddRegionNonOverlappingPanic(10, 15)

	// the deleted region is after the last region
	case4InputRv := NewRegionVector()
	case4InputRv.AddRegionNonOverlappingPanic(10, 15)
	case4InputRegion := &Region{20, 25}
	case4ExpectedRv := NewRegionVector()
	case4ExpectedRv.AddRegionNonOverlappingPanic(10, 15)

	// the deleted region is between any regions
	case5InputRv := NewRegionVector()
	case5InputRv.AddRegionNonOverlappingPanic(10, 15)
	case5InputRv.AddRegionNonOverlappingPanic(20, 25)
	case5InputRegion := &Region{15, 20}
	case5ExpectedRv := NewRegionVector()
	case5ExpectedRv.AddRegionNonOverlappingPanic(10, 15)
	case5ExpectedRv.AddRegionNonOverlappingPanic(20, 25)

	// the deleted region overlaps with the first region left
	case6InputRv := NewRegionVector()
	case6InputRv.AddRegionNonOverlappingPanic(10, 15)
	case6InputRv.AddRegionNonOverlappingPanic(20, 25)
	case6InputRegion := &Region{5, 12}
	case6ExpectedRv := NewRegionVector()
	case6ExpectedRv.AddRegionNonOverlappingPanic(12, 15)
	case6ExpectedRv.AddRegionNonOverlappingPanic(20, 25)

	// the deleted region overlaps with the first region right
	case7InputRv := NewRegionVector()
	case7InputRv.AddRegionNonOverlappingPanic(10, 15)
	case7InputRv.AddRegionNonOverlappingPanic(20, 25)
	case7InputRegion := &Region{12, 18}
	case7ExpectedRv := NewRegionVector()
	case7ExpectedRv.AddRegionNonOverlappingPanic(10, 12)
	case7ExpectedRv.AddRegionNonOverlappingPanic(20, 25)

	// the deleted region overlaps with the first region right and the second region left
	case8InputRv := NewRegionVector()
	case8InputRv.AddRegionNonOverlappingPanic(10, 15)
	case8InputRv.AddRegionNonOverlappingPanic(20, 25)
	case8InputRegion := &Region{12, 22}
	case8ExpectedRv := NewRegionVector()
	case8ExpectedRv.AddRegionNonOverlappingPanic(10, 12)
	case8ExpectedRv.AddRegionNonOverlappingPanic(22, 25)

	// the deleted region overlaps with the first region right, contains the second region and
	// overlaps with the third region left
	case9InputRv := NewRegionVector()
	case9InputRv.AddRegionNonOverlappingPanic(10, 15)
	case9InputRv.AddRegionNonOverlappingPanic(20, 25)
	case9InputRv.AddRegionNonOverlappingPanic(30, 35)
	case9InputRegion := &Region{12, 32}
	case9ExpectedRv := NewRegionVector()
	case9ExpectedRv.AddRegionNonOverlappingPanic(10, 12)
	case9ExpectedRv.AddRegionNonOverlappingPanic(32, 35)

	// the deleted region is contained within a region
	case10InputRv := NewRegionVector()
	case10InputRv.AddRegionNonOverlappingPanic(10, 15)
	case10InputRv.AddRegionNonOverlappingPanic(20, 25)
	case10InputRegion := &Region{22, 24}
	case10ExpectedRv := NewRegionVector()
	case10ExpectedRv.AddRegionNonOverlappingPanic(10, 15)
	case10ExpectedRv.AddRegionNonOverlappingPanic(20, 22)
	case10ExpectedRv.AddRegionNonOverlappingPanic(24, 25)

	testCases := []struct {
		name        string
		inputRv     *RegionVector
		inputRegion *Region
		expectedRv  *RegionVector
	}{
		{
			name:        "Empty region vector",
			inputRv:     case1InputRv,
			inputRegion: case1InputRegion,
			expectedRv:  case1ExpectedRv,
		},
		{
			name:        "Deleted region is empty",
			inputRv:     case2InputRv,
			inputRegion: case2InputRegion,
			expectedRv:  case2ExpectedRv,
		},
		{
			name:        "Deleted region is before the first region",
			inputRv:     case3InputRv,
			inputRegion: case3InputRegion,
			expectedRv:  case3ExpectedRv,
		},
		{
			name:        "Deleted region is after the last region",
			inputRv:     case4InputRv,
			inputRegion: case4InputRegion,
			expectedRv:  case4ExpectedRv,
		},
		{
			name:        "Deleted region is between any regions",
			inputRv:     case5InputRv,
			inputRegion: case5InputRegion,
			expectedRv:  case5ExpectedRv,
		},
		{
			name:        "Deleted region overlaps with the first region left",
			inputRv:     case6InputRv,
			inputRegion: case6InputRegion,
			expectedRv:  case6ExpectedRv,
		},
		{
			name:        "Deleted region overlaps with the first region right",
			inputRv:     case7InputRv,
			inputRegion: case7InputRegion,
			expectedRv:  case7ExpectedRv,
		},
		{
			name:        "Deleted region overlaps with the first region right and the second region left",
			inputRv:     case8InputRv,
			inputRegion: case8InputRegion,
			expectedRv:  case8ExpectedRv,
		},
		{
			name:        "Deleted region overlaps with the first region right, contains the second region and overlaps with the third region left",
			inputRv:     case9InputRv,
			inputRegion: case9InputRegion,
			expectedRv:  case9ExpectedRv,
		},
		{
			name:        "Deleted region is contained within a region",
			inputRv:     case10InputRv,
			inputRegion: case10InputRegion,
			expectedRv:  case10ExpectedRv,
		},
	}

	for _, testCase := range testCases {
		t.Run(testCase.name, func(t *testing.T) {

			testCase.inputRv.RemoveRegion(testCase.inputRegion.Start, testCase.inputRegion.End)

			if !testCase.inputRv.Equals(testCase.expectedRv) {
				t.Errorf("expected: %v, got: %v", testCase.expectedRv, testCase.inputRv)
			}
		})
	}
}
