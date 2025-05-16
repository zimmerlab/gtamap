package mapper

import (
	"fmt"
	"sort"
	"time"
)

// sortedIndicesAsc sorts the indices of a list in ascending order based on the values at those indices.
//
// Example:
// list := []int{3, 1, 2}
// sortedIndices := sortedIndicesAsc(list)
// fmt.Println(sortedIndices) // Output: [1 2 0]
func sortedIndicesDesc(list []int) []int {
	indices := make([]int, len(list))
	for i := range list {
		indices[i] = i
	}
	sort.Slice(indices, func(i, j int) bool {
		return list[indices[i]] > list[indices[j]]
	})
	return indices
}

func scoreSpliceSites(donorFirstBase byte, donorSecondBase byte, acceptorFirstBase byte,
	acceptorSecondBase byte, isForwardStrand bool) int {

	if isForwardStrand {

		if donorFirstBase == byte('G') && donorSecondBase == byte('T') &&
			acceptorFirstBase == byte('A') && acceptorSecondBase == byte('G') {
			// canonical splice site GT/AG
			return 0
		} else if donorFirstBase == byte('G') && donorSecondBase == byte('C') &&
			acceptorFirstBase == byte('A') && acceptorSecondBase == byte('G') {
			// non-canonical splice site GC/AG
			return 1
		} else if donorFirstBase == byte('A') && donorSecondBase == byte('T') &&
			acceptorFirstBase == byte('A') && acceptorSecondBase == byte('C') {
			// non-canonical splice site AT/AC
			return 1
		}

	} else {

		if donorFirstBase == byte('C') && donorSecondBase == byte('T') &&
			acceptorFirstBase == byte('A') && acceptorSecondBase == byte('C') {
			// canonical splice site GT/AG on rev strand
			return 0
		} else if donorFirstBase == byte('C') && donorSecondBase == byte('T') &&
			acceptorFirstBase == byte('G') && acceptorSecondBase == byte('C') {
			// non-canonical splice site GC/AG on rev strand
			return 1
		} else if donorFirstBase == byte('G') && donorSecondBase == byte('T') &&
			acceptorFirstBase == byte('A') && acceptorSecondBase == byte('T') {
			// non-canonical splice site AT/AC on rev strand
			return 1
		}
	}

	// all other non-canonical splice sites
	return 2
}

func FormatDuration(d time.Duration) string {
	if d.Hours() < 1 {
		if d.Minutes() < 1 {
			if d.Seconds() < 1 {
				return fmt.Sprintf("%dms", int(d.Milliseconds()))
			} else {
				return fmt.Sprintf("%ds", int(d.Seconds()))
			}
		} else {
			return fmt.Sprintf("%dm %ds", int(d.Minutes()), int(d.Seconds())%60)
		}
	} else {
		return fmt.Sprintf("%dh %dm %ds", int(d.Hours()), int(d.Minutes())%60, int(d.Seconds())%60)
	}
}
