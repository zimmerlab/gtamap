package mapper

import (
	"sort"
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
