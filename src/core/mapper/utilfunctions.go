package mapper

import "sort"

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

func spliceSiteDonorScore(first byte, second byte, isForwardStrand bool) int {
	if isForwardStrand {
		if first == byte('G') && second == byte('T') {
			// canonical splice donor site dinucleotide GT
			return 0
		} else if first == byte('G') && second == byte('C') {
			// non-canonical splice donor site dinucleotide GC
			return 1
		} else if first == byte('A') && second == byte('T') {
			// non-canonical splice donor site dinucleotide GC
			return 1
		} else {
			// all other splice donor site dinucleotides
			return 2
		}
	} else {
		if first == byte('C') && second == byte('T') {
			// canonical splice acceptor site dinucleotide AG
			return 0
		} else if first == byte('G') && second == byte('T') {
			return 1
		} else {
			return 2
		}
	}
}

func spliceSiteAcceptorScore(first byte, second byte, isForwardStrand bool) int {
	if isForwardStrand {
		if first == byte('A') && second == byte('G') {
			// canonical splice acceptor site dinucleotide AG
			return 0
		} else if first == byte('A') && second == byte('C') {
			return 1
		} else {
			return 2
		}
	} else {
		if first == byte('A') && second == byte('C') {
			// canonical splice donor site dinucleotide GT
			return 0
		} else if first == byte('G') && second == byte('C') {
			// non-canonical splice donor site dinucleotide GC
			return 1
		} else if first == byte('A') && second == byte('T') {
			// non-canonical splice donor site dinucleotide GC
			return 1
		} else {
			// all other splice donor site dinucleotides
			return 2
		}
	}
}
