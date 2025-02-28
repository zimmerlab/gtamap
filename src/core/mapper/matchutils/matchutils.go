package matchutils

import "github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"

type Match struct {
	SequenceIndex int // the index of the sequence in the genome
	FromGenome    int // the start position of the match in the genome
	ToGenome      int // the end position of the match in the genome
	FromRead      int // the start position of the match in the read
	ToRead        int // the end position of the match in the read
	StartGenome   int // the start position of the match in the genome (diagonal)
}

type GlobalMatchResult struct {
	MatchesPerSequence []*SequenceMatchResult
}

type SequenceMatchResult struct {
	MatchesPerDiagonal map[int][]*Match
}

type ReadMatchResult struct {
	SequenceIndex  int                        // the index of the sequence in the genome
	MatchedRead    *regionvector.RegionVector // region vector containing the matched positions in the read
	MatchedGenome  *regionvector.RegionVector // region vector containing the matched positions in the genome
	MismatchesRead []int                      // the positions of the mismatches in the read
}
