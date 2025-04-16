package mapperutils

import (
	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"strconv"
	"strings"
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
	SecondPass     bool                       // true if this result must undergo a second pass
}

func (m ReadMatchResult) GetCigar() string {
	var builder strings.Builder

	isForwardStrand := m.SequenceIndex == 0

	if isForwardStrand {

		index := 0

		for index < len(m.MatchedGenome.Regions) {
			numMatches := m.MatchedGenome.Regions[index].End - m.MatchedGenome.Regions[index].Start
			builder.WriteString(strconv.Itoa(numMatches))
			builder.WriteString("M")

			if index+1 < len(m.MatchedGenome.Regions) {
				numSkipped := m.MatchedGenome.Regions[index+1].Start - m.MatchedGenome.Regions[index].End
				builder.WriteString(strconv.Itoa(numSkipped))

				if numSkipped < config.IntronLengthMin() {
					builder.WriteString("D")
				} else {
					builder.WriteString("N")
				}
			}

			index++
		}
	} else {

		index := len(m.MatchedGenome.Regions) - 1

		for index >= 0 {
			numMatches := m.MatchedGenome.Regions[index].End - m.MatchedGenome.Regions[index].Start
			builder.WriteString(strconv.Itoa(numMatches))
			builder.WriteString("M")

			if index-1 >= 0 {
				numSkipped := m.MatchedGenome.Regions[index].Start - m.MatchedGenome.Regions[index-1].End
				builder.WriteString(strconv.Itoa(numSkipped))

				if numSkipped < config.IntronLengthMin() {
					builder.WriteString("D")
				} else {
					builder.WriteString("N")
				}
			}

			index--
		}
	}

	return builder.String()
}
