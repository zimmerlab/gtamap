package matchutils

import (
	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/sirupsen/logrus"
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
				builder.WriteString("N")
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
				builder.WriteString("N")
			}

			index--
		}
	}

	return builder.String()
}

/* DIAGONAL HANDLER */

type DiagonalHandler struct {
	Diagonals map[int][]*Match
}

func NewDiagonalHandler() *DiagonalHandler {
	return &DiagonalHandler{
		Diagonals: make(map[int][]*Match),
	}
}

func NewDiagonalHandlerWithData(m map[int][]*Match) *DiagonalHandler {
	return &DiagonalHandler{
		Diagonals: m,
	}
}

func (dh *DiagonalHandler) GetBestDiagonal() (int, int) {
	indexMax := -1
	maxMatches := 0
	for key, matches := range dh.Diagonals {
		countUnused := 0
		for _, match := range matches {
			if match.Used {
				continue
			}
			countUnused++
		}
		if countUnused > maxMatches {
			indexMax = key
			maxMatches = countUnused
		}
	}
	return indexMax, maxMatches
}

func (dh *DiagonalHandler) ConsumeDiagonal(diagonal int) {

	logrus.WithFields(logrus.Fields{
		"diagonal": diagonal,
	}).Debug("consuming diagonal")

	for _, match := range dh.Diagonals[diagonal] {
		dh.ConsumeRegionRead(match.FromRead, match.ToRead)
		dh.ConsumeRegionGenome(match.FromGenome, match.ToGenome)
	}
	delete(dh.Diagonals, diagonal)
}

func (dh *DiagonalHandler) ConsumeRegionRead(startRead int, endRead int) {

	logrus.WithFields(logrus.Fields{
		"start": startRead,
		"end":   endRead,
	}).Debug("consuming match based on read region")

	for _, matches := range dh.Diagonals {
		for _, match := range matches {
			if match.Used {
				continue
			}
			// the match overlaps an already used region on the read
			if (match.FromRead >= startRead && match.FromRead < endRead) ||
				(match.ToRead > startRead && match.ToRead <= endRead) {
				match.Used = true

				logrus.WithFields(logrus.Fields{
					"match": match,
				}).Debug("set used")
			}
		}
	}
}
func (dh *DiagonalHandler) ConsumeRegionGenome(startGenome int, endGenome int) {

	logrus.WithFields(logrus.Fields{
		"start": startGenome,
		"end":   endGenome,
	}).Debug("consuming match based on genome region")

	for _, matches := range dh.Diagonals {
		for _, match := range matches {
			if match.Used {
				continue
			}
			// the match overlaps an already used region on the genome
			if (match.FromGenome >= startGenome && match.FromGenome < endGenome) ||
				(match.ToGenome > startGenome && match.ToGenome <= endGenome) {
				match.Used = true

				logrus.WithFields(logrus.Fields{
					"match": match,
				}).Debug("set used")
			}
		}
	}
}
