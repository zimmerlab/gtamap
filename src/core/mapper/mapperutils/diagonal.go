package mapperutils

import (
	"github.com/sirupsen/logrus"
)

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

func (dh *DiagonalHandler) ConsumeKmer(kmerStart int, kmerStop int, kmerStartGenome int, kmerStopGenome int) {

	logrus.WithFields(logrus.Fields{
		"kmerStart": kmerStart,
		"kmerStop":  kmerStop,
	}).Debug("consuming kmer after gap filling")

	dh.ConsumeRegionRead(kmerStart, kmerStop)
	dh.ConsumeRegionGenome(kmerStartGenome, kmerStopGenome)
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
