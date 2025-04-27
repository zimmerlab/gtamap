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

func (dh *DiagonalHandler) IsValidExtension(possibleExtension []*Match, result ReadMatchResult) bool {
	if len(result.MatchedRead.Regions) == 0 {
		return true
	}

	// get first and last match from diagonal
	firstMatch := possibleExtension[0]
	lastMatch := possibleExtension[len(possibleExtension)-1]

	if firstMatch.FromRead >= result.MatchedRead.GetLastRegion().End {
		// right ext
		// the genomic coords of te ext must be grater that the right most genomic coords of the
		// result
		if firstMatch.FromGenome >= result.MatchedGenome.GetLastRegion().End {
			// valid
			return true
		}
		// invalid ext
		return false
	}

	if lastMatch.ToRead <= result.MatchedRead.GetFirstRegion().Start {
		// left ext
		// the genomic coords of te ext must be smaller that the left most genomic coords of the
		if firstMatch.ToGenome <= result.MatchedGenome.GetFirstRegion().Start {
			// valid ext
			return true
		}
		// invalid ext
		return false
	}

	// if there are two or more regions already mapped in result
	// we need to check if the extension is placed in between
	if len(result.MatchedRead.Regions) > 1 {
		for i := 0; i < len(result.MatchedRead.Regions)-1; i++ {
			readGapStart := result.MatchedRead.Regions[i].End + 1
			readGapStop := result.MatchedRead.Regions[i+1].Start - 1
			geneGapStart := result.MatchedGenome.Regions[i].End + 1
			geneGapStop := result.MatchedGenome.Regions[i+1].Start - 1

			// do the read coordinates of the possible extension fit in between
			// the already used kmers?
			if firstMatch.FromRead >= readGapStart && lastMatch.ToRead <= readGapStop {
				// if yes, then the gene coords need to also be in between the already consumed
				// gene coords
				if firstMatch.FromGenome >= geneGapStart && lastMatch.ToGenome <= geneGapStop {
					// we dont expect this to happen.
					return true
				}
			}
			// but this happens a lot where an extension would be a mid extension and we
			// most of the time don't want this (since it is not possible)
			return false

		}
	}
	return true
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
