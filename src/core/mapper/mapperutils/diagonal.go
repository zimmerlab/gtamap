package mapperutils

import (
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
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

// IsValidExtension checks if the diagonal is a valid extension of the already mapped regions and return true if it is.
// A diagonal is not valid if its start and end positions are not consistent with the already mapped regions.
// The diagonal must be between the already mapped regions in the read and genome.
// An error is thrown in the given diagonal contains a kmer that was already used.
//
// An example of a valid extension is:
// Already mapped regions: [0, 10], [20, 30] to genome [1000, 1010], [1020, 1030]
// New diagonal: [10, 20] to genome [1010, 1020]
// An example of an invalid extension is:
// Already mapped regions: [0, 10], [20, 30] to genome [1000, 1010], [1020, 1030]
// New diagonal: [10, 20] to genome [900, 910]
func (dh *DiagonalHandler) IsValidExtension(possibleExtension []*Match, result ReadMatchResult, read *fastq.Read) bool {

	// if there is nothing matched in the read, then any extension is valid
	if len(result.MatchedRead.Regions) == 0 {
		return true
	}

	// the min and max positions of the new diagonal in both the read and the genome
	minRead := -1
	maxRead := -1
	minGenome := -1
	maxGenome := -1

	// only use regions that are not used yet
	for _, match := range possibleExtension {
		if match.Used {
			logrus.WithFields(logrus.Fields{
				"match": match.String(),
				"read":  read.Header,
			}).Warn("tried to use a diagonal containing a kmer that was already used")
			continue
		}

		if minRead == -1 || match.FromRead < minRead {
			minRead = match.FromRead
		}
		if maxRead == -1 || match.ToRead > maxRead {
			maxRead = match.ToRead
		}
		if minGenome == -1 || match.FromGenome < minGenome {
			minGenome = match.FromGenome
		}
		if maxGenome == -1 || match.ToGenome > maxGenome {
			maxGenome = match.ToGenome
		}
	}

	// check if the diagonal position within the genome is still consistent with the already mapped regions
	for i, resultMatch := range result.MatchedRead.Regions {
		if resultMatch.End < minRead {
			// the existing match is before the new match in the read
			// then its diagonal must be before the new diagonal
			if minGenome < result.MatchedGenome.Regions[i].End {
				return false
			}
		} else if resultMatch.Start > maxRead {
			// the existing match is after the new match in the read
			// then its diagonal must be after the new diagonal
			if maxGenome > result.MatchedGenome.Regions[i].Start {
				return false
			}
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
