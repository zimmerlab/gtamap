package mapperutils

import (
	"math"
	"sort"

	"github.com/KleinSamuel/gtamap/src/formats/fastq"
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

func NewDiagonalHandlerWithDataCopy(m map[int][]*Match) *DiagonalHandler {
	mapCopy := make(map[int][]*Match, len(m))

	// for key, matches := range m {
	// 	mapCopy[key] = make([]*Match, len(matches))
	// 	for i, match := range matches {
	// 		mapCopy[key][i] = &Match{
	// 			FromRead:   match.FromRead,
	// 			ToRead:     match.ToRead,
	// 			FromGenome: match.FromGenome,
	// 			ToGenome:   match.ToGenome,
	// 			Used:       match.Used,
	// 		}
	// 	}
	// }
	for key, matches := range m {
		copiedMatches := make([]*Match, len(matches))
		for i, match := range matches {
			copiedMatches[i] = &Match{
				FromRead:   match.FromRead,
				ToRead:     match.ToRead,
				FromGenome: match.FromGenome,
				ToGenome:   match.ToGenome,
				Used:       match.Used,
			}
		}
		mapCopy[key] = copiedMatches
	}

	return &DiagonalHandler{
		Diagonals: mapCopy,
	}
}

// GetBestDiagonal returns the diagonal with the most unused matches.
// It returns the position of the diagonal in the genome and the number of unused matches if a diagonal is found.
// If not, then it returns -1, -1 and false.
func (dh *DiagonalHandler) GetBestDiagonal() (int, int, bool) {
	indexMax := math.MaxInt32
	maxMatches := 0
	found := false

	for key, matches := range dh.Diagonals {

		// number of exact kmer matches in this diagonal that were not used by any other match yet
		countUnused := 0
		for _, match := range matches {
			if match.Used {
				continue
			}
			countUnused++
		}

		if countUnused >= maxMatches {
			if countUnused == maxMatches && key < indexMax {
				indexMax = key
				maxMatches = countUnused
				found = true
			} else {
				indexMax = key
				maxMatches = countUnused
				found = true
			}
		}

	}

	return indexMax, maxMatches, found
}

func (dh *DiagonalHandler) GetAvailableDiagonals() []int {
	// get sorted list of dh.Diagonals keys
	diagonalKeys := make([]int, 0, len(dh.Diagonals))
	for key := range dh.Diagonals {
		diagonalKeys = append(diagonalKeys, key)
	}

	// sort the keys in ascending order
	sort.Slice(diagonalKeys, func(i, j int) bool {
		return diagonalKeys[i] < diagonalKeys[j]
	})

	return diagonalKeys
}

func (dh *DiagonalHandler) ConsumeDiagonal(diagonal int) {
	// logrus.WithFields(logrus.Fields{
	// 	"diagonal": diagonal,
	// }).Debug("consuming diagonal")

	for _, match := range dh.Diagonals[diagonal] {
		// dh.ConsumeRegionRead(match.FromRead, match.ToRead)
		// dh.ConsumeRegionGenome(match.FromGenome, match.ToGenome)
		dh.ConsumeKmer(match.FromRead, match.ToRead, match.FromGenome, match.ToGenome)
	}

	// delete the given diagonal from the map
	delete(dh.Diagonals, diagonal)
}

// RemoveDiagonal removes the diagonal from the map without consuming it.
// This is used when a diagonal is not valid anymore and should be removed from the map without affecting
// the other diagonals as its matches are not applied.
func (dh *DiagonalHandler) RemoveDiagonal(diagonal int) {
	// logrus.WithFields(logrus.Fields{
	// 	"diagonal": diagonal,
	// }).Debug("removing diagonal")

	// delete the given diagonal from the map
	delete(dh.Diagonals, diagonal)
}

func (dh *DiagonalHandler) ConsumeKmer(kmerStart int, kmerStop int, kmerStartGenome int, kmerStopGenome int) {
	// logrus.WithFields(logrus.Fields{
	// 	"kmerStart": kmerStart,
	// 	"kmerStop":  kmerStop,
	// }).Debug("consuming kmer after gap filling")

	dh.ConsumeRegionRead(kmerStart, kmerStop)
	dh.ConsumeRegionGenome(kmerStartGenome, kmerStopGenome)
}

func (dh *DiagonalHandler) RemovedConsumedRegionsAndDiagonals() {
	// logrus.Debug("removing consumed regions and diagonals")

	// delete the matches that were marked as used
	for diagonal, matches := range dh.Diagonals {

		n := 0
		for _, match := range matches {
			if match.Used {
				//logrus.WithFields(logrus.Fields{
				//	"match":    match,
				//	"diagonal": diagonal,
				//}).Debug("removed match from diagonal")

				continue
			}
			matches[n] = match
			n++
		}

		// update the slice to only contain the matches that were not used
		dh.Diagonals[diagonal] = matches[:n]
	}

	// delete all other diagonals which do not contain any unused matches
	var keysToDelete []int

	for key, matches := range dh.Diagonals {
		if len(matches) == 0 {
			keysToDelete = append(keysToDelete, key)
		}
	}

	for _, key := range keysToDelete {
		//logrus.WithFields(logrus.Fields{
		//	"diagonal": key,
		//}).Debug("removed diagonal")

		delete(dh.Diagonals, key)
	}
}

func (dh *DiagonalHandler) RemoveInvalidDiagonals(result *ReadMatchResult, read *fastq.Read) {
	// logrus.Debug("removing invalid diagonals")

	var keysToDelete []int

	// remove diagonals that are not valid anymore
	for diagonal, matches := range dh.Diagonals {
		if !dh.IsValidExtension(matches, *result, read) {
			keysToDelete = append(keysToDelete, diagonal)
		}
	}

	for _, key := range keysToDelete {
		//logrus.WithFields(logrus.Fields{
		//	"diagonal": key,
		//}).Debug("removed diagonal")

		delete(dh.Diagonals, key)
	}
}

// IsValidExtension checks if the diagonal is a valid extension of the already mapped regions and return true if it is.
// A diagonal is not valid if its start and end positions are not consistent with the already mapped regions.
// The diagonal must be between the already mapped regions in the read and genome.
// Kmer matches that are already used are ignored for this check.
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
			// logrus.WithFields(logrus.Fields{
			// 	"match": match.String(),
			// 	"read":  read.Header,
			// }).Debug("tried to use a diagonal containing a kmer that was already used")
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

	minReadResult := -1
	maxReadResult := -1
	minGenomeResult := -1
	maxGenomeResult := -1

	// get the min and max positions of the already mapped regions in the read and genome
	for _, resultMatch := range result.MatchedRead.Regions {
		if minReadResult == -1 || resultMatch.Start < minReadResult {
			minReadResult = resultMatch.Start
		}
		if maxReadResult == -1 || resultMatch.End > maxReadResult {
			maxReadResult = resultMatch.End
		}
	}
	for _, resultMatch := range result.MatchedGenome.Regions {
		if minGenomeResult == -1 || resultMatch.Start < minGenomeResult {
			minGenomeResult = resultMatch.Start
		}
		if maxGenomeResult == -1 || resultMatch.End > maxGenomeResult {
			maxGenomeResult = resultMatch.End
		}
	}

	// Commented out since it prevents adding a thrid diag to res, if rank of third diag is 2 (diag inbetween already mapped diags)
	// // the diagonal regions on the read can not overlap any other region that was already mapped
	// if !(maxRead <= minReadResult || minRead >= maxReadResult) {
	// 	return false
	// }
	//
	// // the diagonal regions on the genome can not overlap any other region that was already mapped
	// if !(maxGenome <= minGenomeResult || minGenome >= maxGenomeResult) {
	// 	return false
	// }

	// OLD
	// the diagonal regions on the read can not overlap any other region that was already mapped
	//if !(minRead <= minReadResult && maxRead <= minReadResult) || !(minRead <= maxReadResult && maxRead <= maxReadResult) {
	//	return false
	//}
	//// the diagonal regions on the genome can not overlap any other region that was already mapped
	//if !(minGenome <= minGenomeResult && maxGenome <= minGenomeResult) || !(minGenome <= maxGenomeResult && maxGenome <= maxGenomeResult) {
	//	return false
	//}

	// check if the diagonal position within the genome is still consistent with the already mapped regions
	for i, resultMatch := range result.MatchedRead.Regions {
		if resultMatch.End <= minRead {
			// the existing match is before the new match in the read
			// then its diagonal must be before the new diagonal
			if minGenome < result.MatchedGenome.Regions[i].End {
				return false
			}
		} else if resultMatch.Start >= maxRead {
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
	// logrus.WithFields(logrus.Fields{
	// 	"start": startRead,
	// 	"end":   endRead,
	// }).Debug("consuming match based on read region")

	for _, matches := range dh.Diagonals {
		for _, match := range matches {
			if match.Used {
				continue
			}
			// the match overlaps an already used region on the read
			if (match.FromRead >= startRead && match.FromRead < endRead) ||
				(match.ToRead > startRead && match.ToRead <= endRead) {
				match.Used = true

				//logrus.WithFields(logrus.Fields{
				//	"match": match,
				//}).Debug("set used")
			}
		}
	}
}

func (dh *DiagonalHandler) ConsumeRegionGenome(startGenome int, endGenome int) {
	// logrus.WithFields(logrus.Fields{
	// 	"start": startGenome,
	// 	"end":   endGenome,
	// }).Debug("consuming match based on genome region")

	for _, matches := range dh.Diagonals {
		for _, match := range matches {
			if match.Used {
				continue
			}
			// the match overlaps an already used region on the genome
			if (match.FromGenome >= startGenome && match.FromGenome < endGenome) ||
				(match.ToGenome > startGenome && match.ToGenome <= endGenome) {
				match.Used = true

				//logrus.WithFields(logrus.Fields{
				//	"match": match,
				//}).Debug("set used")
			}
		}
	}
}

func (dh *DiagonalHandler) Copy() *DiagonalHandler {
	return NewDiagonalHandlerWithDataCopy(dh.Diagonals)
}
