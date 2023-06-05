package mapping

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core/datastructure"
	"github.com/KleinSamuel/gtamap/src/dataloader/fastq"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
	"math/rand"
)

type Position struct {
	SequenceIndex int
	Start         int
}
type PositionCount struct {
	Position *Position
	Count    int
}

func drawNextKmerIndex(startIndices []int) int {
	nextIndex := rand.Intn(len(startIndices))

	// removes the drawn index from the list of possible indices
	startIndices[nextIndex] = startIndices[len(startIndices)-1]
	startIndices = startIndices[:len(startIndices)-1]

	return nextIndex
}

// GetAnchorMatches returns a map of positions in the transcriptome where the read could be aligned to.
// Kmers are generated from the read and searched in the suffix tree.
// The kmers are taken from evenly spaced starting positions in the read.
// TODO: not only keep track of the position in the transcript where the kmer was mapped but also the position in the read
func GetAnchorMatches(read *fastq.Read, tree *datastructure.SuffixTree) *map[Position][]int {

	var kmerSize int = 8
	var numKmers int = 5

	var lenRead int = len(read.Sequence)
	startIndices := utils.Arange(0, lenRead-kmerSize, numKmers)

	hits := make(map[Position][]int)

	for i := 0; i < numKmers; i++ {

		nextIndex := startIndices[i]

		kmer := read.Sequence[nextIndex : nextIndex+kmerSize]

		result := tree.Search(&kmer)

		if result == nil {
			continue
		}

		for _, match := range result.Matches {

			matchIndex := match.From - nextIndex

			// read and transcript sequence can not be aligned
			if matchIndex < 0 {
				continue
			}

			// TODO: (performance) use simple pair (array) instead of struct
			hit := Position{
				SequenceIndex: match.SequenceIndex,
				Start:         matchIndex,
			}

			hits[hit] = append(hits[hit], match.From)
		}
	}

	return &hits
}

func LocateReadPositions() {

}

func MapReadPair(readPair *fastq.ReadPair, tree *datastructure.SuffixTree) {

	rvReadHeader := ""
	if readPair.RvRead != nil {
		rvReadHeader = readPair.RvRead.Header
	}

	logrus.WithFields(logrus.Fields{
		"fwReadHeader": readPair.FwRead.Header,
		"rvReadHeader": rvReadHeader,
	}).Info("Map read pair")

	var hitsR1 *map[Position][]int = GetAnchorMatches(readPair.FwRead, tree)
	var hitsR2 *map[Position][]int = nil

	if len(*hitsR1) == 0 {
		fmt.Println("no hits for R1 read")

		if readPair.RvRead == nil {
			fmt.Println("no R2 read -> discard pair")
			return
		}

		fmt.Println("check R2 read")

		hitsR2 = GetAnchorMatches(readPair.RvRead, tree)

		if len(*hitsR2) == 0 {
			fmt.Println("no hits for R2 read -> discard pair")
			return
		}
	}

	fmt.Println("hits R1: ", *hitsR1)
	fmt.Println("hits R2: ", *hitsR2)

}
