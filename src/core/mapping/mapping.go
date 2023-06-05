package mapping

import (
	"github.com/KleinSamuel/gtamap/src/core/datastructure"
	"github.com/KleinSamuel/gtamap/src/core/index"
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

func LocateR1PositionOnStrands(gtaIndex *index.GtaIndex, r1Read *fastq.Read) (*map[Position][]int, bool) {

	// matches the kmers to the forward strand
	var hitsFw *map[Position][]int = GetAnchorMatches(r1Read, gtaIndex.SuffixTreeFw)

	// no kmer matches on the forward strand
	if len(*hitsFw) == 0 {

		// matches the kmers to the reverse strand
		var hitsRv *map[Position][]int = GetAnchorMatches(r1Read, gtaIndex.SuffixTreeRv)

		// no kmer matches on the reverse strand
		if len(*hitsRv) == 0 {
			// discard read pair
			return nil, false
		}

		return hitsRv, false
	}

	return hitsFw, true
}

func MapReadPair(readPair *fastq.ReadPair, gtaIndex *index.GtaIndex) {

	rvReadHeader := ""
	if readPair.RvRead != nil {
		rvReadHeader = readPair.RvRead.Header
	}

	logrus.WithFields(logrus.Fields{
		"fwReadHeader": readPair.FwRead.Header,
		"rvReadHeader": rvReadHeader,
	}).Info("Map new read pair")

	// the kmer matches on the forward strand
	hitsR1, isForwardStrand := LocateR1PositionOnStrands(gtaIndex, readPair.FwRead)

	logrus.Debug("hitsR1: ", hitsR1)
	logrus.Debug("isForwardStrand: ", isForwardStrand)

	//fmt.Println(hitsR1)
	//fmt.Println(isForwardStrand)
}
