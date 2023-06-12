package mapping

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core/datastructure"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/dataloader/fastq"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
	"math/rand"
)

// Info the mapping information about a kmer
type Info struct {
	IndexReference int // the position of the match in den reference sequence
	IndexRead      int // the position of the kmer in the read
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
func GetAnchorMatches(read *fastq.Read, tree *datastructure.SuffixTree) *map[int][]Info {

	var kmerSize int = 8
	var numKmers int = 5

	var lenRead int = len(read.Sequence)
	startIndices := utils.Arange(0, lenRead-kmerSize, numKmers)

	hits := make(map[int][]Info)

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
			info := Info{
				IndexReference: match.From,
				IndexRead:      nextIndex,
			}

			hits[match.SequenceIndex] = append(hits[match.SequenceIndex], info)
		}
	}

	return &hits
}

func LocateR1PositionOnStrands(gtaIndex *index.GtaIndex, r1Read *fastq.Read) (*map[int][]Info, bool) {

	// matches the kmers to the forward strand
	var hitsFw *map[int][]Info = GetAnchorMatches(r1Read, gtaIndex.SuffixTreeForwardStrandForwardDirection)

	// no kmer matches on the forward strand
	// TODO: create a function to check if there are not sufficiently enough hits on this strand
	if len(*hitsFw) == 0 {

		// matches the kmers to the reverse strand
		var hitsRv *map[int][]Info = GetAnchorMatches(r1Read, gtaIndex.SuffixTreeReverseStrandForwardDirection)

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
	if readPair.ReadR2 != nil {
		rvReadHeader = readPair.ReadR2.Header
	}

	logrus.WithFields(logrus.Fields{
		"fwReadHeader": readPair.ReadR1.Header,
		"rvReadHeader": rvReadHeader,
	}).Info("Map new read pair")

	fmt.Println("R1")
	fmt.Println(readPair.ReadR1.Sequence)
	fmt.Println("R2")
	fmt.Println(readPair.ReadR2.Sequence)

	hitsR1, isForwardStrand := LocateR1PositionOnStrands(gtaIndex, readPair.ReadR1)

	if hitsR1 == nil {
		logrus.Info("Discard read pair because R1 does not match")
		return
	}

	fmt.Println(hitsR1)
	fmt.Println(isForwardStrand)

	// TODO: implement check for single end reads

	var hitsR2 *map[int][]Info

	if isForwardStrand {
		hitsR2 = GetAnchorMatches(readPair.ReadR2, gtaIndex.SuffixTreeForwardStrandReverseDirection)
	} else {
		hitsR2 = GetAnchorMatches(readPair.ReadR2, gtaIndex.SuffixTreeReverseStrandReverseDirection)
	}

	if hitsR2 == nil {
		logrus.Info("Discard read pair because R2 does not match")
		return
	}

	fmt.Println(hitsR2)

	for index1, _ := range *hitsR1 {
		for index2, _ := range *hitsR2 {
			if index1 == index2 {
				fmt.Println(gtaIndex.Transcripts[index1].TranscriptIdEnsembl)
				fmt.Println(index1)
			}
		}
	}
}
