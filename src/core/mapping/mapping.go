package mapping

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/config"
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
func GetAnchorMatches(read *fastq.Read, tree *datastructure.SuffixTree) *map[int][]Info {

	var lenRead int = len(read.Sequence)
	startIndices := utils.Arange(0, lenRead-config.GetKmerLength(), config.GetNumKmers())

	hits := make(map[int][]Info)

	for i := 0; i < config.GetNumKmers(); i++ {

		nextIndex := startIndices[i]

		kmer := read.Sequence[nextIndex : nextIndex+config.GetKmerLength()]

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
			// TODO: check if match index is too far right s.t. the fragment would extend beyond the read

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

func AlignRead(read *fastq.Read, transcriptSequence *string, tree *datastructure.SuffixTree, kmerHitList []Info) {

	fmt.Println(read.Sequence)
	fmt.Println(*transcriptSequence)

	for i, hit := range kmerHitList {

		// the first kmer does not start at the beginning of the read
		if i == 0 && hit.IndexRead != 0 {
			fmt.Println("first kmer does not start at the beginning of the read")

			// the region before the first kmer
			start := 0
			end := hit.IndexRead

			fmt.Println(start, end)

			continue
		}

		fmt.Println(hit)
	}
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

	// TODO: implement check for single end reads

	var hitsR2 *map[int][]Info

	if isForwardStrand {
		hitsR2 = GetAnchorMatches(readPair.ReadR2, gtaIndex.SuffixTreeReverseStrandForwardDirection)
	} else {
		hitsR2 = GetAnchorMatches(readPair.ReadR2, gtaIndex.SuffixTreeForwardStrandForwardDirection)
	}

	if hitsR2 == nil {
		logrus.Info("Discard read pair because R2 does not match")
		return
	}

	// determine most likely transcript by using the count of kmer hits

	maxHits := 0
	idTranscript := -1

	for index1, _ := range *hitsR1 {
		for index2, _ := range *hitsR2 {
			if index1 == index2 {

				// TODO: filter unique hits per transcript

				// TODO: find a better way to determine the most likely transcript (maybe use fragment length)

				if len((*hitsR1)[index1])+len((*hitsR2)[index2]) > maxHits {
					maxHits = len((*hitsR1)[index1]) + len((*hitsR2)[index2])
					idTranscript = index1
				}
			}
		}
	}

	fmt.Println("idTranscript", idTranscript)
	fmt.Println("maxHits", maxHits)

	if isForwardStrand {
		AlignRead(readPair.ReadR1, &gtaIndex.Transcripts[idTranscript].SequenceDnaForwardStrandForwardDirection,
			gtaIndex.SuffixTreeForwardStrandForwardDirection, (*hitsR1)[idTranscript])
	} else {
		AlignRead(readPair.ReadR1, &gtaIndex.Transcripts[idTranscript].SequenceDnaReverseStrandForwardDirection,
			gtaIndex.SuffixTreeReverseStrandForwardDirection, (*hitsR1)[idTranscript])
	}
}
