package mapping

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/algorithms"
	"github.com/KleinSamuel/gtamap/src/core/datastructure"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
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

func addKmerMatchToCigar(cigarList *[]rune) *[]rune {
	for i := 0; i < config.GetKmerLength(); i++ {
		*cigarList = append(*cigarList, 'M')
	}
	return cigarList
}

func finalizeCigar(cigarList *[]rune) string {
	cigarString := ""

	lastRune := (*cigarList)[0]
	occurence := 1

	for i := 1; i < len(*cigarList); i++ {

		currentRune := (*cigarList)[i]
		if currentRune == lastRune {
			occurence++
		} else {
			cigarString += fmt.Sprintf("%d%c", occurence, lastRune)
			occurence = 1
			lastRune = currentRune
		}
	}
	cigarString += fmt.Sprintf("%d%c", occurence, lastRune)

	return cigarString
}

func AlignRead(read *fastq.Read, transcriptSequence *string, tree *datastructure.SuffixTree, kmerHitList []Info) {

	leftmostMappingPosition := -1
	cigarList := make([]rune, 0)

	for i := 0; i < len(kmerHitList); i++ {

		hit := kmerHitList[i]

		fmt.Println(i, hit)

		if i == 0 {
			// the first kmer does not start at the beginning of the read
			if hit.IndexRead != 0 {
				fmt.Println("first kmer does not start at the beginning of the read")

				// the region before the first kmer
				startPrefixGapRead := 0
				endPrefixGapRead := hit.IndexRead
				lenPrefixGapRead := endPrefixGapRead - startPrefixGapRead

				startPrefixGapRef := hit.IndexReference - lenPrefixGapRead
				endPrefixGapRef := hit.IndexReference

				score, cigarPart, seq1, seq2 := algorithms.NeedlemanWunsch((*transcriptSequence)[startPrefixGapRef:endPrefixGapRef], read.Sequence[startPrefixGapRead:endPrefixGapRead])

				fmt.Println(score)
				fmt.Println(cigarPart)
				fmt.Println(seq1)
				fmt.Println(seq2)

				// add the alignment of the region before the first kmer to the cigar
				cigarList = append(cigarList, cigarPart...)

				// TODO: determine leftmost mapping position
			}

			// sets the leftmost mapping position if not already set by gap before first kmer
			if leftmostMappingPosition == -1 {
				leftmostMappingPosition = hit.IndexReference
			}

			// add the kmer match to the cigar
			addKmerMatchToCigar(&cigarList)

			continue
		}

		// check gap between current and previous kmer in read
		startGapRead := kmerHitList[i-1].IndexRead + config.GetKmerLength()
		endGapRead := hit.IndexRead
		lenGapRead := endGapRead - startGapRead
		// check gap between current and previous kmer in reference
		startGapRef := kmerHitList[i-1].IndexReference + config.GetKmerLength()
		endGapRef := hit.IndexReference
		lenGapRef := endGapRef - startGapRef

		if lenGapRead == lenGapRef {
			// TODO: (performance) check if there is any smarter way of comparing strings of same length
		}

		score, cigarPart, seq1, seq2 := algorithms.NeedlemanWunsch((*transcriptSequence)[startGapRef:endGapRef], read.Sequence[startGapRead:endGapRead])

		fmt.Println(score)
		fmt.Println(cigarPart)
		fmt.Println(seq1)
		fmt.Println(seq2)

		// add the alignment to the cigar string
		cigarList = append(cigarList, cigarPart...)

		// add the kmer match to the cigar string
		addKmerMatchToCigar(&cigarList)

		// the last kmer does not end at the end of the read
		if i == len(kmerHitList)-1 && hit.IndexRead < len(read.Sequence)-config.GetKmerLength() {

			fmt.Println("last kmer does not end at the end of the read")

			startSuffixGapRead := hit.IndexRead + config.GetKmerLength()
			endSuffixGapRead := len(read.Sequence)
			lenSuffixGapRead := endSuffixGapRead - startSuffixGapRead

			startSuffixGapRef := hit.IndexReference + config.GetKmerLength()
			endSuffixGapRef := startSuffixGapRef + lenSuffixGapRead

			score, cigarPart, seq1, seq2 := algorithms.NeedlemanWunsch((*transcriptSequence)[startSuffixGapRef:endSuffixGapRef], read.Sequence[startSuffixGapRead:endSuffixGapRead])

			fmt.Println(score)
			fmt.Println(cigarPart)
			fmt.Println(seq1)
			fmt.Println(seq2)

			// add the alignment of the suffix to the cigar string
			cigarList = append(cigarList, cigarPart...)
		}
	}

	fmt.Println("leftmostMappingPosition", leftmostMappingPosition)
	fmt.Println("cigarString", cigarList)
	fmt.Println(len(cigarList))

	fmt.Println(finalizeCigar(&cigarList))
}

func MapReadPair(readPair *fastq.ReadPair, gtaIndex *index.GtaIndex) string {

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
		return ""
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
		return ""
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

	return "mapping result dummy"
}
