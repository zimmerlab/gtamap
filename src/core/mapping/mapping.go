package mapping

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/algorithms"
	"github.com/KleinSamuel/gtamap/src/core/datastructure"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
	"math/rand"
	"strings"
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
	startIndices := utils.Arange(0, lenRead-config.KmerLength(), config.NumKmers())

	hits := make(map[int][]Info)

	for i := 0; i < config.NumKmers(); i++ {

		nextIndex := startIndices[i]

		kmer := read.Sequence[nextIndex : nextIndex+config.KmerLength()]

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

		r1Read.Sequence = utils.ReverseComplementDNA(r1Read.Sequence)

		// matches the kmers to the reverse strand
		//var hitsRv *map[int][]Info = GetAnchorMatches(r1Read, gtaIndex.SuffixTreeReverseStrandForwardDirection)

		var hitsRv *map[int][]Info = GetAnchorMatches(r1Read, gtaIndex.SuffixTreeForwardStrandForwardDirection)

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
	for i := 0; i < config.KmerLength(); i++ {
		*cigarList = append(*cigarList, 'M')
	}
	return cigarList
}

func finalizeCigar(cigarList *[]rune, startPositionInTranscript uint32, transcript *index.Transcript) (int, string) {

	fmt.Println("finalizing cigar")
	fmt.Println("startPositionInTranscript", startPositionInTranscript)

	// the final cigar string which is built in this function
	cigarString := ""
	// the start position of the match relative to the parent gene
	startRelative := 0

	// the current position within the reference (relative to the parent gene)
	var posInRef uint32 = 0

	lastCigarElement := (*cigarList)[0]
	occurencesCigarElements := 1

	// the index of the current exon
	currentExonId := 0

	// find the exon where the read starts
	for i := 0; i < len(transcript.Exons); i++ {

		// the length of the current exon
		lenExon := transcript.Exons[i].EndRelative - transcript.Exons[i].StartRelative

		if startPositionInTranscript >= lenExon {
			startPositionInTranscript -= lenExon
			continue
		}

		currentExonId = i

		posInRef = transcript.Exons[i].StartRelative + startPositionInTranscript

		break
	}

	startRelative = int(posInRef)

	fmt.Println("currentExonId", currentExonId)
	fmt.Println("posInRef", posInRef)

	for i := 1; i < len(*cigarList); i++ {

		posInRef += 1

		// read match exceeds current exon
		if posInRef > transcript.Exons[currentExonId].EndRelative {

			cigarString += fmt.Sprintf("%d%c", occurencesCigarElements, lastCigarElement)

			lenIntron := transcript.Exons[currentExonId+1].StartRelative - transcript.Exons[currentExonId].EndRelative - 1

			// add intron to cigar
			cigarString += fmt.Sprintf("%dN", lenIntron)
			// reset current cigar element counter
			lastCigarElement = 'N'
			occurencesCigarElements = 0

			// set the position in the reference to the start of the next exon
			posInRef += lenIntron

			// set the current exon to the next exon
			currentExonId += 1
		}

		currentCigarElement := (*cigarList)[i]

		if currentCigarElement == lastCigarElement {

			occurencesCigarElements++

		} else {

			if lastCigarElement != 'N' {
				cigarString += fmt.Sprintf("%d%c", occurencesCigarElements, lastCigarElement)
			}
			occurencesCigarElements = 1
			lastCigarElement = currentCigarElement

		}
	}
	if lastCigarElement != 'N' {
		cigarString += fmt.Sprintf("%d%c", occurencesCigarElements, lastCigarElement)
	}

	return startRelative, cigarString
}

func AlignRead(read *fastq.Read, transcriptId int, kmerHitList []Info, gtaIndex *index.GtaIndex) *sam.Entry {

	fmt.Println(read.Sequence)
	fmt.Println("transcript", gtaIndex.Transcripts[transcriptId].TranscriptIdEnsembl)
	fmt.Println(kmerHitList)

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

				score, cigarPart, seq1, seq2 := algorithms.NeedlemanWunsch(
					gtaIndex.Transcripts[transcriptId].SequenceDnaForwardStrandForwardDirection[startPrefixGapRef:endPrefixGapRef],
					read.Sequence[startPrefixGapRead:endPrefixGapRead])

				fmt.Println(score)
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
		startGapRead := kmerHitList[i-1].IndexRead + config.KmerLength()
		endGapRead := hit.IndexRead
		lenGapRead := endGapRead - startGapRead
		// check gap between current and previous kmer in reference
		startGapRef := kmerHitList[i-1].IndexReference + config.KmerLength()
		endGapRef := hit.IndexReference
		lenGapRef := endGapRef - startGapRef

		if lenGapRead == lenGapRef {
			// TODO: (performance) check if there is any smarter way of comparing strings of same length
		}

		score, cigarPart, seq1, seq2 := algorithms.NeedlemanWunsch(
			gtaIndex.Transcripts[transcriptId].SequenceDnaForwardStrandForwardDirection[startGapRef:endGapRef],
			read.Sequence[startGapRead:endGapRead])

		fmt.Println(score)
		fmt.Println(seq1)
		fmt.Println(seq2)

		// add the alignment to the cigar string
		cigarList = append(cigarList, cigarPart...)

		// add the kmer match to the cigar string
		addKmerMatchToCigar(&cigarList)

		// the last kmer does not end at the end of the read
		if i == len(kmerHitList)-1 && hit.IndexRead < len(read.Sequence)-config.KmerLength() {

			fmt.Println("last kmer does not end at the end of the read")

			startSuffixGapRead := hit.IndexRead + config.KmerLength()
			endSuffixGapRead := len(read.Sequence)
			lenSuffixGapRead := endSuffixGapRead - startSuffixGapRead

			startSuffixGapRef := hit.IndexReference + config.KmerLength()
			endSuffixGapRef := startSuffixGapRef + lenSuffixGapRead

			score, cigarPart, seq1, seq2 = algorithms.NeedlemanWunsch(
				gtaIndex.Transcripts[transcriptId].SequenceDnaForwardStrandForwardDirection[startSuffixGapRef:endSuffixGapRef],
				read.Sequence[startSuffixGapRead:endSuffixGapRead])

			fmt.Println(score)
			fmt.Println(seq1)
			fmt.Println(seq2)

			// add the alignment of the suffix to the cigar string
			cigarList = append(cigarList, cigarPart...)
		}
	}

	fmt.Println("leftmostMappingPosition", leftmostMappingPosition)

	startPos, cigarString := finalizeCigar(&cigarList, uint32(leftmostMappingPosition), gtaIndex.Transcripts[transcriptId])

	return &sam.Entry{
		Qname:        strings.Split(read.Header[1:], " ")[0],
		Cigar:        cigarString,
		Pos:          startPos + 1, // start position is 1-based in SAM file
		TranscriptId: transcriptId,
		Qual:         read.Quality,
	}
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

	// TODO: check if discard read pair or report as mate unmapped
	if hitsR1 == nil {
		logrus.Info("Discard read pair because R1 does not match")
		return ""
	}

	// TODO: implement check for single end reads

	var hitsR2 *map[int][]Info

	if isForwardStrand {

		readPair.ReadR2.Sequence = utils.ReverseComplementDNA(readPair.ReadR2.Sequence)

		hitsR2 = GetAnchorMatches(readPair.ReadR2, gtaIndex.SuffixTreeForwardStrandForwardDirection)

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

		samEntryR1 := AlignRead(readPair.ReadR1, idTranscript, (*hitsR1)[idTranscript], gtaIndex)
		samEntryR2 := AlignRead(readPair.ReadR2, idTranscript, (*hitsR2)[idTranscript], gtaIndex)

		samEntryR1.Rname = gtaIndex.Gene.Chromosome
		samEntryR1.Seq = readPair.ReadR1.Sequence
		samEntryR1.Pnext = samEntryR2.Pos
		samEntryR1.Rnext = "="

		samEntryR2.Rname = gtaIndex.Gene.Chromosome
		samEntryR2.Seq = readPair.ReadR2.Sequence
		samEntryR2.Pnext = samEntryR1.Pos
		samEntryR2.Rnext = "="

		return samEntryR1.String() + "\n" + samEntryR2.String() + "\n"

	} else {

		fmt.Println("not implemented yet")

		/*
			AlignRead(readPair.ReadR1, &gtaIndex.Transcripts[idTranscript].SequenceDnaReverseStrandForwardDirection,
				gtaIndex.SuffixTreeReverseStrandForwardDirection, (*hitsR1)[idTranscript])
		*/
	}

	return "mapping result dummy\n"
}
