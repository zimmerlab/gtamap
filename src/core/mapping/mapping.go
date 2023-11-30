package mapping

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core"
	"github.com/KleinSamuel/gtamap/src/core/algorithms"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
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
					//gtaIndex.Transcripts[transcriptId].SequenceDnaForwardStrandForwardDirection[startPrefixGapRef:endPrefixGapRef],
					gtaIndex.Transcripts[transcriptId].SequenceDnaForward53[startPrefixGapRef:endPrefixGapRef],
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
			//gtaIndex.Transcripts[transcriptId].SequenceDnaForwardStrandForwardDirection[startGapRef:endGapRef],
			gtaIndex.Transcripts[transcriptId].SequenceDnaForward53[startGapRef:endGapRef],
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
				//gtaIndex.Transcripts[transcriptId].SequenceDnaForwardStrandForwardDirection[startSuffixGapRef:endSuffixGapRef],
				gtaIndex.Transcripts[transcriptId].SequenceDnaForward53[startSuffixGapRef:endSuffixGapRef],
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

// MapReadPair map a read pair to the reference
// the main mapping function to be called from outside
//func MapReadPair(readPair *fastq.ReadPair, gtaIndex *index.GtaIndex) string {
//
//	rvReadHeader := ""
//	if readPair.ReadR2 != nil {
//		rvReadHeader = readPair.ReadR2.Header
//	}
//
//	logrus.WithFields(logrus.Fields{
//		"fwReadHeader": readPair.ReadR1.Header,
//		"rvReadHeader": rvReadHeader,
//	}).Info("Map new read pair")
//
//	fmt.Println("R1")
//	fmt.Println(readPair.ReadR1.Sequence)
//	fmt.Println("R2")
//	fmt.Println(readPair.ReadR2.Sequence)
//
//	hitsR1, isForwardStrand := LocateR1PositionOnStrands(gtaIndex, readPair.ReadR1)
//
//	// TODO: check if discard read pair or report as mate unmapped
//	if hitsR1 == nil {
//		logrus.Info("Discard read pair because R1 does not match")
//		return ""
//	}
//
//	// TODO: implement check for single end reads
//
//	var hitsR2 *map[int][]Info
//
//	if isForwardStrand {
//
//		readPair.ReadR2.Sequence = utils.ReverseComplementDNA(readPair.ReadR2.Sequence)
//
//		hitsR2 = GetAnchorMatches(readPair.ReadR2, gtaIndex.SuffixTree)
//
//	} else {
//		hitsR2 = GetAnchorMatches(readPair.ReadR2, gtaIndex.SuffixTree)
//	}
//
//	if hitsR2 == nil {
//		logrus.Info("Discard read pair because R2 does not match")
//		return ""
//	}
//
//	// determine most likely transcript by using the count of kmer hits
//
//	maxHits := 0
//	idTranscript := -1
//
//	for index1, _ := range *hitsR1 {
//		for index2, _ := range *hitsR2 {
//			if index1 == index2 {
//
//				// TODO: filter unique hits per transcript
//
//				// TODO: find a better way to determine the most likely transcript (maybe use fragment length)
//
//				if len((*hitsR1)[index1])+len((*hitsR2)[index2]) > maxHits {
//					maxHits = len((*hitsR1)[index1]) + len((*hitsR2)[index2])
//					idTranscript = index1
//				}
//			}
//		}
//	}
//
//	fmt.Println("idTranscript", idTranscript)
//	fmt.Println("maxHits", maxHits)
//
//	if isForwardStrand {
//
//		samEntryR1 := AlignRead(readPair.ReadR1, idTranscript, (*hitsR1)[idTranscript], gtaIndex)
//		samEntryR2 := AlignRead(readPair.ReadR2, idTranscript, (*hitsR2)[idTranscript], gtaIndex)
//
//		samEntryR1.Rname = gtaIndex.Gene.Chromosome
//		samEntryR1.Seq = readPair.ReadR1.Sequence
//		samEntryR1.Pnext = samEntryR2.Pos
//		samEntryR1.Rnext = "="
//
//		samEntryR2.Rname = gtaIndex.Gene.Chromosome
//		samEntryR2.Seq = readPair.ReadR2.Sequence
//		samEntryR2.Pnext = samEntryR1.Pos
//		samEntryR2.Rnext = "="
//
//		return samEntryR1.String() + "\n" + samEntryR2.String() + "\n"
//
//	} else {
//
//		fmt.Println("not implemented yet")
//
//		/*
//			AlignRead(readPair.ReadR1, &gtaIndex.Transcripts[idTranscript].SequenceDnaReverseStrandForwardDirection,
//				gtaIndex.SuffixTreeReverseStrandForwardDirection, (*hitsR1)[idTranscript])
//		*/
//	}
//
//	return "mapping result dummy\n"
//}

func MapReadPairDev(readPair *fastq.ReadPair, gtaIndex *index.GtaIndex) string {

	fmt.Println("index num sequences", gtaIndex.NumSequences)

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

	errorRate := 0.05

	// map R1

	// the maximum number of mismatches allowed
	maxMismatches := int(float64(len(readPair.ReadR1.Sequence)) * errorRate)
	// the min length of each kmer
	lenKmer := 8
	// the number of kmers created from the read
	numKmer := len(readPair.ReadR1.Sequence) / lenKmer
	// the number of bases left over after creating the kmers with equal length
	rest := len(readPair.ReadR1.Sequence) % numKmer

	// the current position within the read
	position := 0
	// the state of the discard step
	failed := false
	// the number of mismatches per sequence
	mismatches := make(map[int]core.DiscardStepMatchInformation, gtaIndex.NumSequences)
	for i := 0; i < gtaIndex.NumSequences; i++ {
		mismatches[i] = core.DiscardStepMatchInformation{
			NumMismatches: 0,
			Matches:       make([]core.ExactMatch, 0),
		}
	}

	// TODO: idea: randomize the kmer that is used for filtering

	for position < len(readPair.ReadR1.Sequence) {

		currentLenKmer := lenKmer
		if rest > 0 {
			currentLenKmer++
			rest--
		}

		// the actual kmer
		kmer := readPair.ReadR1.Sequence[position : position+currentLenKmer]

		logrus.WithFields(logrus.Fields{
			"position": position,
			"length":   currentLenKmer,
			"kmer":     kmer,
		}).Debug("create new kmer")

		// the result of the exact matching against the suffix tree
		result := gtaIndex.SuffixTree.Search(&kmer)

		// contains true for each sequence index that has a match
		matches := make([]bool, gtaIndex.NumSequences)

		// there was at least one exact match in any sequence
		if result != nil {

			for _, match := range result.Matches {

				logrus.WithFields(logrus.Fields{
					"sequenceIndex": match.SequenceIndex,
					"fromTarget":    match.FromTarget,
					"toTarget":      match.ToTarget,
				}).Debug("match")

				matches[match.SequenceIndex] = true
			}
		}

		// iterate every sequence index to add to their mismatch count
		for sequenceIndex := 0; sequenceIndex < len(matches); sequenceIndex++ {

			// skip if there was a match in this sequence
			if matches[sequenceIndex] {
				continue
			}

			logrus.WithFields(logrus.Fields{
				"sequenceIndex": sequenceIndex,
			}).Debug("add mismatch")

			// sequence index was not already discarded
			if val, ok := mismatches[sequenceIndex]; ok {

				logrus.WithFields(logrus.Fields{
					"sequenceIndex": sequenceIndex,
					"numMismatches": val.NumMismatches,
					"numMatches":    len(val.Matches),
				}).Debug("add mismatch to sequence index")

				val.NumMismatches += 1
				mismatches[sequenceIndex] = val

				// remove the sequence index if it has too many mismatches
				if val.NumMismatches > maxMismatches {
					delete(mismatches, sequenceIndex)

					logrus.Debug("discard sequence index because too many mismatches")
				}

			} else {
				logrus.Debug("sequence index already discarded")
			}
		}

		// stop if there are no more sequences left and the read can be discarded
		if len(mismatches) == 0 {
			logrus.Debug("no more sequences left")

			failed = true
			break
		}

		position += currentLenKmer
	}

	if failed {
		logrus.Info("discard this read!")
		return "discard"
	}

	logrus.Info("keep this read!")

	for val := range mismatches {
		logrus.Info(val)
	}

	return "end"
}
