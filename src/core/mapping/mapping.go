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
	"sort"
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

func MapReadPairDev(readPair *fastq.ReadPair, index *index.GtaIndex) string {

	logrus.Info("map new read pair")

	logrus.WithFields(logrus.Fields{
		"header": readPair.ReadR1.Header,
	}).Info("got read R1")
	logrus.Info(readPair.ReadR1.Sequence)

	if readPair.ReadR2 != nil {
		logrus.WithFields(logrus.Fields{
			"header": readPair.ReadR2.Header,
		}).Info("got read R2")
		logrus.Info(readPair.ReadR2.Sequence)
	} else {
		logrus.Info("no read R2")
	}

	logrus.WithFields(logrus.Fields{
		"header": readPair.ReadR1.Header,
	}).Info("map R1 read to reference")

	resultString := ""

	resultR1, valid1 := MapRead(readPair.ReadR1, index)

	fmt.Println("result R1")
	fmt.Println(valid1)
	fmt.Println(resultR1)

	if valid1 {
		resultString += readPair.ReadR1.Header + "\tkeep"
	} else {
		resultString += readPair.ReadR1.Header + "\tdiscard"
	}
	resultString += "\n"

	logrus.WithFields(logrus.Fields{
		"header": readPair.ReadR2.Header,
	}).Info("map R2 read to reference")

	resultR2, valid2 := MapRead(readPair.ReadR2, index)

	fmt.Println("result R2")
	fmt.Println(valid2)
	fmt.Println(resultR2)

	if valid2 {
		resultString += readPair.ReadR2.Header + "\tkeep"
	} else {
		resultString += readPair.ReadR2.Header + "\tdiscard"
	}
	resultString += "\n"

	return resultString
}

func MapRead(read *fastq.Read, index *index.GtaIndex) (core.ReadMapResult, bool) {

	// TODO: move this to config
	errorRate := 0.05

	// the final result of this function containing the mapping information of each context for each sequence
	result := core.ReadMapResult{
		SequenceMatches: make(map[int]core.SequenceMapResult),
	}

	// the maximum number of mismatches allowed
	maxMismatches := int(float64(len(read.Sequence)) * errorRate)
	// the min length of each kmer
	lenKmer := 8
	// the number of kmers created from the read
	numKmer := len(read.Sequence) / lenKmer
	// the number of bases left over after creating the kmers with equal length
	rest := len(read.Sequence) % numKmer

	// the current position within the read
	position := 0
	// the state of the discard step
	failed := false
	// the number of mismatches per sequence (sequenceIndex = sequence index)
	matchMismatchCounts := make(map[int]core.DiscardStepMatchInformation, index.NumSequences)
	for i := 0; i < index.NumSequences; i++ {
		matchMismatchCounts[i] = core.DiscardStepMatchInformation{
			NumMismatches: 0,
			Matches:       make([]core.SequenceMatch, 0),
		}
	}

	// TODO: idea: randomize the kmer that is used for filtering
	// TODO: (maybe even overlapping but this requires more postprocessing during extend step)

	for position < len(read.Sequence) {

		currentLenKmer := lenKmer
		if rest > 0 {
			currentLenKmer++
			rest--
		}

		// the actual kmer
		kmer := read.Sequence[position : position+currentLenKmer]

		logrus.WithFields(logrus.Fields{
			"position": position,
			"length":   currentLenKmer,
			"kmer":     kmer,
		}).Debug("create new kmer")

		// the result of the exact matching of this kmer against the suffix tree
		mappingResult := index.SuffixTree.FindPatternExact(&kmer)

		if mappingResult != nil {
			for _, match := range mappingResult.Matches {
				fmt.Println(match)
				fmt.Println(len(index.GetTranscriptSequenceDna(match.SequenceIndex)))
			}
		}

		// contains true for each sequence index that has a match
		matches := make([]bool, index.NumSequences)

		// there was at least one exact match in any sequence
		if mappingResult != nil {

			for _, match := range mappingResult.Matches {

				// skip all sequence index matches of sequences that were already discarded
				if _, ok := matchMismatchCounts[match.SequenceIndex]; !ok {
					continue
				}

				logrus.WithFields(logrus.Fields{
					"sequenceIndex": match.SequenceIndex,
					"fromTarget":    match.FromTarget,
					"toTarget":      match.ToTarget,
				}).Debug("match")

				matches[match.SequenceIndex] = true

				currentCounts := matchMismatchCounts[match.SequenceIndex]
				match.FromSource = position
				match.ToSource = position + currentLenKmer
				currentCounts.Matches = append(currentCounts.Matches, match)
				matchMismatchCounts[match.SequenceIndex] = currentCounts
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
			if val, ok := matchMismatchCounts[sequenceIndex]; ok {

				logrus.WithFields(logrus.Fields{
					"sequenceIndex": sequenceIndex,
					"numMismatches": val.NumMismatches,
					"numMatches":    len(val.Matches),
				}).Debug("add mismatch to sequence index")

				val.NumMismatches += 1
				matchMismatchCounts[sequenceIndex] = val

				// remove the sequence index if it has too many mismatches
				if val.NumMismatches > maxMismatches {
					delete(matchMismatchCounts, sequenceIndex)

					logrus.Debug("discard sequence index because too many mismatches")
				}

			} else {
				logrus.Debug("sequence index already discarded")
			}
		}

		// stop if there are no more sequences left and the read can be discarded
		if len(matchMismatchCounts) == 0 {
			logrus.Debug("no more sequences left")

			failed = true
			break
		}

		position += currentLenKmer
	}

	if failed {
		logrus.Info("discard this read!")
		return result, false
	}

	logrus.Info("keep this read!")

	// iterate the mappings (different contexts) on each sequence (identified by sequence index)
	for sequenceIndex, val := range matchMismatchCounts {

		logrus.WithFields(logrus.Fields{
			"sequenceIndex": sequenceIndex,
			"numMatches":    len(val.Matches),
		}).Info("process non-discarded sequence index")

		// the final result for this sequence index
		sequenceMapResult := core.SequenceMapResult{
			SequenceContextMatches: make(map[int]core.SequenceContextMatch),
		}

		matches := val.Matches
		sort.Sort(matches)

		// assert that the matches list is sorted
		if config.Env() == "development" {
			logrus.Debug("assert that the matches list is sorted")
			currentIndex := 0
			for i := 0; i < len(matches); i++ {
				if matches[i].FromSource < currentIndex {
					panic("matches are not sorted")
				}
				currentIndex = matches[i].FromSource
			}
			logrus.Debug("assertion passed")
		}

		contexts := make(map[int]core.SequenceContextMatch)

		// sort the matches to the contexts on each sequence index and combine consecutive matches
		for i, match := range matches {

			// the start position of the entire read within the reference sequence
			targetIndex := match.FromTarget - match.FromSource

			logrus.WithFields(logrus.Fields{
				"matchIndex":    i,
				"sequenceIndex": sequenceIndex,
				"targetIndex":   targetIndex,
				"fromSource":    match.FromSource,
				"toSource":      match.ToSource,
				"fromTarget":    match.FromTarget,
				"toTarget":      match.ToTarget,
			}).Debug("process match")

			lengthTranscript := len(index.GetTranscriptSequenceDna(sequenceIndex))
			lengthRead := len(read.Sequence)

			// discard this match because the read is too long for the transcript for this target index
			if targetIndex > lengthTranscript-lengthRead {
				logrus.WithFields(logrus.Fields{
					"lengthTranscript": lengthTranscript,
					"lengthRead":       lengthRead,
					"targetIndex":      targetIndex,
				}).Info("discard match because target index is out of bounds")
				continue
			}

			context := core.SequenceContextMatch{
				SequenceIndex: sequenceIndex,
				TargetIndex:   -1,
				Matches:       make([]core.SequenceMatch, 0),
			}

			if val, ok := contexts[targetIndex]; ok {
				context = val
			}

			// no match yet in this context -> add the current match
			if context.Matches.Len() == 0 {
				context.TargetIndex = targetIndex
				context.Matches = append(context.Matches, match)
				contexts[targetIndex] = context

				logrus.Debug("no match yet in this context -> add the current match")
				continue
			}

			// the last match to be checked for extension
			lastMatch := context.Matches[len(context.Matches)-1]

			// check if the current match extends the last match in this context and combine them
			if match.FromSource == lastMatch.ToSource {
				logrus.Debug("match extends last match")

				lastMatch.ToSource = match.ToSource
				lastMatch.ToTarget = match.ToTarget
				lastMatch.Mismatches = append(lastMatch.Mismatches, match.Mismatches...)

				context.Matches[len(context.Matches)-1] = lastMatch
			} else {
				// the current match does not extend the last match
				logrus.Debug("match does not extend last match")

				// add the current match if not overlapping
				context.Matches = append(context.Matches, match)
			}
		}

		logrus.WithFields(logrus.Fields{
			"sequenceIndex": sequenceIndex,
			"numContexts":   len(contexts),
		}).Info("contexts of sequence index determined")

		for _, context := range contexts {
			logrus.Info(context)
		}

		// extend the matches in the context
		//contextloop:
		for _, context := range contexts {

			logrus.WithFields(logrus.Fields{
				"sequenceIndex": context.SequenceIndex,
				"targetIndex":   context.TargetIndex,
			}).Info("map the unmapped regions of this context")

			logrus.Info("matches of this context:")
			for _, match := range context.Matches {
				logrus.Info(match)
			}

			// containing the whole context with exact matches and mapped unmapped regions
			wholeContext := core.SequenceContextMatch{
				SequenceIndex: context.SequenceIndex,
				TargetIndex:   context.TargetIndex,
				Matches:       make([]core.SequenceMatch, 0),
			}

			// the sequence of the read
			sourceSequence := read.Sequence
			// the starting position of the match in the target sequence (transcript)
			transcriptIndex := context.SequenceIndex / 4
			// the index of the direction of the target sequence (FF, FR, RF, RR)
			directionIndex := context.SequenceIndex % 4
			// the target sequence (transcript)
			targetSequence := index.Transcripts[transcriptIndex].GetSequenceDna(directionIndex)

			// map the unmapped regions and combine with exact match regions

			// the number of mismatches in this context
			numMismatches := 0

			// determine the unmapped regions between the matches and map them
			for i, match := range context.Matches {

				if i == 0 {
					// there is an unmapped region before the first exact match
					if match.FromSource > 0 {

						currentMatch, valid := MapBaseByBase(
							0,
							match.FromSource,
							context.TargetIndex,
							context.TargetIndex+match.FromSource,
							sourceSequence,
							targetSequence,
							&numMismatches,
							maxMismatches)

						if !valid {
							logrus.Debug("too many mismatches!")
							break
						}

						currentMatch.SequenceIndex = match.SequenceIndex

						wholeContext.Matches = append(wholeContext.Matches, currentMatch)

					}

				} else {

					prevMatch := context.Matches[i-1]

					// map the unmapped region between the matches
					currentMatch, valid := MapBaseByBase(
						prevMatch.ToSource,
						match.FromSource,
						prevMatch.ToTarget,
						match.FromTarget,
						sourceSequence,
						targetSequence,
						&numMismatches,
						maxMismatches)

					if !valid {
						logrus.Debug("too many mismatches!")
						break
					}

					currentMatch.SequenceIndex = match.SequenceIndex

					wholeContext.Matches = append(wholeContext.Matches, currentMatch)
				}

				// add the exact sequence match
				wholeContext.Matches = append(wholeContext.Matches, match)

				// the last match in this read
				if i == len(context.Matches)-1 {

					// there is an unmapped region at the end of the read
					if match.ToSource < len(read.Sequence) {

						// map the unmapped region between the last
						currentMatch, valid := MapBaseByBase(
							match.ToSource,
							len(read.Sequence),
							match.ToTarget,
							match.ToTarget+len(read.Sequence),
							sourceSequence,
							targetSequence,
							&numMismatches,
							maxMismatches)

						if !valid {
							logrus.Debug("too many mismatches!")
							break
						}

						currentMatch.SequenceIndex = match.SequenceIndex

						wholeContext.Matches = append(wholeContext.Matches, currentMatch)
					}
				}
			}

			// discard this context if there are more mismatches than allowed
			if numMismatches > maxMismatches {
				logrus.WithFields(logrus.Fields{
					"sequenceIndex": context.SequenceIndex,
					"targetIndex":   context.TargetIndex,
					"numMismatches": int(numMismatches),
					"maxMismatches": maxMismatches,
				}).Info("discard this context because of too many mismatches")
				continue
			}

			logrus.WithFields(logrus.Fields{}).Info("add the whole context to the result")
			logrus.Info(wholeContext)

			sequenceMapResult.SequenceContextMatches[wholeContext.TargetIndex] = wholeContext
		}

		result.SequenceMatches[sequenceIndex] = sequenceMapResult
	}

	return result, true
}

func MapBaseByBase(fromSource int, toSource int, fromTarget int, toTarget int,
	sourceSequence string, targetSequence string, numMismatches *int, maxMismatches int) (core.SequenceMatch, bool) {

	sequenceMatch := core.SequenceMatch{
		SequenceIndex: -1,
		FromSource:    fromSource,
		ToSource:      toSource,
		FromTarget:    fromTarget,
		ToTarget:      toTarget,
		Mismatches:    make([]int, 0),
	}

	// make base by base comparison and return if more mismatches than threshold
	for i := fromSource; i < toSource; i++ {

		sourceBase := sourceSequence[i]
		targetBase := targetSequence[fromTarget+i]

		if sourceBase != targetBase {
			*numMismatches++

			logrus.WithFields(logrus.Fields{
				"sourceIndex":   i,
				"targetIndex":   fromTarget + i,
				"sourceBase":    string(sourceBase),
				"targetBase":    string(targetBase),
				"numMismatches": numMismatches,
			}).Debug("mismatch")

			if *numMismatches > maxMismatches {
				logrus.Debug("too many mismatches -> return false")
				return core.SequenceMatch{}, false
			}

			sequenceMatch.Mismatches = append(sequenceMatch.Mismatches, i)
		}
	}

	return sequenceMatch, true
}
