package mapping

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/interval"
	"github.com/KleinSamuel/gtamap/src/core/timer"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
	"github.com/sirupsen/logrus"
	"math/rand"
	"strconv"
	"time"
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

//func addKmerMatchToCigar(cigarList *[]rune) *[]rune {
//	for i := 0; i < config.KmerLength(); i++ {
//		*cigarList = append(*cigarList, 'M')
//	}
//	return cigarList
//}

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

// MapReadPairDev is the core method that maps a read pair to the reference using the index
// Each read (R1 and R2) is mapped to the reference sequences using the index.
// The mapping is done in a sublinear time complexity by discarding reads that have too many mismatches.
// After mapping both reads, a list of possible locations is returned which contains the transcript, the position
// and the locations of the mismatches.
// The results of R1 and R2 are combined to determine the most likely transcript.
func MapReadPairDev(readPair *fastq.ReadPair, index *index.GtaIndex, timerChannel chan<- *timer.Timer) string {

	//logrus.Info("map new read pair")

	logrus.WithFields(logrus.Fields{
		"header": readPair.ReadR1.Header,
	}).Debug("got read R1")
	// TODO: remove after debugging
	//logrus.Info(readPair.ReadR1.Sequence)

	if readPair.ReadR2 != nil {
		logrus.WithFields(logrus.Fields{
			"header": readPair.ReadR2.Header,
		}).Debug("got read R2")
		// TODO: remove after debugging
		//logrus.Info(readPair.ReadR2.Sequence)
	} else {
		//logrus.Info("no read R2")
	}

	//logrus.WithFields(logrus.Fields{
	//	"header": readPair.ReadR1.Header,
	//}).Debug("map R1 read to reference")

	result := &core.ReadMappingPreResult{
		ReadR1:    readPair.ReadR1,
		ReadR2:    readPair.ReadR2,
		ResultsR1: nil,
		ResultsR2: nil,
	}

	resultString := ""

	t := timer.NewTimer()

	timerStart := time.Now()

	resultsR1, wasMappedR1 := MapRead(readPair.ReadR1, index, t)

	t.MapR1 = time.Since(timerStart)

	if wasMappedR1 {
		result.ResultsR1 = &resultsR1
	}

	if readPair.ReadR2 != nil && (wasMappedR1 || config.IncludeReadsImproperlyPaired()) {

		timerStart = time.Now()

		resultsR2, wasMappedR2 := MapRead(readPair.ReadR2, index, t)

		t.MapR2 = time.Since(timerStart)

		if wasMappedR2 {
			result.ResultsR2 = &resultsR2
		}
	}

	timerStart = time.Now()

	recordPairs := DetermineReadLocation(result, index)

	t.DetermineReadLocation = time.Since(timerStart)

	if recordPairs != nil {
		for _, recordPair := range *recordPairs {
			resultString += recordPair.First.String() + "\n" + recordPair.Second.String() + "\n"
		}
	}

	timerChannel <- t

	return resultString
}

func DetermineReadLocation(result *core.ReadMappingPreResult, index *index.GtaIndex) *[]*sam.RecordPair {

	logrus.Debug("determine read location")

	if result.ResultsR1 == nil && result.ResultsR2 == nil {
		logrus.Debug("both reads were not mapped")

		return nil

	} else if result.ResultsR1 != nil && result.ResultsR2 != nil {
		logrus.Debug("both reads have been mapped")

		return DetermineReadLocationPaired(result, index)

	} else {
		logrus.Debug("only one read has been mapped")

		return DetermineReadLocationUnpaired(result, index)
	}
}

type ReadLocationsOnReference struct {
	TranscriptId     int
	LocationsForward []*core.InexactMatchResult
	LocationsReverse []*core.InexactMatchResult
}

func DetermineReadLocationPaired(result *core.ReadMappingPreResult, index *index.GtaIndex) *[]*sam.RecordPair {

	logrus.Debug("determine read location paired")

	properPairs := make([]*core.ProperPairCandidate, 0)

	// find proper pairs
	for sequenceIndexR1, resultsR1 := range *result.ResultsR1 {

		transcriptIndex := index.SequenceIndexToTranscriptIndex(uint32(sequenceIndexR1))
		isForwardStrandR1 := index.SequenceIndexIsForward(uint32(sequenceIndexR1))

		matchingSequenceIndex := int(transcriptIndex) * 2
		if isForwardStrandR1 {
			matchingSequenceIndex += 1
		}

		// no candidate for a proper pair found (no rev comp on same reference)
		if _, ok := (*result.ResultsR2)[matchingSequenceIndex]; !ok {
			continue
		}

		resultsR2 := (*result.ResultsR2)[matchingSequenceIndex]

		for _, resultR1 := range resultsR1 {

			// The position of the read on the reference to compute the fragment length.
			// Start of forward or end of reverse with respect to forward strand such that the two most distant
			// positions are considered.
			positionOnReferenceR1 := resultR1.FromTarget
			if !isForwardStrandR1 {
				positionOnReferenceR1 = index.Transcripts[transcriptIndex].SequenceLength - resultR1.FromTarget
			}

			for _, resultR2 := range resultsR2 {

				positionOnReferenceR2 := resultR2.FromTarget
				if isForwardStrandR1 {
					positionOnReferenceR2 = index.Transcripts[transcriptIndex].SequenceLength - resultR2.FromTarget
				}

				fragmentLength := positionOnReferenceR2 - positionOnReferenceR1
				if fragmentLength < 0 {
					fragmentLength = -fragmentLength
				}

				//fmt.Println(fragmentLength)

				// TODO: use fragment length to filter proper pairs

				// TODO: for now every pair on the same reference but different strands is considered a proper pair

				properPairs = append(properPairs, &core.ProperPairCandidate{
					ReferenceIndex: int(transcriptIndex),
					FragmentLength: fragmentLength,
					R1isForward:    isForwardStrandR1,
					ResultR1:       resultR1,
					ResultR2:       resultR2,
				})
			}
		}
	}

	if !config.IncludeReadsImproperlyPaired() && len(properPairs) == 0 {
		logrus.Debug("no proper pair found")
		return nil
	}

	recordPairs := make([]*sam.RecordPair, 0)

	genomicPositions := make(map[uint32]bool)

	for _, properPair := range properPairs {
		// the genomic position of the read on the reference
		if properPair.R1isForward {
			genomicPositions[index.TranslateRelativeTranscriptPositionToRelativeGenePosition(uint32(properPair.ReferenceIndex), uint32(properPair.ResultR1.FromTarget))] = true
		} else {
			genomicPositions[index.TranslateRelativeTranscriptPositionToRelativeGenePosition(uint32(properPair.ReferenceIndex), index.TranslateRwTransPosToFwTransPos(uint32(properPair.ReferenceIndex), uint32(properPair.ResultR2.FromTarget)))] = true
		}
	}

	if len(genomicPositions) > 1 && !config.IncludeReadsAmbiguouslyMapped() {

		logrus.WithFields(logrus.Fields{
			"genomicPositions": genomicPositions,
		}).Debug("read is ambiguously mapped")

		for _, properPair := range properPairs {

			firstFlag := sam.Flag{}
			firstFlag.SetPaired()
			firstFlag.SetProperlyPaired()
			firstFlag.SetFirstInPair()

			secondFlag := sam.Flag{}
			secondFlag.SetPaired()
			secondFlag.SetProperlyPaired()
			secondFlag.SetSecondInPair()

			if properPair.R1isForward {
				firstFlag.SetMateReverseStrand()
				secondFlag.SetReverseStrand()
			} else {
				firstFlag.SetReverseStrand()
				secondFlag.SetMateReverseStrand()
			}

			recordPairs = append(recordPairs, &sam.RecordPair{
				First: &sam.Record{
					Qname: result.ReadR1.Header,
					Flag:  firstFlag,
					Pos:   int(index.TranslateRelativeTranscriptPositionToRelativeGenePosition(uint32(properPair.ReferenceIndex), uint32(properPair.ResultR1.FromTarget))),
					Rname: index.Gene.Chromosome,
					Mapq:  255,
					Cigar: strconv.Itoa(len(result.ReadR1.Sequence)) + "M",
					Rnext: "=",
					Pnext: 0,
					Tlen:  properPair.FragmentLength,
					Seq:   result.ReadR1.Sequence,
					Qual:  result.ReadR1.Quality,
				},
				Second: &sam.Record{
					Qname: result.ReadR2.Header,
					Flag:  secondFlag,
					Pos:   int(index.TranslateRelativeTranscriptPositionToRelativeGenePosition(uint32(properPair.ReferenceIndex), uint32(properPair.ResultR2.FromTarget))),
					Rname: index.Gene.Chromosome,
					Mapq:  255,
					Cigar: strconv.Itoa(len(result.ReadR2.Sequence)) + "M",
					Rnext: "=",
					Pnext: 0,
					Tlen:  properPair.FragmentLength,
					Seq:   result.ReadR2.Sequence,
					Qual:  result.ReadR2.Quality,
				},
			})

		}

		return &recordPairs
	}

	firstFlag := sam.Flag{}
	firstFlag.SetPaired()
	firstFlag.SetProperlyPaired()
	firstFlag.SetFirstInPair()

	secondFlag := sam.Flag{}
	secondFlag.SetPaired()
	secondFlag.SetProperlyPaired()
	secondFlag.SetSecondInPair()

	if properPairs[0].R1isForward {
		firstFlag.SetMateReverseStrand()
		secondFlag.SetReverseStrand()
	} else {
		firstFlag.SetReverseStrand()
		secondFlag.SetMateReverseStrand()
	}

	transcriptIndicesR1 := make([]int, len(properPairs))
	transcriptIndicesR2 := make([]int, len(properPairs))
	for i, properPair := range properPairs {
		transcriptIndicesR1[i] = int(index.SequenceIndexToTranscriptIndex(uint32(properPair.ResultR1.SequenceIndex)))
		transcriptIndicesR2[i] = int(index.SequenceIndexToTranscriptIndex(uint32(properPair.ResultR2.SequenceIndex)))
	}

	recordPairs = append(recordPairs, &sam.RecordPair{
		First: &sam.Record{
			Qname:         result.ReadR1.Header,
			Flag:          firstFlag,
			Pos:           int(index.TranslateRelativeTranscriptPositionToRelativeGenePosition(uint32(properPairs[0].ReferenceIndex), uint32(properPairs[0].ResultR1.FromTarget))),
			Rname:         index.Gene.Chromosome,
			Mapq:          255,
			Cigar:         strconv.Itoa(len(result.ReadR1.Sequence)) + "M",
			Rnext:         "=",
			Pnext:         0,
			Tlen:          properPairs[0].FragmentLength,
			Seq:           result.ReadR1.Sequence,
			Qual:          result.ReadR1.Quality,
			TranscriptIds: transcriptIndicesR1,
		},
		Second: &sam.Record{
			Qname:         result.ReadR2.Header,
			Flag:          secondFlag,
			Pos:           int(index.TranslateRelativeTranscriptPositionToRelativeGenePosition(uint32(properPairs[0].ReferenceIndex), uint32(properPairs[0].ResultR2.FromTarget))),
			Rname:         index.Gene.Chromosome,
			Mapq:          255,
			Cigar:         strconv.Itoa(len(result.ReadR2.Sequence)) + "M",
			Rnext:         "=",
			Pnext:         0,
			Tlen:          properPairs[0].FragmentLength,
			Seq:           result.ReadR2.Sequence,
			Qual:          result.ReadR2.Quality,
			TranscriptIds: transcriptIndicesR2,
		},
	})

	return &recordPairs
}

func DetermineReadLocationUnpaired(result *core.ReadMappingPreResult, index *index.GtaIndex) *[]*sam.RecordPair {

	// TODO: change to debugging
	logrus.Debug("determine read location unpaired")

	matches := make([]*core.InexactMatchResult, 0)

	if result.ResultsR1 != nil {
		for _, results := range *result.ResultsR1 {
			matches = append(matches, results...)
		}
	} else {
		for _, results := range *result.ResultsR2 {
			matches = append(matches, results...)
		}
	}

	genePositions := make([]uint32, 0)
	matchesPerGenePosition := make(map[uint32][]*core.InexactMatchResult)

	for _, match := range matches {

		if _, ok := matchesPerGenePosition[match.FromGene]; !ok {
			genePositions = append(genePositions, match.FromGene)
			matchesPerGenePosition[match.FromGene] = make([]*core.InexactMatchResult, 0)
			matchesPerGenePosition[match.FromGene] = append(matchesPerGenePosition[match.FromGene], match)
		} else {
			matchesPerGenePosition[match.FromGene] = append(matchesPerGenePosition[match.FromGene], match)
		}
	}

	logrus.Debug("check if read is uniquely mapped")

	isAmbiguouslyMapped := false

	if len(genePositions) > 1 {
		// the read is ambiguously mapped if there are multiple gene positions

		logrus.WithFields(logrus.Fields{
			"genePositions": genePositions,
		}).Debug("read is ambiguously mapped because of multiple gene positions")

		isAmbiguouslyMapped = true
	} else {
		// The read is mapped to a single gene position, but it can still be
		// ambiguously mapped if the equivalent classes are different.

		// use the first set of equivalence classes as the reference
		ecIds := matchesPerGenePosition[genePositions[0]][0].EquivalenceClassIds

		// compare every other set of equivalence classes with the reference
		for _, match := range matchesPerGenePosition[genePositions[0]][1:] {
			// the classes do not match
			if !index.EquivalenceClassIdsMatch(ecIds, match.EquivalenceClassIds) {

				logrus.WithFields(logrus.Fields{
					"other": ecIds,
					"this":  match.EquivalenceClassIds,
				}).Debug("read is ambiguously mapped because of different equivalence classes")

				isAmbiguouslyMapped = true
				break
			}
		}
	}

	recordPairs := make([]*sam.RecordPair, 0)

	if isAmbiguouslyMapped {

		logrus.Debug("read is ambiguously mapped")

		if config.IncludeReadsAmbiguouslyMapped() {

			logrus.Debug("include ambiguously mapped reads")

			// TODO: return all
			// TODO: determine main mapping location (currently the first one)
			isMainMappingLocation := true

			for _, genePosition := range genePositions {

				transcriptIds := make([]int, len(matchesPerGenePosition[genePosition]))
				for i, match := range matchesPerGenePosition[genePosition] {
					transcriptIds[i] = int(index.SequenceIndexToTranscriptIndex(uint32(match.SequenceIndex)))
				}

				if result.ReadR2 != nil {
					// the reads are paired

					logrus.Debug("reads are paired")

					if result.ResultsR1 != nil {
						// only R1 was mapped

						logrus.Debug("only R1 was mapped")

						flagMapped := sam.Flag{}
						flagMapped.SetPaired()
						flagMapped.SetMateUnmapped()
						flagMapped.SetFirstInPair()

						if !isMainMappingLocation {
							flagMapped.SetNotPrimaryAlignment()
						}

						isForward := index.SequenceIndexIsForward(uint32(matchesPerGenePosition[genePosition][0].SequenceIndex))
						if !isForward {
							flagMapped.SetReverseStrand()
						}

						cigarString := strconv.Itoa(len(result.ReadR1.Sequence)) + "M"

						mappedRecord := &sam.Record{
							Qname:         result.ReadR1.Header,
							Flag:          flagMapped,
							Rname:         index.Gene.Chromosome,
							Pos:           int(genePosition),
							Mapq:          255, // TODO: figure this out
							Cigar:         cigarString,
							Rnext:         result.ReadR2.Header,
							Pnext:         0,
							Tlen:          0, // set to 0 because only one read was mapped
							Seq:           result.ReadR1.Sequence,
							Qual:          result.ReadR1.Quality,
							TranscriptIds: transcriptIds,
						}

						flagUnmapped := sam.Flag{}
						flagUnmapped.SetPaired()
						flagUnmapped.SetUnmapped()
						flagUnmapped.SetSecondInPair()

						unmappedRecord := &sam.Record{
							Qname:         result.ReadR2.Header,
							Flag:          flagUnmapped,
							Rname:         "*",
							Pos:           0,
							Mapq:          255,
							Cigar:         "*",
							Rnext:         result.ReadR1.Header,
							Pnext:         0,
							Tlen:          0,
							Seq:           result.ReadR2.Sequence,
							Qual:          result.ReadR2.Quality,
							TranscriptIds: nil,
						}

						recordPairs = append(recordPairs, &sam.RecordPair{
							First:  mappedRecord,
							Second: unmappedRecord,
						})

					} else {
						// only R2 was mapped

						logrus.Debug("only R2 was mapped")

						flagMapped := sam.Flag{}
						flagMapped.SetPaired()
						flagMapped.SetMateUnmapped()
						flagMapped.SetSecondInPair()

						if !isMainMappingLocation {
							flagMapped.SetNotPrimaryAlignment()
						}

						isForward := index.SequenceIndexIsForward(uint32(matchesPerGenePosition[genePosition][0].SequenceIndex))
						if !isForward {
							flagMapped.SetReverseStrand()
						}

						cigarString := strconv.Itoa(len(result.ReadR1.Sequence)) + "M"

						mappedRecord := &sam.Record{
							Qname:         result.ReadR2.Header,
							Flag:          flagMapped,
							Rname:         index.Gene.Chromosome,
							Pos:           int(genePosition),
							Mapq:          255,
							Cigar:         cigarString,
							Rnext:         result.ReadR1.Header,
							Pnext:         0,
							Tlen:          0, // set to 0 because only one read was mapped
							Seq:           result.ReadR2.Sequence,
							Qual:          result.ReadR2.Quality,
							TranscriptIds: transcriptIds,
						}

						flagUnmapped := sam.Flag{}
						flagUnmapped.SetPaired()
						flagUnmapped.SetUnmapped()
						flagUnmapped.SetFirstInPair()

						unmappedRecord := &sam.Record{
							Qname:         result.ReadR1.Header,
							Flag:          flagUnmapped,
							Rname:         "*",
							Pos:           0,
							Mapq:          255,
							Cigar:         "*",
							Rnext:         result.ReadR2.Header,
							Pnext:         0,
							Tlen:          0,
							Seq:           result.ReadR1.Sequence,
							Qual:          result.ReadR1.Quality,
							TranscriptIds: nil,
						}

						recordPairs = append(recordPairs, &sam.RecordPair{
							First:  unmappedRecord,
							Second: mappedRecord,
						})
					}
				} else {
					// the read was not paired so only R1 is mapped

					logrus.Debug("read was not paired so only R1 was mapped")

					flagMapped := sam.Flag{}
					isForward := index.SequenceIndexIsForward(uint32(matchesPerGenePosition[genePosition][0].SequenceIndex))
					if !isForward {
						flagMapped.SetReverseStrand()
					}

					cigarString := strconv.Itoa(len(result.ReadR1.Sequence)) + "M"

					mappedRecord := &sam.Record{
						Qname:         result.ReadR1.Header,
						Flag:          flagMapped,
						Rname:         index.Gene.Chromosome,
						Pos:           int(genePosition),
						Mapq:          255, // TODO: figure this out
						Cigar:         cigarString,
						Rnext:         "*",
						Pnext:         0,
						Tlen:          0, // set to 0 because only one read was mapped
						Seq:           result.ReadR1.Sequence,
						Qual:          result.ReadR1.Quality,
						TranscriptIds: transcriptIds,
					}

					recordPairs = append(recordPairs, &sam.RecordPair{
						First:  mappedRecord,
						Second: nil,
					})
				}

				isMainMappingLocation = false
			}

		} else {
			return nil
		}

	} else {
		// The default case where the read was uniquely mapped to a single gene position
		// and all potential transcript positions span across the same equivalence classes.

		logrus.Debug("read was uniquely mapped")

		if result.ReadR2 != nil {
			// the reads are paired

			logrus.Debug("reads are paired")

			if result.ResultsR1 != nil {
				// only R1 was mapped

				logrus.Debug("only R1 was mapped")

				flagMapped := sam.Flag{}
				flagMapped.SetPaired()
				flagMapped.SetMateUnmapped()
				flagMapped.SetFirstInPair()

				isForward := index.SequenceIndexIsForward(uint32(matchesPerGenePosition[genePositions[0]][0].SequenceIndex))
				if !isForward {
					flagMapped.SetReverseStrand()
				}

				cigarString := strconv.Itoa(len(result.ReadR1.Sequence)) + "M"

				transcriptIds := make([]int, len(matchesPerGenePosition[genePositions[0]]))
				for i, match := range matchesPerGenePosition[genePositions[0]] {
					transcriptIds[i] = int(index.SequenceIndexToTranscriptIndex(uint32(match.SequenceIndex)))
				}

				mappedRecord := &sam.Record{
					Qname:         result.ReadR1.Header,
					Flag:          flagMapped,
					Rname:         index.Gene.Chromosome,
					Pos:           int(genePositions[0]),
					Mapq:          255, // TODO: figure this out
					Cigar:         cigarString,
					Rnext:         result.ReadR2.Header,
					Pnext:         0,
					Tlen:          0, // set to 0 because only one read was mapped
					Seq:           result.ReadR1.Sequence,
					Qual:          result.ReadR1.Quality,
					TranscriptIds: transcriptIds,
				}

				flagUnmapped := sam.Flag{}
				flagUnmapped.SetPaired()
				flagUnmapped.SetUnmapped()
				flagUnmapped.SetSecondInPair()

				unmappedRecord := &sam.Record{
					Qname:         result.ReadR2.Header,
					Flag:          flagUnmapped,
					Rname:         "*",
					Pos:           0,
					Mapq:          255,
					Cigar:         "*",
					Rnext:         result.ReadR1.Header,
					Pnext:         0,
					Tlen:          0,
					Seq:           result.ReadR2.Sequence,
					Qual:          result.ReadR2.Quality,
					TranscriptIds: nil,
				}

				recordPairs = append(recordPairs, &sam.RecordPair{
					First:  mappedRecord,
					Second: unmappedRecord,
				})
			} else {
				// only R2 was mapped

				logrus.Debug("only R2 was mapped")

				flagMapped := sam.Flag{}
				flagMapped.SetPaired()
				flagMapped.SetMateUnmapped()
				flagMapped.SetSecondInPair()

				isForward := index.SequenceIndexIsForward(uint32(matchesPerGenePosition[genePositions[0]][0].SequenceIndex))
				if !isForward {
					flagMapped.SetReverseStrand()
				}

				cigarString := strconv.Itoa(len(result.ReadR1.Sequence)) + "M"

				transcriptIds := make([]int, len(matchesPerGenePosition[genePositions[0]]))
				for i, match := range matchesPerGenePosition[genePositions[0]] {
					transcriptIds[i] = int(index.SequenceIndexToTranscriptIndex(uint32(match.SequenceIndex)))
				}

				mappedRecord := &sam.Record{
					Qname:         result.ReadR2.Header,
					Flag:          flagMapped,
					Rname:         index.Gene.Chromosome,
					Pos:           int(genePositions[0]),
					Mapq:          255,
					Cigar:         cigarString,
					Rnext:         result.ReadR1.Header,
					Pnext:         0,
					Tlen:          0, // set to 0 because only one read was mapped
					Seq:           result.ReadR2.Sequence,
					Qual:          result.ReadR2.Quality,
					TranscriptIds: transcriptIds,
				}

				flagUnmapped := sam.Flag{}
				flagUnmapped.SetPaired()
				flagUnmapped.SetUnmapped()
				flagUnmapped.SetFirstInPair()

				unmappedRecord := &sam.Record{
					Qname:         result.ReadR1.Header,
					Flag:          flagUnmapped,
					Rname:         "*",
					Pos:           0,
					Mapq:          255,
					Cigar:         "*",
					Rnext:         result.ReadR2.Header,
					Pnext:         0,
					Tlen:          0,
					Seq:           result.ReadR1.Sequence,
					Qual:          result.ReadR1.Quality,
					TranscriptIds: nil,
				}

				recordPairs = append(recordPairs, &sam.RecordPair{
					First:  unmappedRecord,
					Second: mappedRecord,
				})
			}

		} else {
			// the read was not paired so only R1 is mapped

			logrus.Debug("read was not paired so only R1 was mapped")

			flagMapped := sam.Flag{}
			isForward := index.SequenceIndexIsForward(uint32(matchesPerGenePosition[genePositions[0]][0].SequenceIndex))
			if !isForward {
				flagMapped.SetReverseStrand()
			}

			cigarString := strconv.Itoa(len(result.ReadR1.Sequence)) + "M"

			transcriptIds := make([]int, len(matchesPerGenePosition[genePositions[0]]))
			for i, match := range matchesPerGenePosition[genePositions[0]] {
				transcriptIds[i] = int(index.SequenceIndexToTranscriptIndex(uint32(match.SequenceIndex)))
			}

			mappedRecord := &sam.Record{
				Qname:         result.ReadR1.Header,
				Flag:          flagMapped,
				Rname:         index.Gene.Chromosome,
				Pos:           int(genePositions[0]),
				Mapq:          255, // TODO: figure this out
				Cigar:         cigarString,
				Rnext:         "*",
				Pnext:         0,
				Tlen:          0, // set to 0 because only one read was mapped
				Seq:           result.ReadR1.Sequence,
				Qual:          result.ReadR1.Quality,
				TranscriptIds: transcriptIds,
			}

			recordPairs = append(recordPairs, &sam.RecordPair{
				First:  mappedRecord,
				Second: nil,
			})
		}
	}

	return &recordPairs
}

// MapRead maps a single read to the reference using the index
// The mapping algorithm achieves a sublinear time complexity by discarding reads that have too many mismatches
// prior to matching every character of the read to the reference.
// The read is split into k-mers, of which at most x are allowed to no match exactly to allow a max of x mismatches.
// If the length of the read can not be divided by k, then there are some bases left over and excluded in the exact
// matching step. This is to ensure that each k-mer has the exact length k.
// These left-over bases have to be extended in the same way that k-mers are extended when they are combined.
// The exact matches are counted per transcript (sequence) and if at any point there is no sequence which has at least
// (number of current k-mer - x) exact matches, the read is discarded.
// Example: Say we have a read of length 150, and we want k-mers of length 8 and allow 5 mismatches:
// 150 / 8 = 18 k-mers and 6 bases left over
// At most 5 of these 18 k-mers are allowed to not match exactly to the reference sequence.
// This does only mean that the read is not discarded directly, but it still can if there are more than 5 mismatches
// within any of the k-mers that could not be matched exactly.
// When matching the 1. k-mer, there must be at least one sequence which has max(0, 1 - 5) exact matches.
// This means that when matching the 6. k-mer, there must be at least one sequence with an exact match.
// When matching the 15. k-mer, there must be at least max(0, 15 - 5) = 10 exact matches in any of the sequences.
// A k-mer can be matched to multiple sequences, for each of which the number of matches is increased but one k-mer
// can never be counted twice for the same sequence even when there are multiple matches.
func MapRead(read *fastq.Read, index *index.GtaIndex, t *timer.Timer) (map[int][]*core.InexactMatchResult, bool) {

	// The final result of this method. Contains positions of the read on the reference sequence within the
	// allowed range of mismatches.
	// The key is the sequence index (transcript) and the value is the result of the inexact matching.
	inexactMatchResults := make(map[int][]*core.InexactMatchResult)

	// the maximum number of mismatches allowed given the length of the sequence and the error rate
	maxMismatches := int(float64(len(read.Sequence)) * config.ErrorRate())

	logrus.WithFields(logrus.Fields{
		"maxMismatches": maxMismatches,
		"numKmers":      len(read.Sequence) / int(index.KeywordTree.KeywordLength),
	}).Debug("Set discard constraints")

	// the current position within the read
	position := 0
	// the state of the discard step, true if read is to be discarded
	failed := false

	// the matches per sequence index, sequence index is the key, the value is a list of matches for this sequence
	matchesPerSequenceMap := make(map[int]*core.SequenceMatches, index.NumSequences)

	for i := 0; i < index.NumSequences; i++ {
		matchesPerSequenceMap[i] = &core.SequenceMatches{
			SequenceIndex: i,
			NumMatches:    0,
			Matches:       make([]*core.SequenceMatch, 0),
		}
	}

	// the number of processed k-mers, used to determine the minimum number of matches a sequence must have
	countKmer := 0

	// generates every k-mer from the read and exact matches it to the reference
	// it is guaranteed that each k-mer has the length k, if the last k-mer is shorter than k, it is ignored for
	// the exact matching and must later be extended in the same way that k-mers are extended when they are combined
	for position < len(read.Sequence)-int(index.KeywordTree.KeywordLength) {

		countKmer += 1

		// the actual kmer sequence
		kmer := read.Sequence[position : position+int(index.KeywordTree.KeywordLength)]

		//logrus.WithFields(logrus.Fields{
		//	"position": position,
		//	"length":   currentLenKmer,
		//	"kmer":     kmer,
		//}).Info("create new kmer")

		// the result of the exact matching of this kmer using the suffix tree
		//mappingResult := index.SuffixTree.FindPatternExact(&kmer)
		timerStart := time.Now()
		mappingResult := index.KeywordTree.FindKeyword(&kmer)
		t.ExactMatch += time.Since(timerStart)

		// contains true for each sequence index that has a match
		// used to count this k-mer only once for each sequence
		sequencesMatched := make(map[int]bool)

		// add the matches to each sequence if there are any
		if mappingResult != nil {

			for _, match := range mappingResult.Matches {

				// the sequence was already discarded and must not be considered anymore
				if _, sequenceIsInMap := matchesPerSequenceMap[match.SequenceIndex]; !sequenceIsInMap {
					continue
				}

				// set the source positions for the match
				match.FromSource = position
				match.ToSource = position + int(index.KeywordTree.KeywordLength)

				// the candidate position of this read based on the matched position on the reference (FromTarget)
				// and the offset within the read (FromSource)
				candidatePosition := match.FromTarget - match.FromSource

				// if the candidate position is negative, discard the match
				if candidatePosition < 0 {
					continue
				}
				// if the candidate position plus the length of the read exceeds the length of the transcript, discard the match
				if candidatePosition+len(read.Sequence) > len(*index.GetSequenceByIndex(match.SequenceIndex)) {
					continue
				}

				// add the match to the sequence
				matchesPerSequenceMap[match.SequenceIndex].Matches = append(matchesPerSequenceMap[match.SequenceIndex].Matches, match)
				// mark the sequence as matched
				sequencesMatched[match.SequenceIndex] = true
			}

			// increase the number of distinct k-mer matches for each sequence that has a match
			for sequenceIndex, _ := range sequencesMatched {
				matchesPerSequenceMap[sequenceIndex].NumMatches += 1
			}
		}

		// the minimum number of matches a sequences has to have to survive this round
		minMatches := countKmer - maxMismatches

		// only discard sequences if there must be at least one match
		if minMatches > 0 {
			// check each sequence if it must be discarded if it has not gathered a match in this round
			for sequenceIndex := 0; sequenceIndex < index.NumSequences; sequenceIndex++ {
				// sequence was already discarded
				if _, sequenceIsInMap := matchesPerSequenceMap[sequenceIndex]; !sequenceIsInMap {
					continue
				}
				// the sequence has a match with the current k-mer and survives another round
				if _, sequenceHasMatch := sequencesMatched[sequenceIndex]; sequenceHasMatch {
					continue
				}
				// the sequence has not gained a match with the current k-mer and is therefore a potential discard candidate
				// discard the sequence if it has less than the required number of matches
				if matchesPerSequenceMap[sequenceIndex].NumMatches < minMatches {

					//logrus.WithFields(logrus.Fields{
					//	"sequenceIndex": sequenceIndex,
					//	"round":         countKmer,
					//	"minMatches":    minMatches,
					//	"matches":       matchesPerSequenceMap[sequenceIndex].NumMatches,
					//}).Debug("discard sequence index because too few matches")

					delete(matchesPerSequenceMap, sequenceIndex)
				}
			}
		}

		//for sequenceIndex, val := range matchesPerSequenceMap {
		//	logrus.WithFields(logrus.Fields{
		//		"sequenceIndex": sequenceIndex,
		//		"numMatches":    val.NumMatches,
		//	}).Debug("sequence info after k-mer matching")
		//}

		// stop if there are no more sequences left and the read can be discarded
		if len(matchesPerSequenceMap) == 0 {
			logrus.Debug("no more sequences left")

			failed = true
			break
		}

		position += int(index.KeywordTree.KeywordLength)
	}

	// discard the read because each sequence exceeded the allowed number of mismatches
	if failed {
		logrus.Debug("discard this read!")
		return inexactMatchResults, false
	} else {
		logrus.Debug("keep this read!")
	}

	// Idea:
	// Each exact k-mer hit provides a candidate position for the read on the reference sequence.
	// When the k-mer starts at position 0 in the read (FromSource), the candidate position of the read on
	// the reference sequence is the position of the exact match in the reference sequence (FromTarget).
	// When the k-mer does not start at position 0 in the read, the candidate position of the read on the reference
	// sequence is the position of the exact match in the reference sequence minus the position of the k-mer in the read.
	// Example: The k-mer starts at position 25 in the read and the exact match is at position 100 in the reference
	// sequence. The candidate position of the read on the reference sequence is 100 - 25 = 75.
	// Compute all unmatched positions for each candidate position and sequence index.
	// Combine consecutive matches to contexts on each sequence index.
	// Extend the matches to the left and right to find the best inexact-match.

	for sequenceIndex, matches := range matchesPerSequenceMap {

		//fmt.Println("sequenceIndex", sequenceIndex)
		//fmt.Println("numMatches", matches.NumMatches)

		candidatePositions := make(map[int][]*interval.Interval)
		candidatePositionsEcIds := make(map[int][]uint32)

		// collect the exact matched intervals per candidate position
		for _, match := range matches.Matches {

			// the candidate position of the read on the reference sequence (transcript)
			candidatePosition := match.FromTarget - match.FromSource

			// create the candidate position if it does not exist yet
			if _, ok := candidatePositions[candidatePosition]; !ok {
				candidatePositions[candidatePosition] = make([]*interval.Interval, 0)
				candidatePositionsEcIds[candidatePosition] = make([]uint32, 0)
			}
			// add the interval for the current candidate position
			candidatePositions[candidatePosition] = append(candidatePositions[candidatePosition], &interval.Interval{
				Start: match.FromTarget,
				End:   match.ToTarget,
			})
			candidatePositionsEcIds[candidatePosition] = append(candidatePositionsEcIds[candidatePosition], match.EquivalenceClassIds...)
		}

		// sort and merge all overlapping intervals for each candidate position
		for candidatePosition, intervals := range candidatePositions {

			// contains consecutive intervals of exactly matched positions
			matchedRegions := interval.MergeIntervals(intervals)

			logrus.Debug("starting inexact matching in unmatched regions")

			// Idea:
			// Determine all unmatched positions for each candidate position by iterating the sorted list of intervals
			// and match them to the reference sequence.
			// There are three possible locations for unmatched regions as there must always be at least one exact matched
			// interval for each candidate position:
			// 1. Between the candidate position and the start of the first exact matched interval
			// 2. Between two exact matched intervals
			// 3. Between the end of the last exact matched interval and the position (candidate position + read length)
			// Each position in an unmatched region is compared to the reference sequence and the number of mismatches
			// is counted. If the number of mismatches exceeds the allowed number of mismatches, the candidate position
			// is discarded.

			// the mismatched source positions (read) for this candidate position
			mismatches := make([]int, 0)
			// true if the maximum number of mismatches is exceeded
			failed := false

			// location 1: between the candidate position and the start of the first exact matched interval
			if matchedRegions[0].Start > candidatePosition {

				logrus.WithFields(logrus.Fields{
					"sequenceIndex":     sequenceIndex,
					"candidatePosition": candidatePosition,
					"start":             matchedRegions[0].Start,
				}).Debug("unmatched region before first match")

				// the target position is the position within the transcript sequence
				for targetPosition := candidatePosition; targetPosition < matchedRegions[0].Start; targetPosition++ {

					// the position within the read sequence
					readOffset := targetPosition - candidatePosition

					if read.Sequence[readOffset] != (*index.GetSequenceByIndex(sequenceIndex))[targetPosition] {

						logrus.WithFields(logrus.Fields{
							"targetPosition": targetPosition,
							"readOffset":     readOffset,
							"readBase":       string(read.Sequence[readOffset]),
							"refBase":        string((*index.GetSequenceByIndex(sequenceIndex))[targetPosition]),
						}).Debug("mismatch")

						mismatches = append(mismatches, readOffset)
					}

					if len(mismatches) > maxMismatches {
						failed = true
						break
					}
				}
			}

			if failed {
				logrus.WithFields(logrus.Fields{
					"sequenceIndex":     sequenceIndex,
					"candidatePosition": candidatePosition,
					"numMismatches":     len(mismatches),
					"mismatches":        mismatches,
				}).Debug("discard this candidate position")
				continue
			}

			// location 2: between two exact matched intervals
			// the intervals are already sorted and merged (non-overlapping) which means that all regions between
			// two consecutive intervals are unmatched regions.
			for i := 0; i < len(matchedRegions)-1; i++ {

				logrus.WithFields(logrus.Fields{
					"sequenceIndex":     sequenceIndex,
					"candidatePosition": candidatePosition,
					"start":             matchedRegions[i].End,
					"end":               matchedRegions[i+1].Start,
				}).Debug("unmatched region between matched regions")

				// the target position is the position within the transcript sequence
				for targetPosition := matchedRegions[i].End; targetPosition < matchedRegions[i+1].Start; targetPosition++ {

					// the position within the read sequence
					readOffset := targetPosition - candidatePosition

					if read.Sequence[readOffset] != (*index.GetSequenceByIndex(sequenceIndex))[targetPosition] {

						logrus.WithFields(logrus.Fields{
							"targetPosition": targetPosition,
							"readOffset":     readOffset,
							"readBase":       string(read.Sequence[readOffset]),
							"refBase":        string((*index.GetSequenceByIndex(sequenceIndex))[targetPosition]),
						}).Debug("mismatch")

						mismatches = append(mismatches, readOffset)
					}

					if len(mismatches) > maxMismatches {
						failed = true
						break
					}
				}

				if failed {
					break
				}
			}

			if failed {
				logrus.WithFields(logrus.Fields{
					"sequenceIndex":     sequenceIndex,
					"candidatePosition": candidatePosition,
					"numMismatches":     len(mismatches),
					"mismatches":        mismatches,
				}).Debug("discard this candidate position")
				continue
			}

			// location 3: between the end of the last exact matched interval and the position (candidate position + read length)
			if matchedRegions[len(matchedRegions)-1].End < candidatePosition+len(read.Sequence) {

				// the target position is the position within the transcript sequence
				for targetPosition := matchedRegions[len(matchedRegions)-1].End; targetPosition < candidatePosition+len(read.Sequence); targetPosition++ {

					readOffset := targetPosition - candidatePosition

					if read.Sequence[readOffset] != (*index.GetSequenceByIndex(sequenceIndex))[targetPosition] {

						logrus.WithFields(logrus.Fields{
							"targetPosition": targetPosition,
							"readOffset":     readOffset,
							"readBase":       string(read.Sequence[readOffset]),
							"refBase":        string((*index.GetSequenceByIndex(sequenceIndex))[targetPosition]),
						}).Debug("mismatch")

						mismatches = append(mismatches, readOffset)
					}

					if len(mismatches) > maxMismatches {
						failed = true
						break
					}
				}
			}

			if failed {
				logrus.WithFields(logrus.Fields{
					"sequenceIndex":     sequenceIndex,
					"candidatePosition": candidatePosition,
					"numMismatches":     len(mismatches),
					"mismatches":        mismatches,
				}).Debug("discard this candidate position")
				continue
			}

			logrus.WithFields(logrus.Fields{
				"sequenceIndex":     sequenceIndex,
				"candidatePosition": candidatePosition,
				"numMismatches":     len(mismatches),
				"mismatches":        mismatches,
			}).Debug("candidate position passed all tests")

			// the candidate position passed all tests and is added to the final result

			// make the equivalence class id list unique
			uniqueEcIds := make(map[uint32]bool)
			for _, ecId := range candidatePositionsEcIds[candidatePosition] {
				uniqueEcIds[ecId] = true
			}
			ecIds := make([]uint32, 0)
			for ecId, _ := range uniqueEcIds {
				ecIds = append(ecIds, ecId)
			}

			// determine the start position relative to the gene
			transcriptIndex := index.SequenceIndexToTranscriptIndex(uint32(sequenceIndex))
			startPosition := index.TranslateRelativeTranscriptPositionToRelativeGenePosition(transcriptIndex, uint32(candidatePosition))

			// create a new list of inexact match results if it does not exist yet for this sequence index
			if _, ok := inexactMatchResults[sequenceIndex]; !ok {
				inexactMatchResults[sequenceIndex] = make([]*core.InexactMatchResult, 0)
			}
			// add the current candidate position to the list of inexactly matched positions for this sequence index
			inexactMatchResults[sequenceIndex] = append(inexactMatchResults[sequenceIndex],
				&core.InexactMatchResult{
					SequenceIndex:       sequenceIndex,
					FromTarget:          candidatePosition,
					ToTarget:            candidatePosition + len(read.Sequence),
					Mismatches:          mismatches,
					EquivalenceClassIds: ecIds,
					FromGene:            startPosition,
				})
		}
	}

	// return the map of inexact matching results and true if there are any results, else false
	return inexactMatchResults, len(inexactMatchResults) > 0
}
