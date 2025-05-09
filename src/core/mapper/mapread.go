package mapper

import (
	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/sirupsen/logrus"
)

func MapRead(read *fastq.Read, genomeIndex *index.GenomeIndex) ([]mapperutils.ReadMatchResult, bool) {

	logrus.WithFields(logrus.Fields{
		"read":   read.Header,
		"length": len(*read.Sequence),
	}).Debug("Mapping read")

	globalMatches := mapperutils.GlobalMatchResult{
		MatchesPerSequence: make([]*mapperutils.SequenceMatchResult, genomeIndex.KeywordTree.NumSequences),
	}

	// generate non-overlapping k-mers for the read pair
	// skip last kmer if length is not divisible by kmer length
	for i := 0; i <= len(*read.Sequence)-(int(config.KmerLength())); i += int(config.KmerLength()) {

		kmer := (*read.Sequence)[i : i+int(config.KmerLength())]

		// TODO: replace by hashmap
		matches := genomeIndex.KeywordTree.FindKeyword(&kmer, i)

		if matches == nil {
			continue
		}

		for _, match := range matches {

			// add new sequence to global matches
			if globalMatches.MatchesPerSequence[match.SequenceIndex] == nil {
				globalMatches.MatchesPerSequence[match.SequenceIndex] = &mapperutils.SequenceMatchResult{
					MatchesPerDiagonal: make(map[int][]*mapperutils.Match),
				}
			}

			// add new diagonal to sequence match
			if globalMatches.MatchesPerSequence[match.SequenceIndex].MatchesPerDiagonal[match.StartGenome] == nil {
				globalMatches.MatchesPerSequence[match.SequenceIndex].MatchesPerDiagonal[match.StartGenome] =
					make([]*mapperutils.Match, 0)
			}

			// add match to diagonal
			globalMatches.MatchesPerSequence[match.SequenceIndex].MatchesPerDiagonal[match.StartGenome] =
				append(globalMatches.MatchesPerSequence[match.SequenceIndex].MatchesPerDiagonal[match.StartGenome], match)
		}
	}

	// list has size of number of sequences in index
	// each element represents the maximum number of kmers that matched exactly on the same diagonal
	// TODO: change to a hashmap or use info from index
	maxDiagonalHitsPerSequence := make([]int, genomeIndex.KeywordTree.NumSequences)

	for seqIndex, sequenceMatches := range globalMatches.MatchesPerSequence {

		maxHits := 0

		if sequenceMatches != nil {
			for _, matches := range sequenceMatches.MatchesPerDiagonal {
				if len(matches) > maxHits {
					maxHits = len(matches)
				}
			}
		}

		maxDiagonalHitsPerSequence[seqIndex] = maxHits
	}

	seqIndexSorted := sortedIndicesDesc(maxDiagonalHitsPerSequence)

	// minimum number of hits any diagonal must have to be considered
	minHits := 5

	for i, seqIndex := range seqIndexSorted {
		if maxDiagonalHitsPerSequence[seqIndex] < minHits {
			seqIndexSorted[i] = -1
		}
	}

	logrus.WithFields(logrus.Fields{
		"maxDiagonalHitsPerSequence": maxDiagonalHitsPerSequence,
		"sortedByIndex":              seqIndexSorted,
	}).Debug("Potential sequences (-1 index not used)")

	results := make([]mapperutils.ReadMatchResult, 0)

	//fmt.Println(read.Header)

sequenceLoop:
	for _, seqIndex := range seqIndexSorted {

		if seqIndex == -1 {
			// -1 means this sequence was disregarded because of some constraints
			continue
		}

		sequenceMatches := globalMatches.MatchesPerSequence[seqIndex]

		if sequenceMatches == nil {
			continue
		}

		logrus.WithFields(logrus.Fields{
			"seqIndex": seqIndex,
		}).Debug("using sequence to compute best match")

		// 1. find the best diagonal (most matches)
		// 2. find gaps in diagonal (unmatched regions via regionvector)
		// 3. fill gaps in digonal (each position for a match or mismatch)
		// 4. repeat 1-4 until there are no fitting diagonals left
		// now there are multiple options left:
		// - the read is fully matched already -> done
		// - there are multiple diagonals with unmatched positions in between
		//   (probably due to a junction -> extend both diagonals as far as possible)
		// - there are unmatched positions in the front or back of the read
		//   (extend but if not possible then search for another diagonal with shorter kmer)

		result := mapperutils.ReadMatchResult{
			SequenceIndex:  seqIndex,
			MatchedRead:    regionvector.NewRegionVector(),
			MatchedGenome:  regionvector.NewRegionVector(),
			MismatchesRead: make([]int, 0),
			FourthPass:     false,
		}

		diagonalHandler := mapperutils.NewDiagonalHandlerWithData(sequenceMatches.MatchesPerDiagonal)

		logrus.Debug("Handler initialized with the following diagonals")
		for k, v := range diagonalHandler.Diagonals {
			logrus.WithFields(logrus.Fields{
				"posGenome": k,
				"matches":   v,
				"length":    len(v),
			}).Debug("diagonal")
		}

		for result.MatchedRead.Length() < len(*read.Sequence) {

			logrus.Debug("searching for next best diagonal")

			bestDiagonal, bestDiagonalLength, foundDiagonal := diagonalHandler.GetBestDiagonal()

			// exit the loop if no diagonal is found to handle the unmatched regions and fill the match
			if !foundDiagonal {
				logrus.Debug("no suitable diagonal found")
				break
			}

			// if the diagonal is negative or larger than the genome length, this could be a sign of clipping
			// because the match exceeds the target sequence as it is in the index
			if bestDiagonal < 0 {
				logrus.Debug("genome index out of bounds left side (clipping?)")
				break
			}
			if bestDiagonal+len(*read.Sequence) > len(*genomeIndex.Sequences[seqIndex]) {
				logrus.Debug("genome index out of bounds right side (clipping?)")
				break
			}

			logrus.WithFields(logrus.Fields{
				"posGenome": bestDiagonal,
				"length":    bestDiagonalLength,
				"matches":   sequenceMatches.MatchesPerDiagonal[bestDiagonal],
			}).Debug("considering diagonal")

			isDiagonalValid := diagonalHandler.IsValidExtension(sequenceMatches.MatchesPerDiagonal[bestDiagonal], result, read)

			if !isDiagonalValid {

				logrus.Debug("diagonal is not valid")

				// if the extension is not valid, remove from diags
				// but do not consume the kmers, since they could be placed at an other spot maybe
				delete(diagonalHandler.Diagonals, bestDiagonal)
				continue
			}

			logrus.Debug("diagonal is valid")

			diagonalRead := regionvector.NewRegionVector()
			diagonalGenome := regionvector.NewRegionVector()

			for _, match := range sequenceMatches.MatchesPerDiagonal[bestDiagonal] {
				// the region can not be part of another diagonal that is already used
				if match.Used {
					continue
				}
				diagonalRead.AddRegionNonOverlappingPanic(match.FromRead, match.ToRead)
				diagonalGenome.AddRegionNonOverlappingPanic(match.FromGenome, match.ToGenome)
			}

			logrus.WithFields(logrus.Fields{
				"posGenome":  bestDiagonal,
				"numMatches": bestDiagonalLength,
				"read":       diagonalRead,
				"genome":     diagonalGenome,
				"gaps":       len(diagonalRead.Regions) - 1,
			}).Debug("found best diagonal")

			// find gaps in diagonal and fill them (matches and mismatches)
			mismatches := make([]int, 0)

			foundGap := diagonalRead.HasGaps()

			// resolve gaps on the same diagonal
			// gaps can occur because of mismatches in the read and genome which prevent exact kmer matching
			// this is done by filling the gaps in the read and genome and counting the mismatches
			// the regionvectors of read and genome should have the same length as they are coupled
			// because they are part of the same diagonal (no indels, otherwise not on the same diagonal)
			for diagonalRead.HasGaps() {

				gapRead := diagonalRead.GetFirstGap()
				gapGenome := diagonalGenome.GetFirstGap()

				logrus.WithFields(logrus.Fields{
					"read":   gapRead,
					"genome": gapGenome,
				}).Debug("found gap")

				// fill the gap in the read by adding the gap as region
				diagonalRead.AddRegionNonOverlappingPanic(gapRead.Start, gapRead.End)

				geneSeqPos := 0
				for i := gapRead.Start - 1; i < gapRead.End-1; i++ {

					readByte := (*read.Sequence)[i]
					gIndex := gapGenome.Start + geneSeqPos - 1
					genomeByte := (*genomeIndex.Sequences[seqIndex])[gIndex]
					geneSeqPos++

					if readByte != genomeByte {
						mismatches = append(mismatches, i)

						logrus.Debug("mismatch at read position: ", i)

						// TODO: make max allowed mismatches configurable
						if len(mismatches) > 5 {
							logrus.Debug("too many mismatches")
							continue sequenceLoop
						}
					}
				}

				diagonalGenome.AddRegionNonOverlappingPanic(gapGenome.Start, gapGenome.End)

				// also update used status in kmers that were not part of the best diag (but existed as geps inside the best diag)
				// this way, these kmers cant be used in other diagonals
				diagonalHandler.ConsumeKmer(gapRead.Start, gapRead.End, gapGenome.Start, gapGenome.End)
			}

			if foundGap {
				logrus.WithFields(logrus.Fields{
					"read":       diagonalRead,
					"genome":     diagonalGenome,
					"mismatches": mismatches,
				}).Debug("filled gap")
			}

			// add digonal to result (diagonal rv should have length 1 after gap filling)
			result.MatchedRead.AddRegionNonOverlappingPanic(diagonalRead.GetFirstRegion().Start,
				diagonalRead.GetLastRegion().End)

			result.MatchedGenome.AddRegionNonOverlappingPanic(diagonalGenome.GetFirstRegion().Start,
				diagonalGenome.GetLastRegion().End)

			result.MismatchesRead = append(result.MismatchesRead, mismatches...)

			// remove diagonal from handler
			diagonalHandler.ConsumeDiagonal(bestDiagonal)
		}

		logrus.WithFields(logrus.Fields{
			"read":   result.MatchedRead,
			"genome": result.MatchedGenome,
		}).Debug("done processing diagonals")

		// if no diagonal was valid, skip this read
		if result.MatchedRead.Length() == 0 {
			logrus.Debug("no diagonal was valid")
			continue sequenceLoop
		}

		if result.MatchedRead.Length() == len(*read.Sequence) {
			logrus.Debug("read fully matched")
		} else {

			if result.MatchedRead.HasGaps() {

				logrus.Debug("unmatched regions within read")

				// used to keep track of the read position for the next gap
				readGapPos := 0
				indexRegionBeforeGap := result.MatchedRead.GetGapIndexAfterPos(readGapPos)

				for indexRegionBeforeGap > -1 {

					gapRead := result.MatchedRead.GetGapAfterRegionIndex(indexRegionBeforeGap)
					gapGenome := result.MatchedGenome.GetGapAfterRegionIndex(indexRegionBeforeGap)

					if gapGenome == nil {
						logrus.WithFields(logrus.Fields{
							"read":    result.MatchedRead,
							"genome":  result.MatchedGenome,
							"gapRead": gapRead,
						}).Debug("insertion found")
					} else {

						logrus.WithFields(logrus.Fields{
							"gapRead":   gapRead,
							"gapGenome": gapGenome,
						}).Debug("found gap to be handled")

						// TODO: handle insertions
						// when the gap in the read is larger than the gap in the genome
						// there is maybe an insertion in the read
						// currently the assumption is that this is a non-mapping read
						if gapRead.Length() > gapGenome.Length() {
							logrus.WithFields(logrus.Fields{
								"gapRead":   gapRead,
								"gapGenome": gapGenome,
							}).Debug("gap read is larger than gap genome")

							continue sequenceLoop
						}

						bestSplit := determineBestSplit(genomeIndex, read, seqIndex, gapRead, gapGenome)

						if bestSplit == -1 {
							logrus.WithFields(logrus.Fields{
								"qname": read.Header,
							}).Fatal("no best split found")
						}

						logrus.WithFields(logrus.Fields{
							"split": bestSplit,
						}).Debug("best split found")

						logrus.WithFields(logrus.Fields{
							"read":   result.MatchedRead,
							"genome": result.MatchedGenome,
						}).Debug("regions before")

						// when bestSplit is 0 then there is nothing to be added to the left side of the gap
						if bestSplit > 0 {
							result.MatchedRead.AddRegionNonOverlappingPanic(gapRead.Start, gapRead.Start+bestSplit)
							result.MatchedGenome.AddRegionNonOverlappingPanic(gapGenome.Start, gapGenome.Start+bestSplit)
						}

						logrus.WithFields(logrus.Fields{
							"read":   result.MatchedRead,
							"genome": result.MatchedGenome,
						}).Debug("regions after left")

						// when bestSplit is equal to the length of the gap then there is nothing
						// to be added to the right side of the gap
						if bestSplit < gapRead.Length() {
							result.MatchedRead.AddRegionNonOverlappingPanic(gapRead.End-(gapRead.Length()-bestSplit), gapRead.End)
							result.MatchedGenome.AddRegionNonOverlappingPanic(gapGenome.End-(gapRead.Length()-bestSplit), gapGenome.End)
						}

						logrus.WithFields(logrus.Fields{
							"read":   result.MatchedRead,
							"genome": result.MatchedGenome,
						}).Debug("regions after right")
					}

					readGapPos = gapRead.End + 1
					indexRegionBeforeGap = result.MatchedRead.GetGapIndexAfterPos(readGapPos)
				}
			}

			if result.MatchedRead.GetFirstRegion().Start > 0 {

				startRead := 0
				endRead := result.MatchedRead.GetFirstRegion().Start
				extensionLength := endRead - startRead

				startGenome := result.MatchedGenome.GetFirstRegion().Start - extensionLength

				readSequence := (*read.Sequence)[startRead:endRead]

				logrus.WithFields(logrus.Fields{
					"startRead":   startRead,
					"endRead":     endRead,
					"startGenome": startGenome,
					"endGenome":   startGenome + len(readSequence),
				}).Debug("unmatched positions in front of read")

				// TODO: this could be a use case for clipping
				// when a read maps to the target sequence but would go out of bounds
				// it could still be a valid mapping but it needs to be clipped
				if startGenome < 0 {
					logrus.Debug("genome index out of bounds")
					continue sequenceLoop
				}

				genomeSequence := (*genomeIndex.Sequences[seqIndex])[startGenome : startGenome+len(readSequence)]

				numMismatches := 0

				// trying to extend the read to the left
				for i := 0; i < extensionLength; i++ {
					if readSequence[i] != genomeSequence[i] {
						numMismatches++
					}
					if numMismatches > 5 {
						break
					}
				}

				// add to fourth pass if there are too many mismatches
				if numMismatches > 5 {
					result.FourthPass = true
					results = append(results, result)
					continue sequenceLoop
				}

				result.MatchedRead.AddRegionNonOverlappingPanic(startRead, endRead)
				result.MatchedGenome.AddRegionNonOverlappingPanic(startGenome, startGenome+extensionLength)
			}

			if result.MatchedRead.GetLastRegion().End < len(*read.Sequence) {

				startRead := result.MatchedRead.GetLastRegion().End
				endRead := len(*read.Sequence)
				extensionLength := endRead - startRead

				startGenome := result.MatchedGenome.GetLastRegion().End

				readSequence := (*read.Sequence)[startRead:endRead]

				logrus.WithFields(logrus.Fields{
					"startRead":   startRead,
					"endRead":     endRead,
					"startGenome": startGenome,
					"endGenome":   startGenome + len(readSequence),
				}).Debug("unmatched positions in back of read")

				// TODO: this could be a use case for clipping
				// when a read maps to the target sequence but would go out of bounds
				// it could still be a valid mapping but it needs to be clipped
				if startGenome+len(readSequence) > len(*genomeIndex.Sequences[seqIndex]) {
					logrus.Debug("genome index out of bounds")
					continue sequenceLoop
				}

				genomeSequence := (*genomeIndex.Sequences[seqIndex])[startGenome : startGenome+len(readSequence)]

				//fmt.Println("seqIndex: ", seqIndex)
				//fmt.Println("startRead: ", startRead)
				//fmt.Println("endRead: ", endRead)
				//fmt.Println("extensionLength: ", extensionLength)
				//fmt.Println("startGenome: ", startGenome)
				//fmt.Println("readSequence:\t", string(readSequence))
				//fmt.Println("genomeSequence:\t", string(genomeSequence))

				numMismatches := 0

				// trying to extend the read to the right
				for i := 0; i < extensionLength; i++ {
					if readSequence[i] != genomeSequence[i] {
						numMismatches++
					}
					if numMismatches > 5 {
						break
					}
				}

				// add to fourth pass if there are too many mismatches
				if numMismatches > 5 {
					result.FourthPass = true
					results = append(results, result)
					continue sequenceLoop
				}

				result.MatchedRead.AddRegionNonOverlappingPanic(startRead, endRead)
				result.MatchedGenome.AddRegionNonOverlappingPanic(startGenome, startGenome+extensionLength)
			}
		}

		results = append(results, result)
	}

	if len(results) == 0 {
		logrus.Debug("no results found")
		return results, false
	}
	return results, true
}

func determineBestSplit(
	genomeIndex *index.GenomeIndex,
	read *fastq.Read,
	seqIndex int,
	gapRead *regionvector.Region,
	gapGenome *regionvector.Region) int {

	logrus.WithFields(logrus.Fields{
		"gapRead":   gapRead,
		"gapGenome": gapGenome,
	}).Debug("determining best split")

	// cululative mismatch count for the left and right extensions
	// for lErrors the index i represents the number of mismatches for the first i positions of the extension
	// for rErrors the index i represents the number of mismatches for the last i positions of the extension
	lErrors := make([]int, gapRead.Length()+1)
	rErrors := make([]int, gapRead.Length()+1)

	lErrors[0] = 0
	rErrors[0] = 0

	for i := 1; i <= gapRead.Length(); i++ {
		lErrors[i] = lErrors[i-1]
		if (*read.Sequence)[gapRead.Start+i-1] != (*genomeIndex.Sequences[seqIndex])[gapGenome.Start+i-1] {
			lErrors[i]++
		}

		rErrors[i] = rErrors[i-1]
		if (*read.Sequence)[gapRead.End-i] != (*genomeIndex.Sequences[seqIndex])[gapGenome.End-i] {
			rErrors[i]++
		}
	}

	logrus.WithFields(logrus.Fields{
		"lErrors": lErrors,
		"rErrors": rErrors,
	}).Debug("determined mismatches")

	// the minimum number of mismatches
	// the +2 is based on the maximum penalty returned by scoreSpliceSites()
	minErrors := lErrors[gapRead.Length()] + rErrors[gapRead.Length()] + 2
	// the position of the split with the minimum number of mismatches
	minSplit := -1

	// TODO: keep track of the actual mismatch positions
	// TODO: if no suitable split is found then:
	// - maybe there is another exon in between if enough bases missing from read
	// - maybe keep the readpair for fourth pass
	for i := 0; i <= gapRead.Length(); i++ {

		lPos := i
		rPos := gapRead.Length() - i

		numMismatches := lErrors[lPos] + rErrors[rPos]

		donorSiteStart := gapGenome.Start + i
		donorSiteSeq := (*genomeIndex.Sequences[seqIndex])[donorSiteStart : donorSiteStart+2]

		splitRev := gapRead.Length() - i
		acceptorSiteStart := gapGenome.End - splitRev
		acceptorSiteSeq := (*genomeIndex.Sequences[seqIndex])[acceptorSiteStart-2 : acceptorSiteStart]

		isForwardStrand := genomeIndex.IsSequenceForward(seqIndex)

		// add a penalty if the splice site is not canonical
		// 2 means that there is no known splice site
		spliceSitePenalty := scoreSpliceSites(donorSiteSeq[0], donorSiteSeq[1],
			acceptorSiteSeq[0], acceptorSiteSeq[1], isForwardStrand)

		numMismatches += spliceSitePenalty

		logrus.WithFields(logrus.Fields{
			"split":               i,
			"splice site penalty": spliceSitePenalty,
			"numMismatches":       numMismatches,
		}).Debug("possible split")

		if numMismatches < minErrors {
			minErrors = numMismatches
			minSplit = i
		}
	}

	return minSplit
}
