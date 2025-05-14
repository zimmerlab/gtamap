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

		// divide each kmer match into sequence and diagonal per sequence
		// a sequence is the forward (as given) and rev-comp sequence of each given input sequence
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

	logrus.WithFields(logrus.Fields{
		"maxDiagonalHitsPerSequence": maxDiagonalHitsPerSequence,
		"sortedByIndex":              seqIndexSorted,
	}).Debug("Potential sequences")

	results := make([]mapperutils.ReadMatchResult, 0)

sequenceLoop:
	for _, seqIndex := range seqIndexSorted {

		if seqIndex == -1 {
			// -1 means this sequence was disregarded because of some constraints
			// that can be defined previously to skip certain sequences
			continue
		}

		sequenceMatches := globalMatches.MatchesPerSequence[seqIndex]

		if sequenceMatches == nil {
			// skip this sequence if there are no matches found at all
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
			SecondPass:     false,
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
				// add the match to the diagonal
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

			gapsRead, gapsGenome := mapperutils.ComputeGapsInDiagonal(diagonalRead, diagonalGenome, &result)

			// resolve gaps on the same diagonal
			// gaps can occur because of mismatches in the read and genome which prevent exact kmer matching
			// this is done by filling the gaps in the read and genome and counting the mismatches
			// the regionvectors of read and genome should have the same length as they are coupled
			// because they are part of the same diagonal (no indels, otherwise not on the same diagonal)
			for i := 0; i < len(gapsRead.Regions); i++ {

				gapRead := gapsRead.Regions[i]
				gapGenome := gapsGenome.Regions[i]

				logrus.WithFields(logrus.Fields{
					"read":   gapRead,
					"genome": gapGenome,
				}).Debug("found gap")

				// fill the gap in the read by adding the gap as region
				diagonalRead.AddRegionNonOverlappingPanic(gapRead.Start, gapRead.End)

				// count mismatches when filling the gap and skip the match if there are too many
				geneSeqPos := 0
				for i := gapRead.Start - 1; i < gapRead.End-1; i++ {

					readByte := (*read.Sequence)[i]
					gIndex := gapGenome.Start + geneSeqPos - 1
					genomeByte := (*genomeIndex.Sequences[seqIndex])[gIndex]
					geneSeqPos++

					// skip matches
					if readByte == genomeByte {
						continue
					}

					// add the mismatche to the result
					result.MismatchesRead = append(result.MismatchesRead, i)

					// skip this match result if there are too many mismatches
					if exceedsMismatchConstraint(read, result) {

						logrus.WithFields(logrus.Fields{
							"mismatchPercentage":    float64(len(result.MismatchesRead)) * 100 / float64(len(*read.Sequence)),
							"maxMismatchPercentage": config.MaxMismatchPercentage(),
							"mismatches":            result.MismatchesRead,
							"numMismatches":         len(result.MismatchesRead),
						}).Debug("too many mismatches in diagonal filling -> skip sequence")

						continue sequenceLoop
					}
				}

				diagonalGenome.AddRegionNonOverlappingPanic(gapGenome.Start, gapGenome.End)

				// also update used status in kmers that were not part of the best diag (but existed as geps inside the best diag)
				// this way, these kmers cant be used in other diagonals
				diagonalHandler.ConsumeKmer(gapRead.Start, gapRead.End, gapGenome.Start, gapGenome.End)
			}

			if len(gapsRead.Regions) > 0 {
				logrus.WithFields(logrus.Fields{
					"read":       diagonalRead,
					"genome":     diagonalGenome,
					"mismatches": result.MismatchesRead,
				}).Debug("filled gap")
			}

			// add all matches in diagonal to the result
			for i := 0; i < len(diagonalRead.Regions); i++ {
				regionRead := diagonalRead.Regions[i]
				regionGenome := diagonalGenome.Regions[i]

				logrus.WithFields(logrus.Fields{
					"read":   regionRead,
					"genome": regionGenome,
				}).Debug("adding diagonal match to result")

				result.MatchedRead.AddRegionNonOverlappingPanic(regionRead.Start, regionRead.End)
				result.MatchedGenome.AddRegionNonOverlappingPanic(regionGenome.Start, regionGenome.End)
			}

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

			// there are gaps (unmatched regions) in the read (between matched regions / diagonals)
			if result.MatchedRead.HasGaps() {

				logrus.Debug("unmatched regions within read")

				// used to keep track of the read position for the next gap
				readGapPos := 0
				// returns the index of the first region after which a gap occurs (-1 if no gap)
				indexRegionBeforeGap := result.MatchedRead.GetGapIndexAfterPos(readGapPos)

				// loop through all gaps in the read (-1 means there is no more gap)
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
							// this should not happen because a split should be found every time
							// even if the split has a bad score (many mismatches, no splice sites, etc)
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

							// add the split to the result
							result.MatchedRead.AddRegionNonOverlappingPanic(gapRead.Start, gapRead.Start+bestSplit)
							result.MatchedGenome.AddRegionNonOverlappingPanic(gapGenome.Start, gapGenome.Start+bestSplit)

							// the read and genome sequences from the start of the gap to the best split (left)
							readByte := (*read.Sequence)[gapRead.Start : gapRead.Start+bestSplit]
							genomeByte := (*genomeIndex.Sequences[seqIndex])[gapGenome.Start : gapGenome.Start+bestSplit]

							// add the mismatches to the result
							for i := 0; i < bestSplit; i++ {
								// add the mismatche to the result
								if readByte[i] != genomeByte[i] {
									result.MismatchesRead = append(result.MismatchesRead, gapRead.Start+i)

									// skip this match result if there are too many mismatches
									if exceedsMismatchConstraint(read, result) {
										logrus.WithFields(logrus.Fields{
											"mismatchPercentage":    float64(len(result.MismatchesRead)) * 100 / float64(len(*read.Sequence)),
											"maxMismatchPercentage": config.MaxMismatchPercentage(),
											"mismatches":            result.MismatchesRead,
											"numMismatches":         len(result.MismatchesRead),
										}).Debug("too many mismatches in middle extension (left) -> skip sequence")
										continue sequenceLoop
									}
								}
							}
						}

						logrus.WithFields(logrus.Fields{
							"read":   result.MatchedRead,
							"genome": result.MatchedGenome,
						}).Debug("regions after left")

						// when bestSplit is equal to the length of the gap then there is nothing
						// to be added to the right side of the gap
						if bestSplit < gapRead.Length() {

							logrus.WithFields(logrus.Fields{
								"gapReadEnd":    gapRead.End,
								"bestSplit":     bestSplit,
								"gapReadLength": gapRead.Length(),
								"left":          gapRead.End - (gapRead.Length() - bestSplit),
								"right":         gapRead.End,
							}).Debug("debug split right")

							// add the split to the result
							result.MatchedRead.AddRegionNonOverlappingPanic(gapRead.End-(gapRead.Length()-bestSplit), gapRead.End)
							result.MatchedGenome.AddRegionNonOverlappingPanic(gapGenome.End-(gapRead.Length()-bestSplit), gapGenome.End)

							// the read and genome sequences from the best split to the end of the gap (right)
							readByte := (*read.Sequence)[gapRead.End-(gapRead.Length()-bestSplit) : gapRead.End]
							genomeByte := (*genomeIndex.Sequences[seqIndex])[gapGenome.End-(gapRead.Length()-bestSplit) : gapGenome.End]

							// add the mismatches to the result
							for i := 0; i < gapRead.Length()-bestSplit; i++ {
								// add the mismatche to the result
								if readByte[i] != genomeByte[i] {
									result.MismatchesRead = append(result.MismatchesRead, gapRead.End-(bestSplit-i))

									// skip this match result if there are too many mismatches
									if exceedsMismatchConstraint(read, result) {
										logrus.WithFields(logrus.Fields{
											"mismatchPercentage":    float64(len(result.MismatchesRead)) * 100 / float64(len(*read.Sequence)),
											"maxMismatchPercentage": config.MaxMismatchPercentage(),
											"mismatches":            result.MismatchesRead,
											"numMismatches":         len(result.MismatchesRead),
										}).Debug("too many mismatches in middle extension (right) -> skip sequence")
										continue sequenceLoop
									}
								}
							}
						}

						logrus.WithFields(logrus.Fields{
							"read":   result.MatchedRead,
							"genome": result.MatchedGenome,
						}).Debug("regions after right")
					}

					// determine the next gap (-1 if there is none)
					readGapPos = gapRead.End + 1
					indexRegionBeforeGap = result.MatchedRead.GetGapIndexAfterPos(readGapPos)
				}
			}

			// there are unmatched positions in front of the read
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

					// skip matches
					if readSequence[i] == genomeSequence[i] {
						continue
					}

					// add the mismatches to the result
					result.MismatchesRead = append(result.MismatchesRead, startRead+i)
					numMismatches++

					// skip this match result if there are too many mismatches
					if exceedsMismatchConstraint(read, result) {
						logrus.WithFields(logrus.Fields{
							"mismatchPercentage":    float64(len(result.MismatchesRead)) * 100 / float64(len(*read.Sequence)),
							"maxMismatchPercentage": config.MaxMismatchPercentage(),
							"mismatches":            result.MismatchesRead,
							"numMismatches":         len(result.MismatchesRead),
						}).Debug("too many mismatches in left extension -> skip sequence")
						continue sequenceLoop
					}
				}

				logrus.WithFields(logrus.Fields{
					"startRead":     startRead,
					"endRead":       endRead,
					"startGenome":   startGenome,
					"endGenome":     startGenome + extensionLength,
					"numMismatches": numMismatches,
				}).Debug("extended match to left (linear)")

				result.MatchedRead.AddRegionNonOverlappingPanic(startRead, endRead)
				result.MatchedGenome.AddRegionNonOverlappingPanic(startGenome, startGenome+extensionLength)
			}

			logrus.WithFields(logrus.Fields{
				"read":       result.MatchedRead,
				"genome":     result.MatchedGenome,
				"mismatches": result.MismatchesRead,
			}).Debug("match after left extension")

			// there are unmatched positions in the back of the read
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

				// trying to extend the read to the right
				for i := 0; i < extensionLength; i++ {

					// add the mismatches to the result
					if readSequence[i] != genomeSequence[i] {
						result.MismatchesRead = append(result.MismatchesRead, startRead+i)
					}

					// skip this match result if there are too many mismatches
					if exceedsMismatchConstraint(read, result) {
						logrus.WithFields(logrus.Fields{
							"mismatchPercentage":    float64(len(result.MismatchesRead)) * 100 / float64(len(*read.Sequence)),
							"maxMismatchPercentage": config.MaxMismatchPercentage(),
							"mismatches":            result.MismatchesRead,
							"numMismatches":         len(result.MismatchesRead),
						}).Debug("too many mismatches in right extension -> skip sequence")
						continue sequenceLoop
					}
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

// exceedsMismatchConstraint checks if the number of mismatches exceeds the maximum allowed percentage
// as configured in the config. This value is an integer between 0 and 100 which represents the percentage
// of allowed mismatches in the read.
// The function returns true if the number of mismatches exceeds the allowed percentage.
func exceedsMismatchConstraint(read *fastq.Read, result mapperutils.ReadMatchResult) bool {
	// the percentage of the mismatches accumulated in the given result relative to the read length
	mismatchPercentage := uint8(float64(len(result.MismatchesRead)) * 100 / float64(len(*read.Sequence)))
	return mismatchPercentage > config.MaxMismatchPercentage()
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
	// - maybe keep the readpair for second pass
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
