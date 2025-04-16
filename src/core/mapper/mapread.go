package mapper

import (
	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/sirupsen/logrus"
	"os"
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

	//for seqIndex, sequenceMatches := range globalMatches.MatchesPerSequence {
	//
	//	if sequenceMatches == nil {
	//		continue
	//	}
	//
	//	dSizes := make([]int, len(sequenceMatches.MatchesPerDiagonal))
	//	dIndices := make([]int, len(sequenceMatches.MatchesPerDiagonal))
	//	counter := 0
	//
	//	for diagonal, matches := range sequenceMatches.MatchesPerDiagonal {
	//		fmt.Println(seqIndex, diagonal, len(matches))
	//
	//		dSizes[counter] = len(matches)
	//		dIndices[counter] = diagonal
	//		counter++
	//	}
	//
	//	fmt.Println(dSizes)
	//	fmt.Println(dIndices)
	//}

	results := make([]mapperutils.ReadMatchResult, 0)

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
			SecondPass:     false,
		}

		diagonalHandler := mapperutils.NewDiagonalHandlerWithData(sequenceMatches.MatchesPerDiagonal)

		for result.MatchedRead.Length() < len(*read.Sequence) {

			logrus.Debug("searching for next best diagonal")

			bestDiagonal, bestDiagonalLength := diagonalHandler.GetBestDiagonal()

			if bestDiagonal == -1 {
				logrus.Debug("no suitable diagonal found")
				break
			}

			diagonalRead := regionvector.NewRegionVector()
			diagonalGenome := regionvector.NewRegionVector()

			for _, match := range sequenceMatches.MatchesPerDiagonal[bestDiagonal] {
				// the region can not be part of another diagonal that is already used
				if match.Used {
					continue
				}
				diagonalRead.AddRegion(match.FromRead, match.ToRead)
				diagonalGenome.AddRegion(match.FromGenome, match.ToGenome)
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

			foundGap := len(diagonalRead.Regions) > 1

			// more than one element means there is a gap
			// the regionvectors of read and genome should have the same length as they are coupled
			for len(diagonalRead.Regions) > 1 {

				gapRead := diagonalRead.GetFirstGap()
				gapGenome := diagonalGenome.GetFirstGap()

				logrus.WithFields(logrus.Fields{
					"read":   gapRead,
					"genome": gapGenome,
				}).Debug("found gap")

				diagonalRead.AddRegion(gapRead.Start, gapRead.End)
				for i := gapRead.Start; i < gapRead.End; i++ {

					readByte := (*read.Sequence)[i]
					gIndex := diagonalGenome.GetFirstRegion().Start + i
					genomeByte := (*genomeIndex.Sequences[seqIndex])[gIndex]

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
				diagonalGenome.AddRegion(gapGenome.Start, gapGenome.End)
			}

			if foundGap {
				logrus.WithFields(logrus.Fields{
					"read":       diagonalRead,
					"genome":     diagonalGenome,
					"mismatches": mismatches,
				}).Debug("filled gap")
			}

			// add digonal to result (diagonal rv should have length 1 after gap filling)
			result.MatchedRead.AddRegion(diagonalRead.GetFirstRegion().Start, diagonalRead.GetFirstRegion().End)
			result.MatchedGenome.AddRegion(diagonalGenome.GetFirstRegion().Start, diagonalGenome.GetFirstRegion().End)
			result.MismatchesRead = append(result.MismatchesRead, mismatches...)

			// remove diagonal from handler
			diagonalHandler.ConsumeDiagonal(bestDiagonal)
		}

		logrus.WithFields(logrus.Fields{
			"read":   result.MatchedRead,
			"genome": result.MatchedGenome,
		}).Debug("done processing diagonals")

		if result.MatchedRead.Length() == len(*read.Sequence) {
			logrus.Debug("read fully matched")
		} else {

			if len(result.MatchedRead.Regions) > 1 {

				logrus.Debug("unmatched regions within read")

				numGaps := len(result.MatchedRead.Regions) - 1
				genomeGapIndex := 0

				// resolve the borders within the read where at least one junction was found
				for len(result.MatchedRead.Regions) > 1 {

					// position in read where the gap starts
					lRead := result.MatchedRead.GetFirstGap().Start
					// position in genome where the gap starts
					lGenome := result.MatchedGenome.GetGap(genomeGapIndex).Start
					// position in read where the gap ends
					rRead := result.MatchedRead.GetFirstGap().End
					// position in genome where the gap ends
					rGenome := result.MatchedGenome.GetGap(genomeGapIndex).End

					logrus.WithFields(logrus.Fields{
						"startRead":   lRead,
						"endRead":     rRead,
						"startGenome": lGenome,
						"endGenome":   rGenome,
					}).Debug("resolving gap within read")

					//fmt.Println("handling gap:")
					//fmt.Println("lRead: ", lRead)
					//fmt.Println("rRead: ", rRead)
					//fmt.Println("lGenome: ", lGenome)
					//fmt.Println("rGenome: ", rGenome)

					//readSequence := (*read.Sequence)[lRead:rRead]
					//genomeSequenceL := (*genomeIndex.Sequences[seqIndex])[lGenome : lGenome+15]
					//genomeSequenceR := (*genomeIndex.Sequences[seqIndex])[rGenome-15 : rGenome]
					//fmt.Println("READ:\t\t", string(readSequence))
					//fmt.Println("GENOME L:\t", string(genomeSequenceL))
					//fmt.Println("GENOME R:\t", string(genomeSequenceR))
					//
					//fmt.Println("READ (rev):\t", string(utils.ReverseComplementDnaBytes(readSequence)))
					//fmt.Println("GENOME L (rev):\t", string(utils.ReverseComplementDnaBytes(genomeSequenceR)))
					//fmt.Println("GENOME R (rev):\t", string(utils.ReverseComplementDnaBytes(genomeSequenceL)))
					//fmt.Println(lRead, lGenome, rRead, rGenome)

					// the length of the extension based on the number if unmapped bases in the read
					extensionLength := rRead - lRead

					//fmt.Println("extension length: ", extensionLength)

					// cululative mismatch count for the left and right extensions
					// for lErrors the index i represents the number of mismatches for the first i positions of the extension
					// for rErrors the index i represents the number of mismatches for the last i positions of the extension
					lErrors := make([]int, extensionLength+1)
					rErrors := make([]int, extensionLength+1)

					lErrors[0] = 0
					rErrors[0] = 0

					for i := 1; i <= extensionLength; i++ {
						lErrors[i] = lErrors[i-1]
						if (*read.Sequence)[lRead+i-1] != (*genomeIndex.Sequences[seqIndex])[lGenome+i-1] {
							lErrors[i]++
						}

						rErrors[i] = rErrors[i-1]
						if (*read.Sequence)[rRead-i] != (*genomeIndex.Sequences[seqIndex])[rGenome-i] {
							rErrors[i]++
						}
					}

					//fmt.Println("lErrors:", lErrors)
					//fmt.Println("rErrors:", rErrors)

					// the minimum number of mismatches
					minErrors := lErrors[extensionLength] + rErrors[extensionLength]
					// the position of the split with the minimum number of mismatches
					minSplit := -1

					// TODO: keep track of the actual mismatch positions
					// TODO: if no suitable split is found then:
					// - maybe there is another exon in between if enough bases missing from read
					// - maybe keep the readpair for second pass
					for i := 0; i <= extensionLength; i++ {

						lPos := i
						rPos := extensionLength - i

						numMismatches := lErrors[lPos] + rErrors[rPos]

						donorSiteStart := result.MatchedGenome.GetGap(genomeGapIndex).Start + i
						donorSiteSeq := (*genomeIndex.Sequences[seqIndex])[donorSiteStart : donorSiteStart+2]

						splitRev := extensionLength - i
						acceptorSiteStart := result.MatchedGenome.GetGap(genomeGapIndex).End - splitRev
						acceptorSiteSeq := (*genomeIndex.Sequences[seqIndex])[acceptorSiteStart-2 : acceptorSiteStart]

						isForwardStrand := genomeIndex.IsSequenceForward(seqIndex)

						// add a penalty if the splice site is not canonical
						// 2 means that there is no known splice site
						numMismatches += scoreSpliceSites(donorSiteSeq[0], donorSiteSeq[1],
							acceptorSiteSeq[0], acceptorSiteSeq[1], isForwardStrand)

						if numMismatches < minErrors {
							minErrors = numMismatches
							minSplit = i
						}
					}

					// apply the best split
					if minSplit == -1 {
						logrus.Error("no best split found!")
						os.Exit(1)
					}

					//fmt.Println("chose: ", minSplit)
					//fmt.Println("read before", result.MatchedRead.String())
					//fmt.Println("genome before", result.MatchedGenome.String())
					//fmt.Println("add l read: ", lRead, lRead+minSplit)
					//fmt.Println("add l genome: ", lGenome, lGenome+minSplit)

					logrus.WithFields(logrus.Fields{
						"read":   result.MatchedRead,
						"genome": result.MatchedGenome,
					}).Debug("regions before")

					result.MatchedRead.AddRegion(lRead, lRead+minSplit)
					result.MatchedGenome.AddRegion(lGenome, lGenome+minSplit)

					logrus.WithFields(logrus.Fields{
						"read":   result.MatchedRead,
						"genome": result.MatchedGenome,
					}).Debug("regions after left")

					//fmt.Println("read after add l: ", result.MatchedRead.String())
					//fmt.Println("genome after add l: ", result.MatchedGenome.String())
					//fmt.Println("add r read: ", rRead-(extensionLength-minSplit), rRead)
					//fmt.Println("add r genome: ", rGenome-(extensionLength-minSplit), rGenome)

					result.MatchedRead.AddRegion(rRead-(extensionLength-minSplit), rRead)
					result.MatchedGenome.AddRegion(rGenome-(extensionLength-minSplit), rGenome)

					logrus.WithFields(logrus.Fields{
						"read":   result.MatchedRead,
						"genome": result.MatchedGenome,
					}).Debug("regions after right")

					//fmt.Println("read after add r: ", result.MatchedRead.String())
					//fmt.Println("genome after add r: ", result.MatchedGenome.String())

					// only increment the gap index if no gap was removed during the split
					if numGaps == len(result.MatchedGenome.Regions)-1 {
						genomeGapIndex++
					}
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

				// add to second pass if there are too many mismatches
				if numMismatches > 5 {
					result.SecondPass = true
					results = append(results, result)
					continue sequenceLoop
				}

				result.MatchedRead.AddRegion(startRead, endRead)
				result.MatchedGenome.AddRegion(startGenome, startGenome+extensionLength)
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

				// add to second pass if there are too many mismatches
				if numMismatches > 5 {
					result.SecondPass = true
					results = append(results, result)
					continue sequenceLoop
				}

				result.MatchedRead.AddRegion(startRead, endRead)
				result.MatchedGenome.AddRegion(startGenome, startGenome+extensionLength)
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
