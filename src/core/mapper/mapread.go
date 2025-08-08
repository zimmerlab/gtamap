package mapper

import (
	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
)

func MapRead(read *fastq.Read, genomeIndex *index.GenomeIndex, greedy bool) ([]*mapperutils.ReadMatchResult, bool) {
	// logrus.WithFields(logrus.Fields{
	// 	"read":   read.Header,
	// 	"length": len(*read.Sequence),
	// }).Debug("Mapping read")

	globalMatches := mapperutils.GlobalMatchResult{
		MatchesPerSequence: make([]*mapperutils.SequenceMatchResult, genomeIndex.NumSequences()*2),
	}

	// fmt.Println(read.Header)

	// TODO: purposely hardcoded the kmer length to 10 because the fixed size array requires manual
	// TODO: intervention when changing the kmer length, having it variable via config could introduce bugs
	// generate non-overlapping k-mers for the read pair
	// skip last kmer if length is not divisible by kmer length
	for i := 0; i <= len(*read.Sequence)-10; i += 10 {

		kmer := [10]byte((*read.Sequence)[i : i+10])

		matches := genomeIndex.FindKeywordMatchesInMap(&kmer, i)

		if matches == nil || len(matches) == 0 {
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
				globalMatches.MatchesPerSequence[match.SequenceIndex].MatchesPerDiagonal[match.StartGenome] = make([]*mapperutils.Match, 0)
			}

			// add match to diagonal
			globalMatches.MatchesPerSequence[match.SequenceIndex].MatchesPerDiagonal[match.StartGenome] = append(globalMatches.MatchesPerSequence[match.SequenceIndex].MatchesPerDiagonal[match.StartGenome], match)
		}
	}

	// list has size of number of sequences in index
	// each element represents the maximum number of kmers that matched exactly on the same diagonal
	maxDiagonalHitsPerSequence := make([]int, genomeIndex.NumSequences()*2)

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

	// logrus.WithFields(logrus.Fields{
	// 	"maxDiagonalHitsPerSequence": maxDiagonalHitsPerSequence,
	// 	"sortedByIndex":              seqIndexSorted,
	// }).Debug("Potential sequences")

	results := make([]*mapperutils.ReadMatchResult, 0)

	// sequenceLoop:
	for _, seqIndex := range seqIndexSorted {

		if seqIndex == -1 {
			// -1 means this sequence was disregarded because of some constraints
			// that can be defined previously to skip certain sequences
			continue
		}

		// logrus.Debug("")
		// logrus.WithFields(logrus.Fields{
		// 	"seqIndex": seqIndex,
		// }).Debug("Considering sequence")

		sequenceMatches := globalMatches.MatchesPerSequence[seqIndex]

		if sequenceMatches == nil {
			// skip this sequence if there are no matches found at all
			continue
		}

		dh := mapperutils.NewDiagonalHandlerWithDataCopy(sequenceMatches.MatchesPerDiagonal)
		initialDepth := 0

		tmpResults := mapReadToSequence(seqIndex, read, genomeIndex, dh, greedy, &initialDepth)
		// fmt.Println(x)

		// var result mapperutils.ReadMatchResult

		//logrus.WithFields(logrus.Fields{
		//	"seqIndex": seqIndex,
		//}).Debug("using sequence to compute best match")

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

		// results = append(results, result)
		results = append(results, tmpResults...)
	}

	if len(results) == 0 {
		// logrus.Debug("no results found")
		return results, false
	}
	return results, true
}

func applyPossibleDiagonals(read *fastq.Read, genomeIndex *index.GenomeIndex, dh *mapperutils.DiagonalHandler,
	result *mapperutils.ReadMatchResult, results *[]*mapperutils.ReadMatchResult, isGreedy bool, currDepth *int,
) {
	// logrus.Debug("")
	// logrus.WithFields(logrus.Fields{
	// 	"seqIndex":      result.SequenceIndex,
	// 	"matchedRead":   result.MatchedRead,
	// 	"matchedGenome": result.MatchedGenome,
	// 	"mismatches":    result.MismatchesRead,
	// 	"length":        result.MatchedRead.Length(),
	// }).Debug("extending read match")

	diagonalsBefore := dh.GetAvailableDiagonals()

	// check if enough kmers can be matched

	// contains the merged region vector of mappable read positions from all remaining diagonals
	// rvRemaining := regionvector.NewRegionVector()
	//
	// for _, d := range diagonalsBefore {
	// 	for _, match := range dh.Diagonals[d] {
	// 		rvRemaining.AddRegionAndMerge(match.FromRead, match.ToRead)
	// 	}
	// }
	//
	// rvUncovered := result.MatchedRead.UncoveredRegionsBySelfAndOther(rvRemaining, 0, len(*read.Sequence))
	//
	// if rvUncovered.Length() >= int(float32(len(*read.Sequence))*0.5) {
	// 	logrus.WithFields(logrus.Fields{
	// 		"uncoveredLength": rvUncovered.Length(),
	// 	}).Debug("uncovered regions are too large")
	// 	return
	// }

	diagonal, score, found := dh.GetBestDiagonal()

	if !found {
		// logrus.Debug("no suitable diagonal found")
		// logrus.Debug("adding partial result to results")

		// check if result is not empty
		if len(result.MatchedGenome.Regions) == 0 || len(result.MatchedRead.Regions) == 0 {
			return
		}

		// INFO: DNA RNA MODE
		// this result should already have a certain length before we append it to results
		if !config.IsOriginRNA {
			// for DNA reads we expect the raw result to already be of a certain length
			if result.MatchedGenome.Length()*10 > len(*read.Sequence)*7 {
				// fmt.Println(read.Header)
				// fmt.Println(result.MatchedGenome.Length())
				*results = append(*results, result)
			}
		} else {
			*results = append(*results, result)
		}

		return
	}

	if len(result.MatchedRead.Regions) == 0 {

		// INFO: DNA RNA MODE
		// here we can expect a larger length for initial diag of map

		l := dh.Diagonals[diagonal][len(dh.Diagonals[diagonal])-1].ToRead - dh.Diagonals[diagonal][0].FromRead
		if config.IsOriginRNA {
			// RNA: check pot length of best initial diag is smaller than 30, if so, don't map read
			if l < 30 {
				return
			}
		} else if l < 50 {
			// DNA: if potential l smaller than 50: return
			return
		}
	}

	// if there are more than 3 remaining diags and our best score < 2, we abort and finish in remap
	// but if we already mapped 3/4 of the read and only lack one more diag, we want to consider it
	// BAD: many remaining diags and best score < 2
	// GOOD/OKAY: few remaining diags and best score < 2 (Happens if we already mapped large portion of reads
	// thus eliminating many of the remaining diags)
	if score < 2 && len(dh.Diagonals) > 3 {
		// logrus.Debug("no suitable diagonal found")
		// logrus.Debug("adding partial result to results")

		// check if result is not empty
		if len(result.MatchedGenome.Regions) == 0 || len(result.MatchedRead.Regions) == 0 {
			return
		}

		// INFO: DNA RNA MODE
		// this result should already have a certain length before we append it to results
		if !config.IsOriginRNA {
			// for DNA reads we expect the raw result to already be of a certain length
			if result.MatchedGenome.Length()*10 > len(*read.Sequence)*7 {
				// fmt.Println(read.Header)
				// fmt.Println(result.MatchedGenome.Length())
				*results = append(*results, result)
			}
		} else {
			*results = append(*results, result)
		}

		return
	}

	// logrus.WithFields(logrus.Fields{
	// 	"diagonal": diagonal,
	// 	"score":    score,
	// }).Debug("found next best diagonal")

	//if score < 2 {
	//	logrus.Debug("no more diagonals with enough matches")
	//	break
	//}

	dhNew := mapperutils.NewDiagonalHandlerWithDataCopy(dh.Diagonals)
	resultNew := result.Copy()

	keepMatch := applyDiagonal(read, genomeIndex, dhNew, diagonal, resultNew)

	if keepMatch && *currDepth <= config.MaxBranchPoints {
		// keep on extending the currently best diagonal
		*currDepth++
		applyPossibleDiagonals(read, genomeIndex, dhNew, resultNew, results, isGreedy, currDepth)
	}

	if isGreedy {
		// exit recursion if greedy
		return
	}

	diagonalsAfter := dhNew.GetAvailableDiagonals()
	excluded := findExcludedDiagonal(diagonalsBefore, diagonalsAfter, diagonal)

	// there is at least one diagonal that was excluded
	wasExcluded := len(excluded) > 0

	if wasExcluded {
		// logrus.WithFields(logrus.Fields{
		// 	"excluded": excluded,
		// }).Debug("excluded diagonals were found")
	} else {
		// logrus.Debug("no diagonals were excluded")
	}

	// the current diagonal contains less than x kmers
	belowScoreThreshold := score <= 1

	// logrus.WithFields(logrus.Fields{
	// 	"belowScoreThreshold": belowScoreThreshold,
	// 	"score":               score,
	// }).Debug("current diagonal score")

	// no condition is met that requires a branching to explore other matches
	if !wasExcluded && !belowScoreThreshold {
		// the currently applied diagonal did not exclude any other diagonals
		return
	}

	// the currently applied diagonal excluded other diagonals
	// remove the currently applied diagonal (without applying) from the diagonal handler
	// and branch to explore other potential mappings

	// logrus.WithFields(logrus.Fields{
	// 	"diagonal": diagonal,
	// }).Debug("exclude the current diagonal from the next search")

	dh.RemoveDiagonal(diagonal)

	if *currDepth <= config.MaxBranchPoints {
		*currDepth++
		applyPossibleDiagonals(read, genomeIndex, dh, result, results, false, currDepth)
	}
}

func findExcludedDiagonal(before []int, after []int, removed int) []int {
	afterSet := make(map[int]struct{})
	for _, a := range after {
		afterSet[a] = struct{}{}
	}

	excluded := make([]int, 0)

	for _, b := range before {
		if b == removed {
			continue
		}
		if _, ok := afterSet[b]; !ok {
			excluded = append(excluded, b)
		}
	}

	return excluded
}

func applyDiagonal(read *fastq.Read, genomeIndex *index.GenomeIndex, dh *mapperutils.DiagonalHandler,
	diagonal int, result *mapperutils.ReadMatchResult,
) bool {
	diagonalRead := regionvector.NewRegionVector()
	diagonalGenome := regionvector.NewRegionVector()

	for _, match := range dh.Diagonals[diagonal] {
		// the region can not be part of another diagonal that is already used
		if match.Used {
			// logrus.Warn("match.Used should not be happening anymore")
			continue
		}
		// add the match to the diagonal
		// diagonalRead.AddRegionNonOverlappingPanic(match.FromRead, match.ToRead)
		// diagonalGenome.AddRegionNonOverlappingPanic(match.FromGenome, match.ToGenome)

		errRead := diagonalRead.AddRegionNonOverlapping(match.FromRead, match.ToRead)
		if errRead != nil {
			logrus.WithFields(logrus.Fields{
				"start": match.FromRead,
				"end":   match.ToRead,
				"rv":    diagonalRead,
				"read":  read.Header,
			}).Fatal("overlapping region in read (apply diagonal match)")
		}

		errGenome := diagonalGenome.AddRegionNonOverlapping(match.FromGenome, match.ToGenome)
		if errGenome != nil {
			logrus.WithFields(logrus.Fields{
				"start": match.FromGenome,
				"end":   match.ToGenome,
				"rv":    diagonalGenome,
				"read":  read.Header,
			}).Fatal("overlapping region in genome (apply diagonal match)")
		}
	}

	// logrus.WithFields(logrus.Fields{
	// 	"posGenome":  diagonal,
	// 	"numMatches": len(dh.Diagonals[diagonal]),
	// 	"read":       diagonalRead,
	// 	"genome":     diagonalGenome,
	// 	//"gaps":       len(diagonalRead.Regions) - 1,
	// }).Debug("applying diagonal")

	gapsRead, gapsGenome := mapperutils.ComputeGapsInDiagonal(diagonalRead, diagonalGenome, result)

	// resolve gaps on the same diagonal
	// gaps can occur because of mismatches in the read and genome which prevent exact kmer matching
	// this is done by filling the gaps in the read and genome and counting the mismatches
	// the regionvectors of read and genome should have the same length as they are coupled
	// because they are part of the same diagonal (no indels, otherwise not on the same diagonal)
	for i := 0; i < len(gapsRead.Regions); i++ {

		gapRead := gapsRead.Regions[i]
		gapGenome := gapsGenome.Regions[i]

		// logrus.WithFields(logrus.Fields{
		// 	"read":   gapRead,
		// 	"genome": gapGenome,
		// }).Debug("found gap")

		// fill the gap in the read by adding the gap as region
		// diagonalRead.AddRegionNonOverlappingPanic(gapRead.Start, gapRead.End)
		errRead := diagonalRead.AddRegionNonOverlapping(gapRead.Start, gapRead.End)
		if errRead != nil {
			logrus.WithFields(logrus.Fields{
				"start": gapRead.Start,
				"end":   gapRead.End,
				"rv":    diagonalRead,
				"read":  read.Header,
			}).Fatal("overlapping region in read (fill diagonal gap)")
		}

		// count mismatches when filling the gap and skip the match if there are too many
		geneSeqPos := 0
		for i := gapRead.Start - 1; i < gapRead.End-1; i++ {

			readByte := (*read.Sequence)[i]
			gIndex := gapGenome.Start + geneSeqPos - 1
			genomeByte := (*genomeIndex.Sequences[result.SequenceIndex])[gIndex]
			geneSeqPos++

			// skip matches
			if readByte == genomeByte {
				continue
			}

			// add the mismatche to the result
			result.MismatchesRead = append(result.MismatchesRead, i)

			// skip this match result if there are too many mismatches
			if exceedsMismatchConstraint(read, result) {
				// logrus.WithFields(logrus.Fields{
				// 	"mismatchPercentage":    float64(len(result.MismatchesRead)) * 100 / float64(len(*read.Sequence)),
				// 	"maxMismatchPercentage": config.MaxMismatchPercentage(),
				// 	"mismatches":            result.MismatchesRead,
				// 	"numMismatches":         len(result.MismatchesRead),
				// }).Debug("too many mismatches in diagonal filling")

				return false
			}
		}

		// diagonalGenome.AddRegionNonOverlappingPanic(gapGenome.Start, gapGenome.End)
		errGenome := diagonalGenome.AddRegionNonOverlapping(gapGenome.Start, gapGenome.End)
		if errGenome != nil {
			logrus.WithFields(logrus.Fields{
				"start": gapGenome.Start,
				"end":   gapGenome.End,
				"rv":    diagonalGenome,
				"read":  read.Header,
			}).Fatal("overlapping region in genome (fill diagonal gap)")
		}
	}

	if len(gapsRead.Regions) > 0 {
		// logrus.WithFields(logrus.Fields{
		// 	"read":       diagonalRead,
		// 	"genome":     diagonalGenome,
		// 	"mismatches": result.MismatchesRead,
		// }).Debug("filled gap")
	}

	// add all matches in diagonal to the result
	for i := 0; i < len(diagonalRead.Regions); i++ {
		regionRead := diagonalRead.Regions[i]
		regionGenome := diagonalGenome.Regions[i]

		// logrus.WithFields(logrus.Fields{
		// 	"read":   regionRead,
		// 	"genome": regionGenome,
		// }).Debug("adding diagonal match to result")

		// result.MatchedRead.AddRegionNonOverlappingPanic(regionRead.Start, regionRead.End)
		// result.MatchedGenome.AddRegionNonOverlappingPanic(regionGenome.Start, regionGenome.End)

		errRead := result.MatchedRead.AddRegionNonOverlapping(regionRead.Start, regionRead.End)
		if errRead != nil {
			logrus.WithFields(logrus.Fields{
				"start": regionRead.Start,
				"end":   regionRead.End,
				"rv":    result.MatchedRead,
				"read":  read.Header,
			}).Fatal("overlapping region in read (finalize diagonal match)")
		}

		errGenome := result.MatchedGenome.AddRegionNonOverlapping(regionGenome.Start, regionGenome.End)
		if errGenome != nil {
			logrus.WithFields(logrus.Fields{
				"start": regionGenome.Start,
				"end":   regionGenome.End,
				"rv":    result.MatchedGenome,
				"read":  read.Header,
			}).Fatal("overlapping region in genome (finalize diagonal match)")
		}
	}

	// mark all used regions
	for i := 0; i < len(diagonalRead.Regions); i++ {
		dh.ConsumeKmer(diagonalRead.Regions[i].Start, diagonalRead.Regions[i].End,
			diagonalGenome.Regions[i].Start, diagonalGenome.Regions[i].End)
	}

	// remove used matches and empty diagonals
	dh.RemovedConsumedRegionsAndDiagonals()

	// remove invalid diagonals based on the currently applied diagonal
	dh.RemoveInvalidDiagonals(result, read)

	return true
}

func annotateSpliceSites(read *fastq.Read, genomeIndex *index.GenomeIndex, result *mapperutils.ReadMatchResult) {
	// used to keep track of the read position for the next gap
	readGapPos := 0
	// returns the index of the first region after which a gap occurs (-1 if no gap)
	indexRegionBeforeGap := result.MatchedGenome.GetGapIndexAfterPos(readGapPos)

	// loop through all gaps in the read (-1 means there is no more gap)
	seqIndex := result.SequenceIndex
	for indexRegionBeforeGap > -1 {

		gapGenome, _ := result.MatchedGenome.GetGapAfterRegionIndex(indexRegionBeforeGap)

		donorSiteStart := gapGenome.Start
		donorSiteSeq := (*genomeIndex.Sequences[seqIndex])[donorSiteStart : donorSiteStart+2]

		acceptorSiteStart := gapGenome.End
		acceptorSiteSeq := (*genomeIndex.Sequences[seqIndex])[acceptorSiteStart-2 : acceptorSiteStart]

		var lookOnPlusStrand bool
		if genomeIndex.GetSequenceInfo(seqIndex / 2).IsForwardStrand {
			if genomeIndex.IsSequenceForward(seqIndex) {
				lookOnPlusStrand = true
			} else {
				lookOnPlusStrand = false
			}
		} else {
			if genomeIndex.IsSequenceForward(seqIndex) {
				lookOnPlusStrand = false
			} else {
				lookOnPlusStrand = true
			}
		}

		_, isKnownSpliceSite := utils.ScoreSpliceSites(donorSiteSeq[0], donorSiteSeq[1],
			acceptorSiteSeq[0], acceptorSiteSeq[1], lookOnPlusStrand)

		// annotate split with spliceSite info
		if result.SpliceSitesInfo == nil {
			result.SpliceSitesInfo = make([]bool, 0)
			result.SpliceSitesInfo = append(result.SpliceSitesInfo, isKnownSpliceSite)
		} else {
			result.SpliceSitesInfo = append(result.SpliceSitesInfo, isKnownSpliceSite)
		}

		readGapPos = gapGenome.End + 1
		indexRegionBeforeGap = result.MatchedGenome.GetGapIndexAfterPos(readGapPos)
	}
}

func extendDiagonals(read *fastq.Read, genomeIndex *index.GenomeIndex, result *mapperutils.ReadMatchResult) {
	if result.MatchedRead.Length() == len(*read.Sequence) {
		// logrus.Debug("read fully matched")
	} else {

		// there are gaps (unmatched regions) in the read (between matched regions / diagonals)
		if result.MatchedRead.HasGaps() {

			// logrus.Debug("unmatched regions within read")

			// used to keep track of the read position for the next gap
			readGapPos := 0
			// returns the index of the first region after which a gap occurs (-1 if no gap)
			indexRegionBeforeGap := result.MatchedRead.GetGapIndexAfterPos(readGapPos)

			// loop through all gaps in the read (-1 means there is no more gap)
			for indexRegionBeforeGap > -1 {

				gapRead, _ := result.MatchedRead.GetGapAfterRegionIndex(indexRegionBeforeGap)
				gapGenome, gapGenomeOk := result.MatchedGenome.GetGapAfterRegionIndex(indexRegionBeforeGap)

				if !gapGenomeOk {
					// logrus.WithFields(logrus.Fields{
					// 	"read":    result.MatchedRead,
					// 	"genome":  result.MatchedGenome,
					// 	"gapRead": gapRead,
					// }).Debug("insertion found")
				} else {

					// logrus.WithFields(logrus.Fields{
					// 	"gapRead":   gapRead,
					// 	"gapGenome": gapGenome,
					// }).Debug("found gap to be handled")

					// This case is covered in remap read (fillGaps)
					if gapRead.Length() > gapGenome.Length() {
						// logrus.WithFields(logrus.Fields{
						// 	"gapRead":   gapRead,
						// 	"gapGenome": gapGenome,
						// }).Debug("gap read is larger than gap genome")

						result.IncompleteMap = true
						return
					}

					bestSplit := determineBestSplit(genomeIndex, read, result.SequenceIndex, gapRead, gapGenome)

					if bestSplit == -1 {
						// this should not happen because a split should be found every time
						// even if the split has a bad score (many mismatches, no splice sites, etc)
						logrus.WithFields(logrus.Fields{
							"qname": read.Header,
						}).Fatal("no best split found")
					}

					// logrus.WithFields(logrus.Fields{
					// 	"split": bestSplit,
					// }).Debug("best split found")

					// logrus.WithFields(logrus.Fields{
					// 	"read":   result.MatchedRead,
					// 	"genome": result.MatchedGenome,
					// }).Debug("regions before")

					// when bestSplit is 0 then there is nothing to be added to the left side of the gap
					if bestSplit > 0 {

						// add the split to the result
						result.MatchedRead.AddRegionNonOverlappingPanic(gapRead.Start, gapRead.Start+bestSplit)
						result.MatchedGenome.AddRegionNonOverlappingPanic(gapGenome.Start, gapGenome.Start+bestSplit)

						// the read and genome sequences from the start of the gap to the best split (left)
						readByte := (*read.Sequence)[gapRead.Start : gapRead.Start+bestSplit]
						genomeByte := (*genomeIndex.Sequences[result.SequenceIndex])[gapGenome.Start : gapGenome.Start+bestSplit]

						// add the mismatches to the result
						for i := 0; i < bestSplit; i++ {
							// add the mismatche to the result
							if readByte[i] != genomeByte[i] {
								result.MismatchesRead = append(result.MismatchesRead, gapRead.Start+i)

								// skip this match result if there are too many mismatches
								if exceedsMismatchConstraint(read, result) {
									// logrus.WithFields(logrus.Fields{
									// 	"mismatchPercentage":    float64(len(result.MismatchesRead)) * 100 / float64(len(*read.Sequence)),
									// 	"maxMismatchPercentage": config.MaxMismatchPercentage(),
									// 	"mismatches":            result.MismatchesRead,
									// 	"numMismatches":         len(result.MismatchesRead),
									// }).Debug("too many mismatches in middle extension (left) -> skip sequence")
									// continue sequenceLoop

									result.IncompleteMap = true
									return
								}
							}
						}
					}

					// logrus.WithFields(logrus.Fields{
					// 	"read":   result.MatchedRead,
					// 	"genome": result.MatchedGenome,
					// }).Debug("regions after left")

					// when bestSplit is equal to the length of the gap then there is nothing
					// to be added to the right side of the gap
					if bestSplit < gapRead.Length() {

						// logrus.WithFields(logrus.Fields{
						// 	"gapReadEnd":    gapRead.End,
						// 	"bestSplit":     bestSplit,
						// 	"gapReadLength": gapRead.Length(),
						// 	"left":          gapRead.End - (gapRead.Length() - bestSplit),
						// 	"right":         gapRead.End,
						// }).Debug("debug split right")

						// add the split to the result
						result.MatchedRead.AddRegionNonOverlappingPanic(gapRead.End-(gapRead.Length()-bestSplit), gapRead.End)
						result.MatchedGenome.AddRegionNonOverlappingPanic(gapGenome.End-(gapRead.Length()-bestSplit), gapGenome.End)

						// the read and genome sequences from the best split to the end of the gap (right)
						readByte := (*read.Sequence)[gapRead.End-(gapRead.Length()-bestSplit) : gapRead.End]
						genomeByte := (*genomeIndex.Sequences[result.SequenceIndex])[gapGenome.End-(gapRead.Length()-bestSplit) : gapGenome.End]

						// add the mismatches to the result
						for i := 0; i < gapRead.Length()-bestSplit; i++ {
							// add the mismatche to the result
							if readByte[i] != genomeByte[i] {
								// result.MismatchesRead = append(result.MismatchesRead, gapRead.End-(bestSplit-i))
								result.MismatchesRead = append(result.MismatchesRead, gapRead.Start+bestSplit+i)

								// skip this match result if there are too many mismatches
								if exceedsMismatchConstraint(read, result) {
									// logrus.WithFields(logrus.Fields{
									// 	"mismatchPercentage":    float64(len(result.MismatchesRead)) * 100 / float64(len(*read.Sequence)),
									// 	"maxMismatchPercentage": config.MaxMismatchPercentage(),
									// 	"mismatches":            result.MismatchesRead,
									// 	"numMismatches":         len(result.MismatchesRead),
									// }).Debug("too many mismatches in middle extension (right) -> skip sequence")

									result.IncompleteMap = true
									return
								}
							}
						}
					}

					// logrus.WithFields(logrus.Fields{
					// 	"read":   result.MatchedRead,
					// 	"genome": result.MatchedGenome,
					// }).Debug("regions after right")
				}

				// determine the next gap (-1 if there is none)
				readGapPos = gapRead.End + 1
				indexRegionBeforeGap = result.MatchedRead.GetGapIndexAfterPos(readGapPos)
			}
		}

		// there are unmatched positions in front of the read
		firstRegionRead, _ := result.MatchedRead.GetFirstRegion()
		if firstRegionRead.Start > 0 {

			startRead := 0
			endRead := firstRegionRead.Start
			extensionLength := endRead - startRead

			firstRegionGenome, _ := result.MatchedGenome.GetFirstRegion()
			startGenome := firstRegionGenome.Start - extensionLength

			readSequence := (*read.Sequence)[startRead:endRead]

			// logrus.WithFields(logrus.Fields{
			// 	"startRead":   startRead,
			// 	"endRead":     endRead,
			// 	"startGenome": startGenome,
			// 	"endGenome":   startGenome + len(readSequence),
			// }).Debug("unmatched positions in front of read")

			// TODO: this could be a use case for clipping
			// when a read maps to the target sequence but would go out of bounds
			// it could still be a valid mapping but it needs to be clipped
			if startGenome < 0 {
				// logrus.Debug("genome index out of bounds")

				result.IncompleteMap = true
				return
			}

			genomeSequence := (*genomeIndex.Sequences[result.SequenceIndex])[startGenome : startGenome+len(readSequence)]

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
					// logrus.WithFields(logrus.Fields{
					// 	"mismatchPercentage":    float64(len(result.MismatchesRead)) * 100 / float64(len(*read.Sequence)),
					// 	"maxMismatchPercentage": config.MaxMismatchPercentage(),
					// 	"mismatches":            result.MismatchesRead,
					// 	"numMismatches":         len(result.MismatchesRead),
					// }).Debug("too many mismatches in left extension -> skip sequence")

					result.IncompleteMap = true
					return
				}
			}

			// logrus.WithFields(logrus.Fields{
			// 	"startRead":     startRead,
			// 	"endRead":       endRead,
			// 	"startGenome":   startGenome,
			// 	"endGenome":     startGenome + extensionLength,
			// 	"numMismatches": numMismatches,
			// }).Debug("extended match to left (linear)")

			result.MatchedRead.AddRegionNonOverlappingPanic(startRead, endRead)
			result.MatchedGenome.AddRegionNonOverlappingPanic(startGenome, startGenome+extensionLength)
		}

		// logrus.WithFields(logrus.Fields{
		// 	"read":       result.MatchedRead,
		// 	"genome":     result.MatchedGenome,
		// 	"mismatches": result.MismatchesRead,
		// }).Debug("match after left extension")

		// there are unmatched positions in the back of the read
		lastRegionRead, _ := result.MatchedRead.GetLastRegion()
		if lastRegionRead.End < len(*read.Sequence) {

			startRead := lastRegionRead.End
			endRead := len(*read.Sequence)
			extensionLength := endRead - startRead

			lastRegionGenome, _ := result.MatchedGenome.GetLastRegion()
			startGenome := lastRegionGenome.End

			readSequence := (*read.Sequence)[startRead:endRead]

			// logrus.WithFields(logrus.Fields{
			// 	"startRead":   startRead,
			// 	"endRead":     endRead,
			// 	"startGenome": startGenome,
			// 	"endGenome":   startGenome + len(readSequence),
			// }).Debug("unmatched positions in back of read")

			// TODO: this could be a use case for clipping
			// when a read maps to the target sequence but would go out of bounds
			// it could still be a valid mapping but it needs to be clipped
			if startGenome+len(readSequence) > len(*genomeIndex.Sequences[result.SequenceIndex]) {
				// logrus.Debug("genome index out of bounds")

				result.IncompleteMap = true
				return
			}

			genomeSequence := (*genomeIndex.Sequences[result.SequenceIndex])[startGenome : startGenome+len(readSequence)]

			// trying to extend the read to the right
			for i := 0; i < extensionLength; i++ {

				// add the mismatches to the result
				if readSequence[i] != genomeSequence[i] {
					result.MismatchesRead = append(result.MismatchesRead, startRead+i)
				}

				// skip this match result if there are too many mismatches
				if exceedsMismatchConstraint(read, result) {
					// logrus.WithFields(logrus.Fields{
					// 	"mismatchPercentage":    float64(len(result.MismatchesRead)) * 100 / float64(len(*read.Sequence)),
					// 	"maxMismatchPercentage": config.MaxMismatchPercentage(),
					// 	"mismatches":            result.MismatchesRead,
					// 	"numMismatches":         len(result.MismatchesRead),
					// }).Debug("too many mismatches in right extension -> skip sequence")
					// continue sequenceLoop

					result.IncompleteMap = true
					return
				}
			}

			result.MatchedRead.AddRegionNonOverlappingPanic(startRead, endRead)
			result.MatchedGenome.AddRegionNonOverlappingPanic(startGenome, startGenome+extensionLength)
		}
	}
}

func mapReadToSequence(seqIndex int, read *fastq.Read, genomeIndex *index.GenomeIndex,
	diagonalHandler *mapperutils.DiagonalHandler, greedy bool, currDepth *int,
) []*mapperutils.ReadMatchResult {
	// list of read match results
	results := make([]*mapperutils.ReadMatchResult, 0)

	result := &mapperutils.ReadMatchResult{
		SequenceIndex:  seqIndex,
		MatchedRead:    regionvector.NewRegionVector(),
		MatchedGenome:  regionvector.NewRegionVector(),
		MismatchesRead: make([]int, 0),
	}

	applyPossibleDiagonals(read, genomeIndex, diagonalHandler, result, &results, greedy, currDepth)

	finalResults := make([]*mapperutils.ReadMatchResult, 0)

	for _, res := range results {

		extendDiagonals(read, genomeIndex, res)

		// TODO: DNA RNA
		// Only annotate if RNA
		if res.MatchedGenome.HasGaps() {
			annotateSpliceSites(read, genomeIndex, res)
		}

		finalResults = append(finalResults, res)
	}

	return finalResults
}

// exceedsMismatchConstraint checks if the number of mismatches exceeds the maximum allowed percentage
// as configured in the config. This value is an integer between 0 and 100 which represents the percentage
// of allowed mismatches in the read.
// The function returns true if the number of mismatches exceeds the allowed percentage.
func exceedsMismatchConstraint(read *fastq.Read, result *mapperutils.ReadMatchResult) bool {
	// the percentage of the mismatches accumulated in the given result relative to the read length
	mismatchPercentage := uint8(float64(len(result.MismatchesRead)) * 100 / float64(len(*read.Sequence)))
	return mismatchPercentage > config.MaxMismatchPercentage()
}

func determineBestSplit(
	genomeIndex *index.GenomeIndex,
	read *fastq.Read,
	seqIndex int,
	gapRead regionvector.Region,
	gapGenome regionvector.Region,
) int {
	// logrus.WithFields(logrus.Fields{
	// 	"gapRead":   gapRead,
	// 	"gapGenome": gapGenome,
	// }).Debug("determining best split")

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

	// logrus.WithFields(logrus.Fields{
	// 	"lErrors": lErrors,
	// 	"rErrors": rErrors,
	// }).Debug("determined mismatches")

	// the minimum number of mismatches
	// the +2 is based on the maximum penalty returned by scoreSpliceSites()
	minErrors := lErrors[gapRead.Length()] + rErrors[gapRead.Length()] + 2
	// the position of the split with the minimum number of mismatches
	minSplit := -1

	// TODO: keep track of the actual mismatch positions
	// TODO: if no suitable split is found then:
	// - maybe there is another exon in between if enough bases missing from read
	// - maybe keep the readpair for unmapped pass
	for i := 0; i <= gapRead.Length(); i++ {

		lPos := i
		rPos := gapRead.Length() - i

		numMismatches := lErrors[lPos] + rErrors[rPos]

		donorSiteStart := gapGenome.Start + i
		donorSiteSeq := (*genomeIndex.Sequences[seqIndex])[donorSiteStart : donorSiteStart+2]

		splitRev := gapRead.Length() - i
		acceptorSiteStart := gapGenome.End - splitRev
		acceptorSiteSeq := (*genomeIndex.Sequences[seqIndex])[acceptorSiteStart-2 : acceptorSiteStart]

		var lookOnPlusStrand bool
		if genomeIndex.GetSequenceInfo(seqIndex / 2).IsForwardStrand {
			if genomeIndex.IsSequenceForward(seqIndex) {
				lookOnPlusStrand = true
			} else {
				lookOnPlusStrand = false
			}
		} else {
			if genomeIndex.IsSequenceForward(seqIndex) {
				lookOnPlusStrand = false
			} else {
				lookOnPlusStrand = true
			}
		}

		// add a penalty if the splice site is not canonical
		// 2 means that there is no known splice site

		// INFO: DNA RNA MODE
		// only score splicesites in RNA mode

		spliceSitePenalty, _ := utils.ScoreSpliceSites(donorSiteSeq[0], donorSiteSeq[1],
			acceptorSiteSeq[0], acceptorSiteSeq[1], lookOnPlusStrand)
		numMismatches += spliceSitePenalty

		// logrus.WithFields(logrus.Fields{
		// 	"split":               i,
		// 	"splice site penalty": spliceSitePenalty,
		// 	"numMismatches":       numMismatches,
		// }).Debug("possible split")

		if numMismatches < minErrors {
			minErrors = numMismatches
			minSplit = i
		}
	}

	return minSplit
}
