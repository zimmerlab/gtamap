package secondpass

import (
	"fmt"
	"sort"
	"strings"
	"sync"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/datastructure"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/core/mapper/thirdpass"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
)

func SecondpassMappingWorker(secondPassChan *SecondPassChannel, wgIncompleteMapping *sync.WaitGroup, annotationChan <-chan map[int]*mapperutils.TargetAnnotation, thirdPassChan *thirdpass.ThirdPassChannel, genomeIndex *index.GenomeIndex) {
	// in here we receive all non-confident readpairs. Some have multiple maps for fw and rv and some are
	// not completely mapped yet. Confident readpairs are contained in confidentChan.
	defer wgIncompleteMapping.Done()

	// main seq id -> seq annot
	var annotation map[int]*mapperutils.TargetAnnotation
	for annot := range annotationChan {
		// is only one object
		annotation = annot
	}
	logrus.Info("Started second pass")
	logrus.Info("Started output worker")

	var wgRemap sync.WaitGroup

	for {
		task, ok := secondPassChan.Receive()

		if !ok {
			break
		}

		wgRemap.Add(1)
		go func(t *mapperutils.ReadPairMatchResults) {
			defer wgRemap.Done()

			remapReadPair(task, annotation, genomeIndex)

			thirdPassChan.Send(&thirdpass.ThirdPassTask{
				ReadPairId: task.ReadPair.ReadR1.Header,
				TargetInfo: task,
			})
		}(task)
	}
	wgRemap.Wait()
}

func remapReadPair(readPairMapping *mapperutils.ReadPairMatchResults, annotationMap map[int]*mapperutils.TargetAnnotation, genomeIndex *index.GenomeIndex) {
	fwRemaps := make([]*mapperutils.ReadMatchResult, 0)

	for _, mapping := range readPairMapping.Fw {
		mainSeqId := mapping.SequenceIndex / 2
		remaps := remapRead(mapping, annotationMap[mainSeqId], readPairMapping.ReadPair.ReadR1, genomeIndex)
		fwRemaps = append(fwRemaps, remaps...)
	}

	rvRemaps := make([]*mapperutils.ReadMatchResult, 0)

	for _, mapping := range readPairMapping.Rv {
		mainSeqId := mapping.SequenceIndex / 2
		remaps := remapRead(mapping, annotationMap[mainSeqId], readPairMapping.ReadPair.ReadR2, genomeIndex)
		rvRemaps = append(rvRemaps, remaps...)
	}

	// make mappings uniq
	uniqFwRemaps := getUniqRemaps(fwRemaps)
	uniqRvRemaps := getUniqRemaps(rvRemaps)

	// get valid mappings
	if len(uniqFwRemaps) > 0 {
		validMaps := filterValidMaps(uniqFwRemaps, len(*readPairMapping.ReadPair.ReadR1.Sequence))
		readPairMapping.Fw = validMaps
	} else {
		validMaps := filterValidMaps(readPairMapping.Fw, len(*readPairMapping.ReadPair.ReadR1.Sequence))
		readPairMapping.Fw = validMaps
	}

	// get valid mappings
	if len(uniqRvRemaps) > 0 {
		validMaps := filterValidMaps(uniqRvRemaps, len(*readPairMapping.ReadPair.ReadR2.Sequence))
		readPairMapping.Rv = validMaps
	} else {
		validMaps := filterValidMaps(readPairMapping.Rv, len(*readPairMapping.ReadPair.ReadR2.Sequence))
		readPairMapping.Rv = validMaps
	}
}

func getUniqRemaps(r []*mapperutils.ReadMatchResult) []*mapperutils.ReadMatchResult {
	uniq := make([]*mapperutils.ReadMatchResult, 0)
	seen := make(map[string]bool)
	for _, remap := range r {
		remap.MergeRegions()
		hash := createHash(remap)
		if !seen[hash] {
			seen[hash] = true
			uniq = append(uniq, remap)
		}
	}
	return uniq
}

func filterValidMaps(mappings []*mapperutils.ReadMatchResult, readLength int) []*mapperutils.ReadMatchResult {
	valid := make([]*mapperutils.ReadMatchResult, 0)
	for _, mapping := range mappings {
		if mapping.MatchedRead.Regions[0].Start == 0 && mapping.MatchedRead.Regions[len(mapping.MatchedRead.Regions)-1].End == readLength && mapping.MatchedGenome.Length() == readLength {
			// only append complete remaps (remap can still be shorter than read length but only if insert)
			if uint8(float64(len(mapping.MismatchesRead))*100/float64(readLength)) <= config.MaxMismatchPercentage() {
				// last check to make sure we do not append remaps which exeed max mm prec (this is required again because fillGaps potentially adds more mm in a remap)
				valid = append(valid, mapping)
			}
		}
	}
	// if len(valid) == 0 {
	// 	return mappings, false
	// } else {
	// 	return valid, true
	// }
	return valid
}

func createHash(r *mapperutils.ReadMatchResult) string {
	var sb strings.Builder
	for _, iv := range r.MatchedGenome.Regions {
		sb.WriteString(fmt.Sprintf("%d,", iv.Start))
	}
	sb.WriteString("|")
	for _, iv := range r.MatchedRead.Regions {
		sb.WriteString(fmt.Sprintf("%d,", iv.Start))
	}
	return sb.String()
}

func remapRead(readMapping *mapperutils.ReadMatchResult, annotation *mapperutils.TargetAnnotation, read *fastq.Read, genomeIndex *index.GenomeIndex) []*mapperutils.ReadMatchResult {
	alternativeReadMatchResults := make([]*mapperutils.ReadMatchResult, 0)
	readMapping.NormalizeRegions() // NOTE: this is CRUCIAL and NEEDS to be called before ANY REMAP!!!

	// do we have an annotation?
	if annotation != nil {
		if readMapping.IncompleteMap {
			remaps := fixPointRNARemap(readMapping, annotation.Introns[readMapping.SequenceIndex], read, genomeIndex, annotation.IntronTrees[readMapping.SequenceIndex])
			alternativeReadMatchResults = append(alternativeReadMatchResults, remaps...)
		} else {
			overhangCorrected := correctOverhangs(readMapping, annotation.Introns[readMapping.SequenceIndex], read, genomeIndex, annotation.IntronTrees[readMapping.SequenceIndex])
			alternativeReadMatchResults = append(alternativeReadMatchResults, overhangCorrected...)
		}
	}

	// fill potential gaps in all remaps
	for _, alternativeMap := range alternativeReadMatchResults {
		if alternativeMap.MatchedRead.Regions[0].Start == 0 && alternativeMap.MatchedRead.Regions[len(alternativeMap.MatchedRead.Regions)-1].End == len(*read.Sequence) && alternativeMap.MatchedRead.Length() != len(*read.Sequence) {
			// For some reason go decides to share arr space of these objects in alternativeReadMatchResults which
			// leads to unexpected behavior when appending to the mismatch slices of the objects
			// this seems to have fixed it. Before the behavior was that when i appened a mm to the last obj, the last-1 obj
			// somehow also appended this mm and this lead to duplicated mms and a mismatch between query len and cigar len....
			originalMismatches := alternativeMap.MismatchesRead
			alternativeMap.MismatchesRead = make([]int, len(originalMismatches))
			copy(alternativeMap.MismatchesRead, originalMismatches)
			fillGaps(alternativeMap, genomeIndex, read)
		}
	}

	// fill gaps in original map
	if readMapping.MatchedRead.Regions[0].Start == 0 && readMapping.MatchedRead.Regions[len(readMapping.MatchedRead.Regions)-1].End == len(*read.Sequence) && readMapping.MatchedRead.Length() != len(*read.Sequence) {
		// before we do anything, check if the missing part of the map is a gap in the mapping (lets say 10 bases are missing
		// in the middle of the read). We want to use the already mapped part of the read before dicarding them in the remap
		r := readMapping.Copy()
		originalMismatches := r.MismatchesRead
		r.MismatchesRead = make([]int, len(originalMismatches))
		copy(r.MismatchesRead, originalMismatches)
		fillGaps(r, genomeIndex, read)
		alternativeReadMatchResults = append(alternativeReadMatchResults, r)
	}

	if !readMapping.IncompleteMap && readMapping.MatchedRead.Length() == len(*read.Sequence) {
		// after all remap work, append original readMapping if it is a complete remap
		alternativeReadMatchResults = append(alternativeReadMatchResults, readMapping)
	}

	if config.IsOriginRNA && annotation != nil {
		for _, alternative := range alternativeReadMatchResults {
			symIntronCorrected := correctSymmetricIntronErrorsEnhanced(alternative, annotation.Introns[readMapping.SequenceIndex])
			if symIntronCorrected != nil {
				alternativeReadMatchResults = append(alternativeReadMatchResults, symIntronCorrected...)
			}
		}
	}
	return alternativeReadMatchResults
}

// rightRemapAlignmentBlockFromPos remaps a given readSequenceToRemap from a given remapStart position.
func rightRemapAlignmentBlockFromPos(remapStart int, anchorRegion regionvector.Region, targetSeqIntronSet *regionvector.RegionSet, readSequenceToRemap []byte, refSeq *[]byte, startInRead int) ([]*regionvector.Region, int, []int) {
	remap := make([]*regionvector.Region, 0)

	// get next intron
	nextIntronFromAnchor := targetSeqIntronSet.GetNextIntron(anchorRegion.End) // changed .Start to End
	if nextIntronFromAnchor != nil && nextIntronFromAnchor.Start == anchorRegion.End && nextIntronFromAnchor.Rank != len(targetSeqIntronSet.Regions)-1 {
		nextIntronFromAnchor = targetSeqIntronSet.Regions[nextIntronFromAnchor.Rank+1]
	}

	refPos := remapStart
	score := 0
	mm := make([]int, 0)

	for j := 0; j < len(readSequenceToRemap); j++ {
		if refPos == len(*refSeq) {
			// abort in that case
			score = -1
			break
		}

		// jump to next exon in this case
		//        |anchor                  | here we jump
		// (READ) +++++++--------------****####-------------->
		//  (REF) +++++++****--------------------------####++>
		//                  | boundary
		if nextIntronFromAnchor != nil && refPos == nextIntronFromAnchor.Start {
			if remapStart != refPos {
				remap = append(remap, &regionvector.Region{Start: remapStart, End: refPos})
			}

			refPos = nextIntronFromAnchor.End
			remapStart = nextIntronFromAnchor.End
			// update intron if not last intron
			if nextIntronFromAnchor.Rank != len(targetSeqIntronSet.Regions)-1 {
				nextIntronFromAnchor = targetSeqIntronSet.Regions[nextIntronFromAnchor.Rank+1]
			}
			// continue
		}

		if readSequenceToRemap[j] == (*refSeq)[refPos] {
			score++
		} else {
			mm = append(mm, startInRead+j)
		}
		refPos++
	}

	if remapStart != refPos {
		remap = append(remap, &regionvector.Region{Start: remapStart, End: refPos})
	}

	return remap, score, mm
}

// leftRemapAlignmentBlockFromPos remaps a given readSequenceToRemap from a given remapStart position.
func leftRemapAlignmentBlockFromPos(remapStart int, anchorRegion regionvector.Region, targetSeqIntronSet *regionvector.RegionSet, readSequenceToRemap []byte, refSeq *[]byte) ([]*regionvector.Region, int, []int) {
	remap := make([]*regionvector.Region, 0)

	// get next intron
	nextIntronFromAnchor := targetSeqIntronSet.GetPrevIntron(anchorRegion.Start)

	if nextIntronFromAnchor != nil && nextIntronFromAnchor.End == anchorRegion.Start && nextIntronFromAnchor.Rank != 0 {
		nextIntronFromAnchor = targetSeqIntronSet.Regions[nextIntronFromAnchor.Rank-1]
	}
	refPos := remapStart
	mm := make([]int, 0)
	explainedBases := 0
	score := 0

	for j := len(readSequenceToRemap) - 1; j >= 0; j-- {
		if refPos == -1 {
			// if map goes outside ref bounds set score to -1
			score = -1
			break
		}

		// jump to next exon in this case
		//                      | jump here
		// (READ) ----------####****-------+++++++++++---->
		//  (REF) -####----------------****+++++++++++---->
		//
		if nextIntronFromAnchor != nil && refPos == nextIntronFromAnchor.End {
			if refPos != remapStart+1 {
				remap = append(remap, &regionvector.Region{Start: refPos, End: remapStart + 1})
			}

			if readSequenceToRemap[j] != (*refSeq)[refPos] {
				mm = append(mm, j)
			} else {
				score++
			}

			refPos = nextIntronFromAnchor.Start - 1
			remapStart = nextIntronFromAnchor.Start - 1
			// update intron if intron not first intron
			if nextIntronFromAnchor.Rank != 0 {
				nextIntronFromAnchor = targetSeqIntronSet.Regions[nextIntronFromAnchor.Rank-1]
			}
			continue
		}

		if readSequenceToRemap[j] != (*refSeq)[refPos] {
			mm = append(mm, j)
		} else {
			score++
		}

		refPos--
		explainedBases++
	}

	if refPos != remapStart {
		remap = append(remap, &regionvector.Region{Start: refPos + 1, End: remapStart + 1})
	}

	return remap, score, mm
}

// OLD FUNC WHICH IS NO NONGER REQUIRED
// func refineMapping(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, genomeIndex *index.GenomeIndex) *mapperutils.ReadMatchResult {
// 	if readMatchResult.MatchedGenome.HasGaps() {
// 		// get larges anchor in map for right remap
// 		mainAnchor, mainAnchorRank, mainAnchorIndex := readMatchResult.GetLargestAnchor(targetSeqIntronSet)
//
// 		// Check if there are gaps to the right direction of the anchor
// 		// if already corrected the mainAnchor to the right, we don't need to perform a right remap again
// 		if mainAnchorIndex < len(readMatchResult.MatchedGenome.Regions)-1 {
// 			for mainAnchorIndex != len(readMatchResult.MatchedGenome.Regions)-1 {
// 				//////////////////////////////
// 				//////   RIGHT REMAP   ///////
// 				//////////////////////////////
// 				// mainAnchor is left anchor extend right as far as possible
// 				// there are gaps to the right
// 				rGap, ok := readMatchResult.MatchedGenome.GetGap(mainAnchorIndex)
// 				if !ok {
// 					break
// 				}
// 				weakAnchorIndex := mainAnchorIndex + 1
// 				weakAnchor := readMatchResult.MatchedGenome.Regions[weakAnchorIndex]
//
// 				overlappingIntrons := targetSeqIntronSet.GetIntersectingIntrons(rGap)
//
// 				if len(overlappingIntrons) == 0 {
// 					// Case I (Deletion)
// 					// (READ) +++++++-------------+++++++
// 					// (REF) ++++++++++++++++++++++++++++++++
// 					// logrus.WithFields(logrus.Fields{
// 					// 	"Implied gap in read": rGap,
// 					// 	"Length of Deletion":  rGap.Length(),
// 					// 	"Read Blocks":         readMatchResult.MatchedGenome.Regions,
// 					// }).Debug("Found deletion in read.")
// 					// Just log / use for report later
// 					// check if left or right diag can be extended
// 					// startInRead := regionvector.GenomicCoordToReadCoord(readMatchResult.MatchedRead.GetFirstRegion().Start, weakAnchor.Start, readMatchResult.MatchedGenome.Regions)
// 					startInRead := readMatchResult.MatchedRead.Regions[mainAnchorIndex+1].Start
// 					remapOptions := rightRemap(readMatchResult, targetSeqIntronSet, read, weakAnchor, mainAnchor, genomeIndex, mainAnchorRank, startInRead)
// 					if remapOptions != nil {
// 						bestOption := chooseBestRemap(remapOptions)
//
// 						if uint8(float64(len(remapOptions[bestOption].MismatchesRead))*100/float64(len(*read.Sequence))) <= config.MaxMismatchPercentage() {
// 							// remove weakAnchor and everything right to it
// 							lastRegionGenome, _ := readMatchResult.MatchedGenome.GetLastRegion()
// 							readMatchResult.MatchedGenome.RemoveRegion(weakAnchor.Start, lastRegionGenome.End)
// 							// update map
// 							for _, region := range remapOptions[bestOption].RemapVector {
// 								readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(region.Start, region.End)
// 							}
// 							// update mm
// 							readMatchResult.MismatchesRead = remapOptions[bestOption].MismatchesRead
// 						}
// 						readMatchResult.NormalizeRegions()
// 					}
//
// 					mainAnchor = weakAnchor
// 					mainAnchorIndex++
// 					weakAnchorIndex++
// 					continue
// 				}
//
// 				// (READ)      +++-----------------+++++++
// 				// (REF)  ++++++++--------+++------+++++++
// 				//                lIntron    rIntron
// 				// if only one intron is overlapping rIntron == lIntron
// 				lIntron := overlappingIntrons[0]                         // left most intron -> end of alignment block
// 				rIntron := overlappingIntrons[len(overlappingIntrons)-1] // right most intron -> start of next alignment block
//
// 				if lIntron.Start == rGap.Start && rIntron.End == rGap.End {
// 					// if gap already complies with intron boundaries,don't do anything
// 					// and skip to next anchor
// 					mainAnchor = weakAnchor
// 					mainAnchorIndex++
// 					weakAnchorIndex++
// 					continue
// 				}
//
// 				// logrus.WithFields(logrus.Fields{
// 				// 	"Implied gap in read": rGap,
// 				// 	"Read Blocks":         readMatchResult.MatchedGenome.Regions,
// 				// 	"Overlapping Introns": overlappingIntrons,
// 				// }).Debug("Found inconsistent read junction not following inferred intron boundaries")
//
// 				correctedL := lIntron.Start - rGap.Start
// 				correctedR := rIntron.End - rGap.End
//
// 				// if correctedL != correctedR -> remapDiagonal
// 				// this would equal CASE II right remap
// 				// (READ) +++++--------------------++++++
// 				// (REF)  ++++++++-------+++++++++++++++++++++++
// 				// or
// 				// (READ) +++++---------------++++
// 				// (REF)  ++++++++------------------------++++++
// 				if correctedR != correctedL {
//
// 					// startInRead := regionvector.GenomicCoordToReadCoord(readMatchResult.MatchedRead.GetFirstRegion().Start, weakAnchor.Start, readMatchResult.MatchedGenome.Regions)
// 					startInRead := readMatchResult.MatchedRead.Regions[mainAnchorIndex+1].Start
//
// 					// we check for padding if part of the main anchor is reaching into intron
// 					if mainAnchor.End > lIntron.Start && mainAnchor.Start < lIntron.Start {
// 						// check for potential overhangs into neighbor intron of anchor seed
// 						// before:               |remap|
// 						// READ: ++++++++++*-----#######--------->
// 						//  REF: ++++++++++---------------*######>
// 						//                 |lIntron
// 						padding := mainAnchor.End - lIntron.Start
// 						// after:          |   remap   |
// 						// READ: ++++++++++*-----#######--------->
// 						//  REF: ++++++++++---------------*######>
// 						mainAnchor.End -= padding
// 						weakAnchor.Start -= padding
// 						readMatchResult.MatchedGenome.Regions[mainAnchorIndex].End -= padding
// 						readMatchResult.MatchedGenome.Regions[weakAnchorIndex].Start -= padding
// 						startInRead -= padding
// 					}
// 					remapOptions := rightRemap(readMatchResult, targetSeqIntronSet, read, weakAnchor, mainAnchor, genomeIndex, mainAnchorRank, startInRead)
//
// 					if remapOptions != nil {
// 						// add best matching region back to result
// 						// we want the remap with the highest score and if there are remaps with same score, we start with
// 						bestOption := chooseBestRemap(remapOptions)
// 						if uint8(float64(len(remapOptions[bestOption].MismatchesRead))*100/float64(len(*read.Sequence))) <= config.MaxMismatchPercentage() {
//
// 							// remove weakAnchor and everything right to it
// 							lastRegionGenome, _ := readMatchResult.MatchedGenome.GetLastRegion()
// 							readMatchResult.MatchedGenome.RemoveRegion(weakAnchor.Start, lastRegionGenome.End)
//
// 							for _, region := range remapOptions[bestOption].RemapVector {
// 								readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(region.Start, region.End)
// 							}
//
// 							// update mm in readMatchResult
// 							readMatchResult.MismatchesRead = remapOptions[bestOption].MismatchesRead
// 							readMatchResult.NormalizeRegions()
// 						}
// 					}
// 				} else {
// 					// Case III
// 					// (READ) ++++++++++-------+++++++++++
// 					// (READ) ++++++-------+++++++++++++++ (or)
// 					// (REF)  ++++++++-------++++++++++++
// 					// to
// 					// (READ) ++++++++-------+++++++++++++
// 					// (READ) ++++++++-------+++++++++++++
// 					// (REF)  ++++++++-------++++++++++++
// 					if mainAnchor.Length() >= abs(correctedL) && weakAnchor.Length() >= abs(correctedR) {
// 						mainAnchor.End = mainAnchor.End + correctedL
// 						weakAnchor.Start = weakAnchor.Start + correctedR
// 						readMatchResult.MatchedGenome.Regions[mainAnchorIndex].End += correctedL
// 						readMatchResult.MatchedGenome.Regions[weakAnchorIndex].Start += correctedR
//
// 					}
// 				}
// 				// go to next anchor
// 				mainAnchor = weakAnchor
// 				weakAnchorIndex++
// 				mainAnchorIndex++
// 			}
// 		}
//
// 		// get larges anchor in map for left remap
// 		mainAnchor, mainAnchorRank, mainAnchorIndex = readMatchResult.GetLargestAnchor(targetSeqIntronSet)
//
// 		// if already corrected the mainAnchor to the left, we don't need to perform a left remap again
// 		if mainAnchorIndex > 0 {
// 			for mainAnchorIndex > 0 {
// 				//////////////////////////////
// 				//////   LEFT REMAP    ///////
// 				//////////////////////////////
// 				// mainAnchor is right anchor -> extend as far left as possible
// 				// there are gaps to the left
// 				lGap, ok := readMatchResult.MatchedGenome.GetGap(mainAnchorIndex - 1)
// 				if !ok {
// 					break
// 				}
// 				weakAnchorIndex := mainAnchorIndex - 1
// 				weakAnchor := readMatchResult.MatchedGenome.Regions[weakAnchorIndex]
//
// 				overlappingIntrons := targetSeqIntronSet.GetIntersectingIntrons(lGap)
//
// 				if len(overlappingIntrons) == 0 {
// 					// Case I (Deletion)
// 					//                            |mainAnchor
// 					// (READ) ++++----------------+++++++
// 					// (REF) ++++++++++++++++++++++++++++++++
// 					// logrus.WithFields(logrus.Fields{
// 					// 	"Implied gap in read": lGap,
// 					// 	"Length of Deletion":  lGap.Length(),
// 					// 	"Read Blocks":         readMatchResult.MatchedGenome.Regions,
// 					// }).Debug("Found deletion in read.")
// 					// Just log / use for report later
// 					// check if left or right diag can be extended
// 					regionReadEnd, err := regionvector.GenomicCoordToReadCoord(0, weakAnchor.End, readMatchResult.MatchedGenome.Regions, len(*read.Sequence))
// 					if err != nil {
// 						println(read.Header)
// 						fmt.Println("LEFT REMAP")
//
// 						fmt.Println(weakAnchor)
// 						fmt.Println(readMatchResult.MatchedGenome.Regions)
// 						fmt.Println(readMatchResult.MatchedRead.Regions)
//
// 						logrus.Errorf("Error while converting genomic coord to read coord in read %s", read.Header)
// 						logrus.Fatal(err)
// 					}
// 					remapOptions := leftRemap(readMatchResult, targetSeqIntronSet, read, weakAnchor, mainAnchor, genomeIndex, mainAnchorRank, regionReadEnd)
// 					if remapOptions != nil {
//
// 						bestOption := chooseBestRemap(remapOptions)
//
// 						if uint8(float64(len(remapOptions[bestOption].MismatchesRead))*100/float64(len(*read.Sequence))) <= config.MaxMismatchPercentage() {
// 							// remove weakAnchor and everything left from it
//
// 							firstRegionGenome, _ := readMatchResult.MatchedGenome.GetFirstRegion()
// 							readMatchResult.MatchedGenome.RemoveRegion(firstRegionGenome.Start, weakAnchor.End)
// 							for _, region := range remapOptions[bestOption].RemapVector {
// 								readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(region.Start, region.End)
// 							}
//
// 							readMatchResult.MismatchesRead = remapOptions[bestOption].MismatchesRead
// 							readMatchResult.NormalizeRegions()
// 						}
// 					}
//
// 					mainAnchor = weakAnchor
// 					mainAnchorIndex--
// 					weakAnchorIndex--
// 					continue
// 				}
//
// 				//             |weakAnchor             |mainAnchor
// 				// (READ)      +++---------------------+++++++
// 				// (REF)  ++++++++---------+++---------+++++++
// 				//                |lIntron|   |rIntron|
// 				// if only one intron is overlapping rIntron == lIntron
// 				lIntron := overlappingIntrons[0]
// 				rIntron := overlappingIntrons[len(overlappingIntrons)-1]
//
// 				if lIntron.Start == lGap.Start && rIntron.End == lGap.End {
// 					// if gap already complies with intron boundaries, don't do anything and skip to next anchor
// 					mainAnchor = weakAnchor
// 					mainAnchorIndex--
// 					weakAnchorIndex--
// 					continue
// 				}
//
// 				// logrus.WithFields(logrus.Fields{
// 				// 	"Implied gap in read": lGap,
// 				// 	"Read Blocks":         readMatchResult.MatchedGenome.Regions,
// 				// 	"Overlapping Introns": overlappingIntrons,
// 				// }).Debug("Found inconsistent read junction not following inferred intron boundaries")
//
// 				correctedL := lIntron.Start - lGap.Start
// 				correctedR := rIntron.End - lGap.End
//
// 				// if correctedL != correctedR -> remapDiagonal
// 				// this would equal CASE II left remap
// 				//                              |mainAnchor
// 				// (READ) +++++-----------------++++++++++++
// 				// (REF)  ++++++++-------+++++++++++++++++++++++
// 				// or                                    |mainAnchor
// 				// (READ)            +++++---------------++++++++++
// 				// (REF)  ++++++++-----------------------++++++++++
// 				if correctedR != correctedL {
// 					// we check for padding if part of the main anchor is reaching into intron
// 					if mainAnchor.Start < rIntron.End && mainAnchor.End > rIntron.End {
// 						// before:      |remap|   mainAnchor
// 						//(READ) -------####------*+++++---------->
// 						//(REF)  ####*-------------+++++---------->
// 						//            |  rIntron  |
// 						padding := rIntron.End - mainAnchor.Start
// 						// after:       |  remap  |
// 						//(READ) -------####------*+++++---------->
// 						//(REF)  ####*-------------+++++---------->
// 						mainAnchor.Start += padding
// 						weakAnchor.End += padding
// 						readMatchResult.MatchedGenome.Regions[mainAnchorIndex].Start += padding
// 						readMatchResult.MatchedGenome.Regions[weakAnchorIndex].End += padding
//
// 					}
//
// 					regionReadEnd, err := regionvector.GenomicCoordToReadCoord(0, weakAnchor.End, readMatchResult.MatchedGenome.Regions, len(*read.Sequence))
// 					if err != nil {
// 						println(read.Header)
// 						fmt.Println("LEFT REMAP")
// 						fmt.Println(weakAnchor)
// 						fmt.Println(readMatchResult.MatchedGenome.Regions)
// 						fmt.Println(readMatchResult.MatchedRead.Regions)
// 						logrus.Errorf("Error while converting genomic coord to read coord in read %s", read.Header)
// 						logrus.Fatal(err)
// 					}
// 					remapOptions := leftRemap(readMatchResult, targetSeqIntronSet, read, weakAnchor, mainAnchor, genomeIndex, mainAnchorRank, regionReadEnd)
//
// 					if remapOptions != nil {
// 						bestOption := chooseBestRemap(remapOptions)
//
// 						if uint8(float64(len(remapOptions[bestOption].MismatchesRead))*100/float64(len(*read.Sequence))) <= config.MaxMismatchPercentage() {
// 							// remove weakAnchor and everything left from it
// 							firstRegionGenome, _ := readMatchResult.MatchedGenome.GetFirstRegion()
// 							readMatchResult.MatchedGenome.RemoveRegion(firstRegionGenome.Start, weakAnchor.End)
//
// 							// use remap
// 							for _, region := range remapOptions[bestOption].RemapVector {
// 								readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(region.Start, region.End)
// 							}
// 							// update mm in readMatchResult
// 							readMatchResult.MismatchesRead = remapOptions[bestOption].MismatchesRead
// 							readMatchResult.NormalizeRegions()
// 						}
// 					}
// 				} else {
// 					// Case III
// 					// (READ) ++++++++++-------+++++++++++
// 					// (READ) ++++++-------+++++++++++++++ (or)
// 					// (REF)  ++++++++-------++++++++++++
// 					// to
// 					// (READ) ++++++++-------+++++++++++++
// 					// (READ) ++++++++-------+++++++++++++
// 					// (REF)  ++++++++-------++++++++++++
// 					if mainAnchor.Length() > abs(correctedL) && weakAnchor.Length() > abs(correctedR) {
// 						mainAnchor.Start = mainAnchor.Start + correctedL
// 						weakAnchor.End = weakAnchor.End + correctedR
// 						readMatchResult.MatchedGenome.Regions[mainAnchorIndex].Start += correctedL
// 						readMatchResult.MatchedGenome.Regions[weakAnchorIndex].End += correctedR
// 					}
// 				}
// 				// go to next anchor
// 				mainAnchor = weakAnchor
// 				mainAnchorIndex--
// 				weakAnchorIndex--
// 			}
// 		}
// 	}
//
// 	////////////////////////////////////////////////////////
// 	///// CHECK IF WE NEED TO CORRECT START OR END REGION //
// 	////////////////////////////////////////////////////////
// 	correctOverhangsOLD(readMatchResult, targetSeqIntronSet, read, genomeIndex) // check for first region
//
// 	// in any case we want to reprocess the remapped read again in order to annotate splicesites and new mm (caused by padding)
// 	mm := make([]int, 0)
// 	// readMatchResult.SyncRegions()
// 	readMatchResult.NormalizeRegions()
// 	for i := 0; i < len(readMatchResult.MatchedRead.Regions); i++ {
// 		readRegion := readMatchResult.MatchedRead.Regions[i]
// 		genomeRegion := readMatchResult.MatchedGenome.Regions[i]
// 		for j := 0; j < readRegion.Length(); j++ {
// 			posRead := readRegion.Start + j
// 			posGenome := genomeRegion.Start + j
// 			if (*read.Sequence)[posRead] != (*genomeIndex.Sequences[readMatchResult.SequenceIndex])[posGenome] {
// 				mm = append(mm, readRegion.Start+j)
// 			}
// 		}
// 	}
//
// 	if readMatchResult.MatchedGenome.Length() != len(*read.Sequence) {
// 		return nil
// 	}
//
// 	// check if read could be mapped/remapped
// 	if uint8(float64(len(readMatchResult.MismatchesRead))*100/float64(len(*read.Sequence))) > config.MaxMismatchPercentage() {
// 		readMatchResult.IncompleteMap = true
// 		return nil
// 	} else {
// 		readMatchResult.IncompleteMap = false
// 		return readMatchResult
// 	}
// }

func correctSymmetricIntronErrors(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, genomeIndex *index.GenomeIndex) *mapperutils.ReadMatchResult {
	correctedReadMatchResult := readMatchResult.Copy()
	hasCorrection := false

	// iterate over all ali blocks and check if sym intron error
	for i := 0; i <= len(readMatchResult.MatchedGenome.Regions)-2; i++ {
		leftMappedRegion := readMatchResult.MatchedGenome.Regions[i]
		rightMappedRegion := readMatchResult.MatchedGenome.Regions[i+1]
		gap, ok := readMatchResult.MatchedGenome.GetGap(i)
		if !ok {
			break
		}

		overlappingIntrons := targetSeqIntronSet.GetIntersectingIntrons(gap)

		if len(overlappingIntrons) == 0 {
			// go to next gap
			continue
		}

		if len(overlappingIntrons) == 1 {
			intron := overlappingIntrons[0]

			correctedL := intron.Start - gap.Start
			correctedR := intron.End - gap.End

			if correctedL == correctedR && correctedL != 0 {
				if leftMappedRegion.Length() >= utils.Abs(correctedL) && rightMappedRegion.Length() >= utils.Abs(correctedR) {
					correctedReadMatchResult.MatchedGenome.Regions[i].End += correctedL
					correctedReadMatchResult.MatchedGenome.Regions[i+1].Start += correctedR
					hasCorrection = true
					correctedReadMatchResult.IsSymInErr = true
					correctedReadMatchResult.SymInErrLen = append(correctedReadMatchResult.SymInErrLen, utils.Abs(correctedL))
				}
			}
		} else {
			preferedIntron := -1
			score := -1
			count := 0
			for i, overlapIntron := range overlappingIntrons {
				if score == -1 {
					preferedIntron = i
					score = overlapIntron.SpliceSiteScore
					continue
				}
				if overlapIntron.SpliceSiteScore > score {
					preferedIntron = i
					score = overlapIntron.SpliceSiteScore
					continue
				}
				if overlapIntron.SpliceSiteScore == score {
					count++
				}
			}

			// if score == 0 {
			// 	continue // dont bother correction this junction since spliceSite score is 0
			// }

			if preferedIntron == -1 {
				continue
			}

			if count > 1 {
				logrus.Warnf("Potentially missed padding in sym intron err with several introns having score greater than 0 %s", read.Header)
			}

			correctedL := overlappingIntrons[preferedIntron].Start - gap.Start
			correctedR := overlappingIntrons[preferedIntron].End - gap.End
			if correctedL == correctedR && correctedL != 0 {
				if leftMappedRegion.Length() >= utils.Abs(correctedL) && rightMappedRegion.Length() >= utils.Abs(correctedR) {
					correctedReadMatchResult.MatchedGenome.Regions[i].End += correctedL
					correctedReadMatchResult.MatchedGenome.Regions[i+1].Start += correctedR
					hasCorrection = true
					correctedReadMatchResult.IsSymInErr = true
					correctedReadMatchResult.SymInErrLen = append(correctedReadMatchResult.SymInErrLen, utils.Abs(correctedL))
				}
			}
		}
	}

	if hasCorrection {
		return correctedReadMatchResult
	} else {
		return nil
	}
}

func cartesianProduct(corrections map[int][]int) [][]int {
	if len(corrections) == 0 {
		return nil
	}

	keys := make([]int, 0, len(corrections))
	for k := range corrections {
		keys = append(keys, k)
	}
	sort.Ints(keys)

	var results [][]int

	var backtrack func(idx int, current []int)
	backtrack = func(idx int, current []int) {
		if idx == len(keys) {
			combo := make([]int, len(current))
			copy(combo, current)
			results = append(results, combo)
			return
		}

		k := keys[idx]
		for _, v := range corrections[k] {
			backtrack(idx+1, append(current, v))
		}
	}

	backtrack(0, []int{})
	return results
}

func correctSymmetricIntronErrorsEnhanced(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet) []*mapperutils.ReadMatchResult {
	// collect all sym paddings into this map
	corrections := make(map[int][]int)

	// iterate over all ali blocks and check if sym intron error
	for i := 0; i <= len(readMatchResult.MatchedGenome.Regions)-2; i++ {
		leftMappedRegion := readMatchResult.MatchedGenome.Regions[i]
		rightMappedRegion := readMatchResult.MatchedGenome.Regions[i+1]
		gap, ok := readMatchResult.MatchedGenome.GetGap(i)
		if !ok {
			break
		}

		overlappingIntrons := targetSeqIntronSet.GetIntersectingIntrons(gap)
		for _, overlapIntron := range overlappingIntrons {
			correctedL := overlapIntron.Start - gap.Start
			correctedR := overlapIntron.End - gap.End
			if correctedL == correctedR && correctedL != 0 {
				if leftMappedRegion.Length() >= utils.Abs(correctedL) && rightMappedRegion.Length() >= utils.Abs(correctedR) {
					corrections[i] = append(corrections[i], correctedL)
				}
			}
		}

	}

	combinations := cartesianProduct(corrections) // in combinations, each sub slice is a comb of paddings for the readMatchResult. Each index is a gap
	correctedRes := make([]*mapperutils.ReadMatchResult, 0)
	for _, comb := range combinations {
		corrected := readMatchResult.Copy()
		for idxGap, padding := range comb {
			corrected.MatchedGenome.Regions[idxGap].End += padding
			corrected.MatchedGenome.Regions[idxGap+1].Start += padding
			corrected.IsSymInErr = true
			corrected.SymInErrLen = append(corrected.SymInErrLen, utils.Abs(padding))
		}
		correctedRes = append(correctedRes, corrected)
	}

	if len(correctedRes) != 0 {
		return correctedRes
	} else {
		return nil
	}
}

func getPossiblePaddingCombinations(leftMappedRegion, rightMappedRegion regionvector.Region, lPaddings, rPaddings map[int]struct{}) []int {
	paddingsToCheck := make([]int, 0)
	for lPad := range lPaddings {
		for rPad := range rPaddings {
			if lPad == rPad {
				if leftMappedRegion.Length() >= utils.Abs(lPad) && rightMappedRegion.Length() >= utils.Abs(rPad) {
					paddingsToCheck = append(paddingsToCheck, lPad) // since lPad == rPad, we can just append one of them to our list
				}
			}
		}
	}
	return paddingsToCheck
}

func correctOverhangs(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, genomeIndex *index.GenomeIndex, intronTree *datastructure.INTree) []*mapperutils.ReadMatchResult {
	remaps := make([]*mapperutils.ReadMatchResult, 0)

	if len(readMatchResult.MatchedGenome.Regions) > 1 {
		// we need to check left most anchor and right most anchor
		firstRegion := readMatchResult.MatchedGenome.Regions[0]
		firstRegionRead := readMatchResult.MatchedRead.Regions[0]

		lastRegion := readMatchResult.MatchedGenome.Regions[len(readMatchResult.MatchedGenome.Regions)-1]
		lastRegionRead := readMatchResult.MatchedRead.Regions[len(readMatchResult.MatchedRead.Regions)-1]

		// here we need to do it manually
		lPaddings, _ := getPadding(firstRegion, targetSeqIntronSet, intronTree)
		_, rPaddings := getPadding(lastRegion, targetSeqIntronSet, intronTree)

		adjustedAnchorsRight := make([]regionvector.Region, 0)
		adjustedReadAnchorsRight := make([]regionvector.Region, 0)

		adjustedAnchorsLeft := make([]regionvector.Region, 0)
		adjustedReadAnchorsLeft := make([]regionvector.Region, 0)

		if len(rPaddings) > 0 {
			for r := range rPaddings {
				adjustedAnchor := lastRegion.Copy()
				adjustedReadAnchor := lastRegionRead.Copy()
				if lastRegionRead.End-r > lastRegionRead.Start {
					adjustedAnchor.End -= r
					adjustedReadAnchor.End -= r
				}

				adjustedAnchorsRight = append(adjustedAnchorsRight, adjustedAnchor)
				adjustedReadAnchorsRight = append(adjustedReadAnchorsRight, adjustedReadAnchor)
			}
		}

		if len(lPaddings) > 0 {
			for l := range lPaddings {
				adjustedAnchor := firstRegion.Copy()
				adjustedReadAnchor := firstRegionRead.Copy()
				if firstRegionRead.Start+l < firstRegionRead.End {
					adjustedAnchor.Start += l
					adjustedReadAnchor.Start += l
				}

				adjustedAnchorsLeft = append(adjustedAnchorsLeft, adjustedAnchor)
				adjustedReadAnchorsLeft = append(adjustedReadAnchorsLeft, adjustedReadAnchor)
			}
		}

		leftSections := make([]*RemapSection, 0)
		for i := 0; i < len(adjustedAnchorsLeft); i++ {
			// correct first region
			missingBases := adjustedReadAnchorsLeft[i].Start
			startPos := adjustedAnchorsLeft[i].Start
			leftPaths := targetSeqIntronSet.TranscriptomeGraph.FindPathsLeft(startPos, missingBases)
			if len(leftPaths) != 0 {
				remapSectionsLeft := scoreLeftOptions(leftPaths, adjustedReadAnchorsLeft[i].Start, read.Sequence, genomeIndex.Sequences[readMatchResult.SequenceIndex], adjustedAnchorsLeft[i])
				leftSections = append(leftSections, remapSectionsLeft...)
			}
		}
		// correct last region
		rightSections := make([]*RemapSection, 0)
		for i := 0; i < len(adjustedReadAnchorsRight); i++ {
			missingBases := len(*read.Sequence) - adjustedReadAnchorsRight[i].End
			startPos := adjustedAnchorsRight[i].End
			rightPaths := targetSeqIntronSet.TranscriptomeGraph.FindPathsRight(startPos, missingBases)
			if len(rightPaths) != 0 {
				remapSectionsRight := scoreRightOptions(rightPaths, adjustedReadAnchorsRight[i].End, read.Sequence, genomeIndex.Sequences[readMatchResult.SequenceIndex], adjustedAnchorsRight[i])
				rightSections = append(rightSections, remapSectionsRight...)
			}
		}

		rightRemaps := chooseBestRemapSection(rightSections)
		leftRemaps := chooseBestRemapSection(leftSections)

		// pair remaining sections in leftRemaps and rightRemaps
		finalRemaps := make([]*mapperutils.ReadMatchResult, 0)

		if rightRemaps != nil && leftRemaps != nil {
			for _, lSection := range leftRemaps {
				lExtStopRead := lSection.MatchedRead[0].End
				for _, rSection := range rightRemaps {
					rExtStartRead := rSection.MatchedRead[0].Start

					// create template correction once (readMatchResult with regions removed which got remapped)
					correctedTemplate := readMatchResult.Copy()
					// remove left regions
					correctedTemplate.MatchedGenome.RemoveRegion(0, lSection.MainAnchor.Start)
					// remove right regions
					correctedTemplate.MatchedGenome.RemoveRegion(rSection.MainAnchor.End, len(*genomeIndex.Sequences[readMatchResult.SequenceIndex]))

					templateMappedRead := regionvector.Region{Start: lExtStopRead, End: rExtStartRead}
					templateMM := extractMMofAnchor(templateMappedRead, correctedTemplate.MismatchesRead)
					corrected := correctedTemplate.Copy()

					for _, region := range lSection.MatchedGenome {
						corrected.MatchedGenome.AddRegionNonOverlappingPanic(region.Start, region.End)
					}

					for _, region := range rSection.MatchedGenome {
						corrected.MatchedGenome.AddRegionNonOverlappingPanic(region.Start, region.End)
					}
					corrected.MismatchesRead = templateMM
					corrected.MismatchesRead = append(corrected.MismatchesRead, lSection.Mm...)
					corrected.MismatchesRead = append(corrected.MismatchesRead, rSection.Mm...)

					corrected.TotalLeftOptions = lSection.TotalLeftPaths
					corrected.TotalRightOptions = rSection.TotalRightPaths
					corrected.ValidLeftOptions = len(leftRemaps)
					corrected.ValidRightOptions = len(rightRemaps)
					finalRemaps = append(finalRemaps, corrected)
				}
			}
		} else if rightRemaps != nil {
			for _, rSection := range rightRemaps {
				rExtStartRead := rSection.MatchedRead[len(rSection.MatchedRead)-1].Start
				// create template correction once (readMatchResult with regions removed which got remapped)
				correctedTemplate := readMatchResult.Copy()
				// remove right regions
				correctedTemplate.MatchedGenome.RemoveRegion(rSection.MainAnchor.End, len(*genomeIndex.Sequences[readMatchResult.SequenceIndex]))
				templateMappedRead := regionvector.Region{Start: 0, End: rExtStartRead}
				templateMM := extractMMofAnchor(templateMappedRead, correctedTemplate.MismatchesRead)
				corrected := correctedTemplate.Copy()
				for _, region := range rSection.MatchedGenome {
					corrected.MatchedGenome.AddRegionNonOverlappingPanic(region.Start, region.End)
				}
				corrected.MismatchesRead = templateMM
				corrected.MismatchesRead = append(corrected.MismatchesRead, rSection.Mm...)
				corrected.TotalRightOptions = rSection.TotalRightPaths
				corrected.ValidRightOptions = len(rightRemaps)
				finalRemaps = append(finalRemaps, corrected)
			}
		} else if leftRemaps != nil {
			for _, lSection := range leftRemaps {
				lExtStopRead := lSection.MatchedRead[0].End
				// create template correction once (readMatchResult with regions removed which got remapped)
				correctedTemplate := readMatchResult.Copy()
				// remove left regions
				correctedTemplate.MatchedGenome.RemoveRegion(0, lSection.MainAnchor.Start)
				templateMappedRead := regionvector.Region{Start: lExtStopRead, End: len(*read.Sequence)}
				templateMM := extractMMofAnchor(templateMappedRead, correctedTemplate.MismatchesRead)
				corrected := correctedTemplate.Copy()
				for _, region := range lSection.MatchedGenome {
					corrected.MatchedGenome.AddRegionNonOverlappingPanic(region.Start, region.End)
				}
				corrected.MismatchesRead = templateMM
				corrected.MismatchesRead = append(corrected.MismatchesRead, lSection.Mm...)
				corrected.TotalLeftOptions = lSection.TotalLeftPaths
				corrected.ValidLeftOptions = len(leftRemaps)
				finalRemaps = append(finalRemaps, corrected)
			}
		}

		for _, r := range finalRemaps {
			r.IsOverhangCorrected = true
		}
		return finalRemaps

	} else {
		// we can do same strategy as in fixPointRNARemap since there is only one anchor
		overhangCorrected := fixPointRNARemap(readMatchResult, targetSeqIntronSet, read, genomeIndex, intronTree)
		for _, r := range overhangCorrected {
			r.IsOverhangCorrected = true
		}
		remaps = append(remaps, overhangCorrected...)
	}
	return remaps
}

func chooseBestRemapSection(sections []*RemapSection) []*RemapSection {
	minSections := make([]*RemapSection, 0)
	minMM := 50000
	for _, s := range sections {
		if len(s.Mm) < minMM {
			minMM = len(s.Mm)
		}
	}

	for _, s := range sections {
		if len(s.Mm) == minMM {
			minSections = append(minSections, s)
		}
	}

	return minSections
}

// func correctOverhangsOLD(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, genomeIndex *index.GenomeIndex) {
// 	firstRegion, _ := readMatchResult.MatchedGenome.GetFirstRegion()
// 	overlappingIntrons := targetSeqIntronSet.GetIntersectingIntrons(firstRegion)
// 	if len(overlappingIntrons) != 0 {
// 		// solid map, no remap possible
// 		mapStart := firstRegion.Start
// 		mapEnd := firstRegion.End
// 		// either
// 		// +++++++++ (READ)
// 		// ----+++++(REF)
// 		// or
// 		// +++++++++ (READ)
// 		// +++++---- (REF)
// 		// or  (unlikely i think)
// 		// ++++++++++++ (READ)
// 		// ----+++++---(REF)
// 		// or  (unlikely i think)
// 		// +++++++++++++++++++++++++++++ (READ)
// 		// +++++++++++++++----++++++++++(REF) (deletion) -> we do not want to correct here
// 		lIntron := overlappingIntrons[0]
// 		rIntron := overlappingIntrons[len(overlappingIntrons)-1]
//
// 		// handle deletion
// 		if lIntron.Start > mapStart && rIntron.End < mapEnd {
// 			return
// 		}
//
// 		// check if read start needs left remap
// 		if lIntron.End > mapStart && mapEnd > lIntron.End {
// 			//         |mapStart
// 			//         +++++++++ (READ)
// 			// +++---------+++++(REF)
// 			//            |lIntron.End
//
// 			// adjust anchors
// 			nextAnchor := regionvector.Region{Start: firstRegion.Start, End: lIntron.End}
// 			firstRegion.Start = lIntron.End
// 			readMatchResult.MatchedGenome.Regions[0].Start = lIntron.End
// 			mainAnchorRank := lIntron.Rank + 1
// 			readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(nextAnchor.Start, nextAnchor.End)
// 			regionReadEnd, err := regionvector.GenomicCoordToReadCoord(0, nextAnchor.End, readMatchResult.MatchedGenome.Regions, len(*read.Sequence))
// 			if err != nil {
// 				println(read.Header)
// 				fmt.Println("CORRECT OVERHANGS")
// 				fmt.Println(nextAnchor)
// 				fmt.Println(readMatchResult.MatchedGenome.Regions)
// 				fmt.Println(readMatchResult.MatchedRead.Regions)
// 				logrus.Errorf("Error while converting genomic coord to read coord in read %s", read.Header)
// 				logrus.Fatal(err)
// 			}
// 			remapOptions := leftRemap(readMatchResult, targetSeqIntronSet, read, nextAnchor, firstRegion, genomeIndex, mainAnchorRank, regionReadEnd)
//
// 			if remapOptions != nil {
// 				bestOption := chooseBestRemap(remapOptions)
// 				if uint8(float64(len(remapOptions[bestOption].MismatchesRead))*100/float64(len(*read.Sequence))) <= config.MaxMismatchPercentage() {
//
// 					// remove weakAnchor
// 					firstRegionGenome, _ := readMatchResult.MatchedGenome.GetFirstRegion()
// 					readMatchResult.MatchedGenome.RemoveRegion(firstRegionGenome.Start, nextAnchor.End)
//
// 					// use remap
// 					for _, r := range remapOptions[bestOption].RemapVector {
// 						readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(r.Start, r.End)
// 					}
// 					// update mm in readMatchResult
// 					readMatchResult.MismatchesRead = remapOptions[bestOption].MismatchesRead
// 				}
// 			}
// 		}
// 	}
//
// 	lastRegion, _ := readMatchResult.MatchedGenome.GetLastRegion()
// 	overlappingIntrons = targetSeqIntronSet.GetIntersectingIntrons(lastRegion)
// 	if len(overlappingIntrons) == 0 {
// 		// solid map, no remap possible
// 		return
// 	}
// 	mapStart := lastRegion.Start
// 	mapEnd := lastRegion.End
// 	rIntron := overlappingIntrons[len(overlappingIntrons)-1]
//
// 	// check if read end needs right remap
// 	if rIntron.Start < mapEnd && mapStart < rIntron.Start {
// 		// remap
// 		// +++++++++ (READ)
// 		// +++++---------++++++++ (REF)
// 		// adjust anchor
// 		nextAnchor := regionvector.Region{Start: rIntron.Start, End: mapEnd}
// 		lastRegion.End = rIntron.Start
// 		readMatchResult.MatchedGenome.Regions[len(readMatchResult.MatchedGenome.Regions)-1].End = rIntron.Start // adjust main anchor overhang
// 		readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(nextAnchor.Start, nextAnchor.End)            // add region to remap back as separate block
//
// 		mainAnchorRank := rIntron.Rank - 1
// 		firstRegionRead, _ := readMatchResult.MatchedRead.GetFirstRegion()
//
// 		startInRead, err := regionvector.GenomicCoordToReadCoord(firstRegionRead.Start, nextAnchor.Start, readMatchResult.MatchedGenome.Regions, len(*read.Sequence))
// 		if err != nil {
// 			println(read.Header)
// 			fmt.Println("CORRECT OVERHANGS")
// 			fmt.Println(nextAnchor)
// 			fmt.Println(readMatchResult.MatchedGenome.Regions)
// 			fmt.Println(readMatchResult.MatchedRead.Regions)
// 			logrus.Errorf("Error while converting genomic coord to read coord")
// 			logrus.Fatal(err)
// 		}
// 		remapOptions := rightRemap(readMatchResult, targetSeqIntronSet, read, nextAnchor, lastRegion, genomeIndex, mainAnchorRank, startInRead)
// 		if remapOptions != nil {
// 			bestOption := chooseBestRemap(remapOptions)
// 			if uint8(float64(len(remapOptions[bestOption].MismatchesRead))*100/float64(len(*read.Sequence))) <= config.MaxMismatchPercentage() {
// 				// remove weakAnchor
// 				lastRegionGenome, _ := readMatchResult.MatchedGenome.GetLastRegion()
// 				readMatchResult.MatchedGenome.RemoveRegion(nextAnchor.Start, lastRegionGenome.End)
//
// 				// use remap
// 				for _, r := range remapOptions[bestOption].RemapVector {
// 					readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(r.Start, r.End)
// 				}
// 				// update mm in readMatchResult
// 				readMatchResult.MismatchesRead = remapOptions[bestOption].MismatchesRead
// 			}
// 		}
// 	}
// }

// func rightRemap(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, weakAnchor regionvector.Region, anchorRegion regionvector.Region, genomeIndex *index.GenomeIndex, anchorRank int, startInRead int) map[int]*RemapOption {
// 	refSeq := genomeIndex.Sequences[readMatchResult.SequenceIndex]
// 	if weakAnchor.End > len(*refSeq) {
// 		return nil
// 	}
// 	remapOptions := make(map[int]*RemapOption)
//
// 	regionReadEnd := len(*read.Sequence)
//
// 	if regionReadEnd > len(*read.Sequence) || startInRead < 0 {
// 		logrus.Infof("Aborting right remap to prevent index error in %s, sindex=%d", read.Header, readMatchResult.SequenceIndex)
// 		logrus.Infof("Anchor region end: %d", regionReadEnd)
// 		fmt.Println(readMatchResult.MatchedGenome)
// 		fmt.Println(readMatchResult.MatchedRead)
// 		return nil
// 	}
//
// 	readSequenceToRemap := (*read.Sequence)[startInRead:regionReadEnd]
//
// 	// get next intron
// 	intronRight := targetSeqIntronSet.GetNextIntron(anchorRegion.Start)
// 	if intronRight != nil {
// 		if anchorRegion.End < intronRight.Start {
// 			// first check if end of anchorRegion is at intron boundary
// 			// if now we need to also remap from anchorRegion.End
// 			//               |anchorRegion.End
// 			// (READ) +++++++--------------++---------->
// 			//  (REF) +++++++++++------------------++++>
// 			//                   | intronRight.Start
// 			remapPos := anchorRegion.End
// 			regions, score, mm := rightRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq, startInRead)
// 			mmCombined := alignMismatches(readMatchResult.MismatchesRead, mm, true, startInRead)
// 			remapOptions[remapPos] = &RemapOption{
// 				RemapStart:       remapPos,
// 				RemapVector:      regions,
// 				MismatchesRead:   mmCombined,
// 				Score:            score,
// 				DistToMainAnchor: 0, // since remap starts within same exon as anchor
// 			}
// 		}
//
// 		// remap for all remaining exons (in right direction)
// 		for i := intronRight.Rank; i < len(targetSeqIntronSet.Regions); i++ {
// 			remapPos := targetSeqIntronSet.Regions[i].End
// 			regions, score, mm := rightRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq, startInRead)
// 			mmCombined := alignMismatches(readMatchResult.MismatchesRead, mm, true, startInRead)
// 			remapOptions[remapPos] = &RemapOption{
// 				RemapStart:       remapPos,
// 				RemapVector:      regions,
// 				MismatchesRead:   mmCombined,
// 				Score:            score,
// 				DistToMainAnchor: utils.Abs(i - anchorRank), // since remap starts after intron
// 			}
// 		}
// 	} else {
// 		// if there are no introns left in right direction, start remap from anchorregion end
// 		remapPos := anchorRegion.End
// 		regions, score, mm := rightRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq, startInRead)
// 		mmCombined := alignMismatches(readMatchResult.MismatchesRead, mm, true, startInRead)
// 		remapOptions[remapPos] = &RemapOption{
// 			RemapStart:       remapPos,
// 			RemapVector:      regions,
// 			MismatchesRead:   mmCombined,
// 			Score:            score,
// 			DistToMainAnchor: utils.Abs(targetSeqIntronSet.Regions[len(targetSeqIntronSet.Regions)-1].Rank - anchorRank), // remap is after last intron
// 		}
// 	}
//
// 	return remapOptions
// }

// func leftRemap(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, weakAnchor regionvector.Region, anchorRegion regionvector.Region, genomeIndex *index.GenomeIndex, anchorRank, regionReadEnd int) map[int]*RemapOption {
// 	if weakAnchor.Start < 0 {
// 		// if this is the case, theres not enough gene left to remap the read
// 		// GENOME [2,122]
// 		// READ   [30,150] -> 30 bases missing but cant do left remap since only 2 bases in gene are left
// 		// then weakAnchor == [-28, 2]
// 		// retrun and dont do anything
// 		return nil
// 	}
// 	remapOptions := make(map[int]*RemapOption)
//
// 	regionReadStart := 0 // remap entire read from 0 to end of weakAnchor
// 	if regionReadEnd > len(*read.Sequence) {
// 		logrus.Infof("Aborting left remap to prevent index error in %s, sindex=%d", read.Header, readMatchResult.SequenceIndex)
// 		logrus.Infof("Anchor region end: %d", regionReadEnd)
// 		fmt.Println(readMatchResult.MatchedGenome)
// 		fmt.Println(readMatchResult.MatchedRead)
// 		return nil
// 	}
//
// 	readSequenceToRemap := (*read.Sequence)[regionReadStart:regionReadEnd]
// 	refSeq := genomeIndex.Sequences[readMatchResult.SequenceIndex]
//
// 	// get next intron
// 	intronLeft := targetSeqIntronSet.GetPrevIntron(anchorRegion.Start)
// 	if intronLeft != nil {
// 		if anchorRegion.Start > intronLeft.End {
// 			// check if start of anchorRegion is at intron boundary
// 			//                                 |anchor
// 			// (READ) -----++------------------+++++++->
// 			//  (REF) +++--------------+++++++++++++++->
// 			//                        |intronLeft.End
// 			remapPos := anchorRegion.Start - 1
// 			regions, score, mm := leftRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
// 			mmCombined := alignMismatches(readMatchResult.MismatchesRead, mm, false, regionReadEnd)
//
// 			remapOptions[remapPos] = &RemapOption{
// 				RemapStart:       remapPos,
// 				RemapVector:      regions,
// 				MismatchesRead:   mmCombined,
// 				Score:            score,
// 				DistToMainAnchor: 0, // since remap starts within same exon as anchor
// 			}
// 		}
//
// 		// remap for remaining exons (in left direction)
// 		for i := intronLeft.Rank; i >= 0; i-- {
// 			remapPos := targetSeqIntronSet.Regions[i].Start - 1
// 			regions, score, mm := leftRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
// 			mmCombined := alignMismatches(readMatchResult.MismatchesRead, mm, false, regionReadEnd)
// 			remapOptions[remapPos] = &RemapOption{
// 				RemapStart:       remapPos,
// 				RemapVector:      regions,
// 				MismatchesRead:   mmCombined,
// 				Score:            score,
// 				DistToMainAnchor: utils.Abs(anchorRank-i) + 1,
// 			}
// 		}
// 	} else {
// 		// if there are no introns left in left direction, start remap from anchorregion start - 1
// 		remapPos := anchorRegion.Start - 1
// 		regions, score, mm := leftRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
// 		mmCombined := alignMismatches(readMatchResult.MismatchesRead, mm, false, regionReadEnd)
// 		remapOptions[remapPos] = &RemapOption{
// 			RemapStart:       remapPos,
// 			RemapVector:      regions,
// 			MismatchesRead:   mmCombined,
// 			Score:            score,
// 			DistToMainAnchor: anchorRank, // remap is before first intron
// 		}
// 	}
// 	return remapOptions
// }

type RemapOption struct {
	RemapStart       int
	RemapVector      []*regionvector.Region
	MismatchesRead   []int
	Score            int
	DistToMainAnchor int // this allows us to prioritize remaps if they have the same score
}

type Remap struct {
	SequenceIndex    int
	Mm               []int
	MainAnchorRead   regionvector.Region
	MainAnchorGenome regionvector.Region
	LeftSections     []*RemapSection
	RightSections    []*RemapSection
	TotalLeftPaths   int
	TotalRightPaths  int
	LeftLength       int
	RightLength      int
}

type RemapSection struct {
	Mm              []int
	MatchedRead     []regionvector.Region
	MatchedGenome   []regionvector.Region
	MainAnchor      regionvector.Region
	TotalLeftPaths  int
	TotalRightPaths int
}

func fixPointRNARemap(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, genomeIndex *index.GenomeIndex, intronTree *datastructure.INTree) []*mapperutils.ReadMatchResult {
	// init list of remaps
	remaps := make([]*Remap, 0)

	// get main anchor of original map
	oriMainAnchor, _, mainAnchorIndex := readMatchResult.GetLargestAnchor(targetSeqIntronSet)
	oriReadMainAnchor := readMatchResult.MatchedRead.Regions[mainAnchorIndex]
	oriMmMainAnchor := extractMMofAnchor(oriReadMainAnchor, readMatchResult.MismatchesRead)

	if oriMainAnchor.Length()*10 < len(*read.Sequence)*3 {
		return nil // if the main anchor is super short, dont bother remapping
	}

	// get possible padding of l and r (we expect only one padding each)
	lPaddings, rPaddings := getPadding(oriMainAnchor, targetSeqIntronSet, intronTree)

	adjustedAnchors := make([]regionvector.Region, 0)
	adjustedReadAnchors := make([]regionvector.Region, 0)
	adjustedMms := make([][]int, 0)

	adjustedAnchors = append(adjustedAnchors, oriMainAnchor)
	adjustedReadAnchors = append(adjustedReadAnchors, oriReadMainAnchor)
	adjustedMms = append(adjustedMms, oriMmMainAnchor)

	if len(lPaddings) > 0 && len(rPaddings) > 0 {
		for l := range lPaddings {
			for r := range rPaddings {
				adjustedAnchor := oriMainAnchor.Copy()
				adjustedReadAnchor := oriReadMainAnchor.Copy()

				// try left pad
				appliedL := false
				if adjustedAnchor.Start+l < adjustedReadAnchor.End {
					adjustedAnchor.Start += l
					adjustedReadAnchor.Start += l
					appliedL = true
				}

				// try right pad
				appliedR := false
				if adjustedAnchor.End-r > adjustedAnchor.Start {
					adjustedAnchor.End -= r
					adjustedReadAnchor.End -= r
					appliedR = true
				}

				adjustedMM := extractMMofAnchor(adjustedReadAnchor, oriMmMainAnchor)

				// we need to account for case where applying l pad makes it impossible to also apply r pad

				// case 1: both paddings applied successfully
				if appliedL && appliedR {
					adjustedAnchors = append(adjustedAnchors, adjustedAnchor)
					adjustedReadAnchors = append(adjustedReadAnchors, adjustedReadAnchor)
					adjustedMms = append(adjustedMms, adjustedMM)
					continue
				}

				// case 2: only left applied and r pad could not be applied
				if appliedL && !appliedR {
					anchorL := adjustedAnchor.Copy()
					readL := adjustedReadAnchor.Copy()
					mmL := extractMMofAnchor(readL, oriMmMainAnchor)

					adjustedAnchors = append(adjustedAnchors, anchorL)
					adjustedReadAnchors = append(adjustedReadAnchors, readL)
					adjustedMms = append(adjustedMms, mmL)
				}

				// case 3: only right applied
				if appliedR && !appliedL {
					anchorR := oriMainAnchor.Copy()
					readR := oriReadMainAnchor.Copy()

					anchorR.End -= r
					readR.End -= r
					mmR := extractMMofAnchor(readR, oriMmMainAnchor)

					adjustedAnchors = append(adjustedAnchors, anchorR)
					adjustedReadAnchors = append(adjustedReadAnchors, readR)
					adjustedMms = append(adjustedMms, mmR)
				}
			}
		}
	} else if len(rPaddings) > 0 {
		for r := range rPaddings {
			adjustedAnchor := oriMainAnchor.Copy()
			adjustedReadAnchor := oriReadMainAnchor.Copy()
			if oriMainAnchor.End-r > oriMainAnchor.Start {
				adjustedAnchor.End -= r
				adjustedReadAnchor.End -= r
			}
			adjustedMM := extractMMofAnchor(adjustedReadAnchor, oriMmMainAnchor)

			adjustedAnchors = append(adjustedAnchors, adjustedAnchor)
			adjustedReadAnchors = append(adjustedReadAnchors, adjustedReadAnchor)
			adjustedMms = append(adjustedMms, adjustedMM)
		}
	} else if len(lPaddings) > 0 {
		for l := range lPaddings {
			adjustedAnchor := oriMainAnchor.Copy()
			adjustedReadAnchor := oriReadMainAnchor.Copy()
			if oriReadMainAnchor.Start+l < oriReadMainAnchor.End {
				adjustedAnchor.Start += l
				adjustedReadAnchor.Start += l
			}
			adjustedMM := extractMMofAnchor(adjustedReadAnchor, oriMmMainAnchor)

			adjustedAnchors = append(adjustedAnchors, adjustedAnchor)
			adjustedReadAnchors = append(adjustedReadAnchors, adjustedReadAnchor)
			adjustedMms = append(adjustedMms, adjustedMM)
		}
	}

	for i := 0; i < len(adjustedReadAnchors); i++ {
		remap := Remap{
			SequenceIndex:    readMatchResult.SequenceIndex,
			Mm:               adjustedMms[i],
			MainAnchorRead:   adjustedReadAnchors[i],
			MainAnchorGenome: adjustedAnchors[i],
			LeftSections:     make([]*RemapSection, 0),
			RightSections:    make([]*RemapSection, 0),
		}

		remaps = append(remaps, &remap)
	}

	for _, anchorToRemap := range remaps {
		// is left remap possible from mainAnchor?
		if anchorToRemap.MainAnchorRead.Start != 0 {
			missingBases := anchorToRemap.MainAnchorRead.Start
			startPos := anchorToRemap.MainAnchorGenome.Start
			leftPaths := targetSeqIntronSet.TranscriptomeGraph.FindPathsLeft(startPos, missingBases)
			// we can also add the previously mapped region to possible leftPaths
			// lets say a read is already mapped with two segments and misses some bases at the end
			// we also want to see how good of a remap we can get by using the already mapped short anchor(s) to the left of main anchor as a remap possibility
			if mainAnchorIndex != 0 {
				path := make([]regionvector.Region, 0)
				l := 0
				for k := mainAnchorIndex - 1; k >= 0; k-- {
					path = append(path, readMatchResult.MatchedGenome.Regions[k])
					l += readMatchResult.MatchedGenome.Regions[k].Length()
				}
				if l == missingBases {
					leftPaths = append(leftPaths, path)
				}
			}
			if len(leftPaths) != 0 {
				remapSectionsLeft := scoreLeftOptions(leftPaths, anchorToRemap.MainAnchorRead.Start, read.Sequence, genomeIndex.Sequences[anchorToRemap.SequenceIndex], anchorToRemap.MainAnchorGenome)
				anchorToRemap.TotalLeftPaths = len(leftPaths)
				anchorToRemap.LeftSections = remapSectionsLeft
			}
		}

		// is right remap possible from mainAnchor?
		if anchorToRemap.MainAnchorRead.End != len(*read.Sequence) {
			missingBases := len(*read.Sequence) - anchorToRemap.MainAnchorRead.End
			startPos := anchorToRemap.MainAnchorGenome.End
			rightPaths := targetSeqIntronSet.TranscriptomeGraph.FindPathsRight(startPos, missingBases)
			// we can also add the next mapped region to possible rightPaths
			// lets say a read is already mapped with two segments and misses some bases at the start
			// we also want to see how good of a remap we can get by using the already mapped short anchor(s) to the right of main anchor as a remap possibility
			if mainAnchorIndex != len(readMatchResult.MatchedGenome.Regions)-1 {
				path := make([]regionvector.Region, 0)
				l := 0
				for k := mainAnchorIndex + 1; k < len(readMatchResult.MatchedGenome.Regions); k++ {
					path = append(path, readMatchResult.MatchedGenome.Regions[k])
					l += readMatchResult.MatchedGenome.Regions[k].Length()
				}
				if l == missingBases {
					rightPaths = append(rightPaths, path)
				}
			}
			if len(rightPaths) != 0 {
				remapSectionsRight := scoreRightOptions(rightPaths, anchorToRemap.MainAnchorRead.End, read.Sequence, genomeIndex.Sequences[anchorToRemap.SequenceIndex], anchorToRemap.MainAnchorGenome)
				anchorToRemap.TotalRightPaths = len(rightPaths)
				anchorToRemap.RightSections = remapSectionsRight
			}
		}
	}
	readMatchResultRemaps := extractCandidates(remaps, len(*read.Sequence))
	return readMatchResultRemaps
}

func (r *RemapSection) reverseRegions() {
	// WARN:only call for left sections
	lenGenome := len(r.MatchedGenome)
	revGenome := make([]regionvector.Region, lenGenome)
	for i, reg := range r.MatchedGenome {
		revGenome[lenGenome-1-i] = reg
	}
	r.MatchedGenome = revGenome

	lenRead := len(r.MatchedRead)
	revRead := make([]regionvector.Region, lenRead)
	for i, reg := range r.MatchedRead {
		revRead[lenRead-1-i] = reg
	}
	r.MatchedRead = revRead

	lenMM := len(r.Mm)
	revMm := make([]int, lenMM)
	for i, mm := range r.Mm {
		revMm[lenMM-1-i] = mm
	}
	r.Mm = revMm
}

func extractCandidates(anchorRemaps []*Remap, readLength int) []*mapperutils.ReadMatchResult {
	// for each anchor, there can be several combinations. This means 1 Anchor -> N readMatchResults
	finalResults := make([]*mapperutils.ReadMatchResult, 0)

	for _, anchorRemap := range anchorRemaps {
		// get lowest score for left sections
		var leftCandidates []*RemapSection
		var rightCandidates []*RemapSection
		if len(anchorRemap.LeftSections) != 0 {
			leftCandidates = getBestCandidateSections(anchorRemap.LeftSections)
			// for these candidates we need to do an extra step since their regions and mms are in rev order
			for _, l := range leftCandidates {
				l.reverseRegions()
			}
		}
		if len(anchorRemap.RightSections) != 0 {
			rightCandidates = getBestCandidateSections(anchorRemap.RightSections)
		}

		if rightCandidates != nil && leftCandidates != nil {
			// create all pairwise combos :)
			for _, lCandidate := range leftCandidates {
				for _, rCandidate := range rightCandidates {
					alternativeReadMatchResult := &mapperutils.ReadMatchResult{
						SequenceIndex:       anchorRemap.SequenceIndex,
						MatchedRead:         &regionvector.RegionVector{},
						MatchedGenome:       &regionvector.RegionVector{},
						MismatchesRead:      make([]int, 0),
						IsFixPoint:          true,
						MainAnchorMM:        len(anchorRemap.Mm),
						MainAnchorLength:    anchorRemap.MainAnchorRead.Length(),
						LeftFixpointLength:  anchorRemap.MainAnchorRead.Start,
						RightFixpointLength: readLength - anchorRemap.MainAnchorRead.End,
						TotalLeftOptions:    anchorRemap.TotalLeftPaths,
						TotalRightOptions:   anchorRemap.TotalRightPaths,
						ValidLeftOptions:    len(leftCandidates),
						ValidRightOptions:   len(rightCandidates),
					}
					alternativeReadMatchResult.MismatchesRead = append(alternativeReadMatchResult.MismatchesRead, lCandidate.Mm...)
					alternativeReadMatchResult.MismatchesRead = append(alternativeReadMatchResult.MismatchesRead, anchorRemap.Mm...)
					alternativeReadMatchResult.MismatchesRead = append(alternativeReadMatchResult.MismatchesRead, rCandidate.Mm...)
					if uint8(float64(len(alternativeReadMatchResult.MismatchesRead))*100/float64(readLength)) > config.MaxMismatchPercentage() {
						continue
					} else {
						alternativeReadMatchResult.IncompleteMap = false
					}

					alternativeReadMatchResult.MatchedRead.Regions = append(alternativeReadMatchResult.MatchedRead.Regions, lCandidate.MatchedRead...)
					alternativeReadMatchResult.MatchedRead.Regions = append(alternativeReadMatchResult.MatchedRead.Regions, anchorRemap.MainAnchorRead)
					alternativeReadMatchResult.MatchedRead.Regions = append(alternativeReadMatchResult.MatchedRead.Regions, rCandidate.MatchedRead...)

					alternativeReadMatchResult.MatchedGenome.Regions = append(alternativeReadMatchResult.MatchedGenome.Regions, lCandidate.MatchedGenome...)
					alternativeReadMatchResult.MatchedGenome.Regions = append(alternativeReadMatchResult.MatchedGenome.Regions, anchorRemap.MainAnchorGenome)
					alternativeReadMatchResult.MatchedGenome.Regions = append(alternativeReadMatchResult.MatchedGenome.Regions, rCandidate.MatchedGenome...)

					finalResults = append(finalResults, alternativeReadMatchResult)
				}
			}
		} else if leftCandidates != nil {
			for _, lCandidate := range leftCandidates {
				alternativeReadMatchResult := &mapperutils.ReadMatchResult{
					SequenceIndex:       anchorRemap.SequenceIndex,
					MatchedRead:         &regionvector.RegionVector{},
					MatchedGenome:       &regionvector.RegionVector{},
					MismatchesRead:      make([]int, 0),
					IsFixPoint:          true,
					MainAnchorMM:        len(anchorRemap.Mm),
					MainAnchorLength:    anchorRemap.MainAnchorRead.Length(),
					LeftFixpointLength:  anchorRemap.MainAnchorRead.Start,
					RightFixpointLength: readLength - anchorRemap.MainAnchorRead.End,
					TotalLeftOptions:    anchorRemap.TotalLeftPaths,
					TotalRightOptions:   anchorRemap.TotalRightPaths,
					ValidLeftOptions:    len(leftCandidates),
					ValidRightOptions:   len(rightCandidates),
				}
				alternativeReadMatchResult.MismatchesRead = append(alternativeReadMatchResult.MismatchesRead, lCandidate.Mm...)
				alternativeReadMatchResult.MismatchesRead = append(alternativeReadMatchResult.MismatchesRead, anchorRemap.Mm...)

				if uint8(float64(len(alternativeReadMatchResult.MismatchesRead))*100/float64(readLength)) > config.MaxMismatchPercentage() {
					continue
				} else {
					alternativeReadMatchResult.IncompleteMap = false
				}

				alternativeReadMatchResult.MatchedRead.Regions = append(alternativeReadMatchResult.MatchedRead.Regions, lCandidate.MatchedRead...)
				alternativeReadMatchResult.MatchedRead.Regions = append(alternativeReadMatchResult.MatchedRead.Regions, anchorRemap.MainAnchorRead)

				alternativeReadMatchResult.MatchedGenome.Regions = append(alternativeReadMatchResult.MatchedGenome.Regions, lCandidate.MatchedGenome...)
				alternativeReadMatchResult.MatchedGenome.Regions = append(alternativeReadMatchResult.MatchedGenome.Regions, anchorRemap.MainAnchorGenome)

				finalResults = append(finalResults, alternativeReadMatchResult)
			}
		} else if rightCandidates != nil {
			for _, rCandidate := range rightCandidates {
				alternativeReadMatchResult := &mapperutils.ReadMatchResult{
					SequenceIndex:       anchorRemap.SequenceIndex,
					MatchedRead:         &regionvector.RegionVector{},
					MatchedGenome:       &regionvector.RegionVector{},
					MismatchesRead:      make([]int, 0),
					IsFixPoint:          true,
					MainAnchorMM:        len(anchorRemap.Mm),
					MainAnchorLength:    anchorRemap.MainAnchorRead.Length(),
					LeftFixpointLength:  anchorRemap.MainAnchorRead.Start,
					RightFixpointLength: readLength - anchorRemap.MainAnchorRead.End,
					TotalLeftOptions:    anchorRemap.TotalLeftPaths,
					TotalRightOptions:   anchorRemap.TotalRightPaths,
					ValidLeftOptions:    len(leftCandidates),
					ValidRightOptions:   len(rightCandidates),
				}
				alternativeReadMatchResult.MismatchesRead = append(alternativeReadMatchResult.MismatchesRead, anchorRemap.Mm...)
				alternativeReadMatchResult.MismatchesRead = append(alternativeReadMatchResult.MismatchesRead, rCandidate.Mm...)

				if uint8(float64(len(alternativeReadMatchResult.MismatchesRead))*100/float64(readLength)) > config.MaxMismatchPercentage() {
					continue
				} else {
					alternativeReadMatchResult.IncompleteMap = false
				}

				alternativeReadMatchResult.MatchedRead.Regions = append(alternativeReadMatchResult.MatchedRead.Regions, anchorRemap.MainAnchorRead)
				alternativeReadMatchResult.MatchedRead.Regions = append(alternativeReadMatchResult.MatchedRead.Regions, rCandidate.MatchedRead...)

				alternativeReadMatchResult.MatchedGenome.Regions = append(alternativeReadMatchResult.MatchedGenome.Regions, anchorRemap.MainAnchorGenome)
				alternativeReadMatchResult.MatchedGenome.Regions = append(alternativeReadMatchResult.MatchedGenome.Regions, rCandidate.MatchedGenome...)

				finalResults = append(finalResults, alternativeReadMatchResult)
			}
		}

	}
	return finalResults
}

func getBestCandidateSections(sections []*RemapSection) []*RemapSection {
	// determine min score and return all candidates with min score
	minScore := 10000000
	candidates := make([]*RemapSection, 0)
	for _, remapSection := range sections {
		curr := len(remapSection.Mm)
		if curr < minScore {
			minScore = curr
		}
	}
	for _, remapSection := range sections {
		curr := len(remapSection.Mm)
		if curr == minScore {
			candidates = append(candidates, remapSection)
		}
	}
	return candidates
}

func scoreRightOptions(rightPaths [][]regionvector.Region, startInRead int, readSeq *[]byte, refSeq *[]byte, mainAnchor regionvector.Region) []*RemapSection {
	remapSectionsRight := make([]*RemapSection, 0)
	lenReadSeq := len(*readSeq)
PathLoop:
	for _, path := range rightPaths {
		currSection := RemapSection{
			MainAnchor:    mainAnchor,
			MatchedGenome: path,
			MatchedRead:   make([]regionvector.Region, 0),
			Mm:            make([]int, 0),
		}
		lastStart := startInRead
		for _, genomicRegion := range path {
			readRegion := regionvector.Region{
				Start: lastStart,
				End:   lastStart + genomicRegion.Length(),
			}

			// score mm
			for refPos := genomicRegion.Start; refPos < genomicRegion.End; refPos++ {
				if refPos >= len(*refSeq) || lastStart >= len(*readSeq) || lastStart < 0 {
					continue PathLoop // dont append this path to remapSectionsRight since remap path is out of gene bounds
				}
				if (*refSeq)[refPos] != (*readSeq)[lastStart] && lastStart != startInRead {
					currSection.Mm = append(currSection.Mm, lastStart)
				}
				if uint8(float64(len(currSection.Mm))*100/float64(lenReadSeq)) >= config.MaxMismatchPercentage() {
					continue PathLoop // skip this path and don't append to remapSectionsLeft since mms already too many
				}
				lastStart++
			}

			lastStart = readRegion.End
			currSection.MatchedRead = append(currSection.MatchedRead, readRegion)
		}
		remapSectionsRight = append(remapSectionsRight, &currSection)
	}
	return remapSectionsRight
}

func scoreLeftOptions(leftPaths [][]regionvector.Region, startInRead int, readSeq *[]byte, refSeq *[]byte, mainAnchor regionvector.Region) []*RemapSection {
	remapSectionsLeft := make([]*RemapSection, 0)
	lenReadSeq := len(*readSeq)
PathLoop:
	for _, path := range leftPaths {
		// WARN: left pasths are inverted ;)
		currSection := RemapSection{
			MainAnchor:    mainAnchor,
			MatchedGenome: path,
			MatchedRead:   make([]regionvector.Region, 0),
			Mm:            make([]int, 0),
		}
		lastStart := startInRead
		for _, genomicRegion := range path {
			readRegion := regionvector.Region{
				Start: lastStart - genomicRegion.Length(),
				End:   lastStart,
			}

			// score mm
			for refPos := genomicRegion.End - 1; refPos >= genomicRegion.Start; refPos-- {
				if refPos < 0 || lastStart >= lenReadSeq || lastStart < 0 {
					continue PathLoop // skip this path and don't append to remapSectionsLeft since out of gene bounds
				}
				if (*refSeq)[refPos] != (*readSeq)[lastStart-1] && startInRead != lastStart {
					currSection.Mm = append(currSection.Mm, lastStart-1)
				}

				if uint8(float64(len(currSection.Mm))*100/float64(lenReadSeq)) >= config.MaxMismatchPercentage() {
					continue PathLoop // skip this path and don't append to remapSectionsLeft since mms already too many
				}

				lastStart--
			}

			lastStart = readRegion.Start
			currSection.MatchedRead = append(currSection.MatchedRead, readRegion)
		}
		remapSectionsLeft = append(remapSectionsLeft, &currSection)
	}
	return remapSectionsLeft
}

func getPadding(mainAnchor regionvector.Region, targetSeqIntronSet *regionvector.RegionSet, intronTree *datastructure.INTree) (map[int]struct{}, map[int]struct{}) { // l, r
	// overlappingIntrons := targetSeqIntronSet.GetIntersectingIntrons(mainAnchor)
	overlappingIntrons := make([]*regionvector.Intron, 0)
	leftIntersections := intronTree.Including(float64(mainAnchor.Start))
	rightIntersections := intronTree.Including(float64(mainAnchor.End))
	for _, inx := range leftIntersections {
		overlappingIntrons = append(overlappingIntrons, targetSeqIntronSet.Regions[inx])
	}
	for _, inx := range rightIntersections {
		overlappingIntrons = append(overlappingIntrons, targetSeqIntronSet.Regions[inx])
	}

	// try every combination of start coords
	rPadding := make(map[int]struct{})
	lPadding := make(map[int]struct{})
	for _, intron := range overlappingIntrons {
		if intron.LeftOverlap(mainAnchor) {
			padding := intron.End - mainAnchor.Start
			if utils.Abs(padding) < mainAnchor.Length() {
				lPadding[padding] = struct{}{}
			}
		}
		if intron.RightOverlap(mainAnchor) {
			padding := mainAnchor.End - intron.Start
			if utils.Abs(padding) < mainAnchor.Length() {
				rPadding[padding] = struct{}{}
			}
		}
	}
	return lPadding, rPadding
}

func extractMMofAnchor(anchor regionvector.Region, mms []int) []int {
	extracted := make([]int, 0)
	for _, mm := range mms {
		if mm >= anchor.Start && mm < anchor.End {
			extracted = append(extracted, mm)
		}
	}
	return extracted
}

func fillGaps(readMatchResult *mapperutils.ReadMatchResult, genomeIndex *index.GenomeIndex, read *fastq.Read) {
	// used to keep track of the read position for the next gap
	readGapPos := 0
	// returns the index of the first region after which a gap occurs (-1 if no gap)
	indexRegionBeforeGap := readMatchResult.MatchedRead.GetGapIndexAfterPos(readGapPos)

	// loop through all gaps in the read (-1 means there is no more gap)
	for indexRegionBeforeGap > -1 {

		gapRead, _ := readMatchResult.MatchedRead.GetGapAfterRegionIndex(indexRegionBeforeGap)
		gapGenome, _ := readMatchResult.MatchedGenome.GetGapAfterRegionIndex(indexRegionBeforeGap)
		gapStart := gapRead.Start
		gapEnd := gapRead.End
		gapGenomeStart := gapGenome.Start
		gapGenomeEnd := gapGenome.End

		if gapRead.Start != gapRead.End {
			// if we end up in here, it means our map looks like this
			// [0,97], [113, 150] -> mid block is missing

			// is there enough space in the genome gap to fill in the missing read portion
			if gapEnd-gapStart <= gapGenome.End-gapGenome.Start && gapRead.Length() > 0 {
				// is insert
				bestSplit := determineBestSplit(genomeIndex, read, readMatchResult.SequenceIndex, &gapRead, &gapGenome)

				if bestSplit == -1 {
					// this should not happen because a split should be found every time
					// even if the split has a bad score (many mismatches, no splice sites, etc)
					logrus.WithFields(logrus.Fields{
						"qname": read.Header,
					}).Fatal("no best split found")
				}

				// when bestSplit is 0 then there is nothing to be added to the left side of the gap
				if bestSplit > 0 {

					// add the split to the readMatchResult
					readMatchResult.MatchedRead.AddRegionNonOverlappingPanic(gapRead.Start, gapRead.Start+bestSplit)
					readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(gapGenome.Start, gapGenome.Start+bestSplit)

					// the read and genome sequences from the start of the gap to the best split (left)
					readByte := (*read.Sequence)[gapRead.Start : gapRead.Start+bestSplit]
					genomeByte := (*genomeIndex.Sequences[readMatchResult.SequenceIndex])[gapGenome.Start : gapGenome.Start+bestSplit]

					// add the mismatches to the readMatchResult
					for j := 0; j < bestSplit; j++ {
						// add the mismatche to the readMatchResult
						if readByte[j] != genomeByte[j] {
							readMatchResult.MismatchesRead = append(readMatchResult.MismatchesRead, gapRead.Start+j)
						}
					}
				}

				// when bestSplit is equal to the length of the gap then there is nothing
				// to be added to the right side of the gap
				if bestSplit < gapRead.Length() {

					// add the split to the readMatchResult
					readMatchResult.MatchedRead.AddRegionNonOverlappingPanic(gapRead.End-(gapRead.Length()-bestSplit), gapRead.End)
					readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(gapGenome.End-(gapRead.Length()-bestSplit), gapGenome.End)

					// the read and genome sequences from the best split to the end of the gap (right)
					readByte := (*read.Sequence)[gapRead.End-(gapRead.Length()-bestSplit) : gapRead.End]
					genomeByte := (*genomeIndex.Sequences[readMatchResult.SequenceIndex])[gapGenome.End-(gapRead.Length()-bestSplit) : gapGenome.End]

					// add the mismatches to the readMatchResult
					for j := 0; j < gapRead.Length()-bestSplit; j++ {
						// add the mismatche to the readMatchResult
						if readByte[j] != genomeByte[j] {
							// readMatchResult.MismatchesRead = append(readMatchResult.MismatchesRead, gapRead.End-(bestSplit-i))
							readMatchResult.MismatchesRead = append(readMatchResult.MismatchesRead, gapRead.Start+bestSplit+j) // NEW
						}
					}
					readMatchResult.NormalizeRegions()
					readMatchResult.IsGapFillOverflow = true
					readMatchResult.GapsFilledOverflow = append(readMatchResult.GapsFilledOverflow, gapGenome.Length())
				}
			} else {
				if gapGenome.Length() == 0 {
					// there's no gap at all in genome
					indexRegionBeforeGap = readMatchResult.MatchedRead.GetGapIndexAfterPos(gapRead.End + 1)
					continue
				}
				// gap in genome smaller that read gap, we can minimize mm in that region
				// genomeGapLen < readGapLen
				seqIndex := readMatchResult.SequenceIndex
				if gapRead.Length() <= 0 {
					// fmt.Println(readMatchResult.MatchedGenome.Regions)
					// fmt.Println(readMatchResult.MatchedRead.Regions)
					// fmt.Println(read.Header)
					return
				}

				lErrors := make([]int, gapGenome.Length()+1)
				rErrors := make([]int, gapGenome.Length()+1)

				lErrors[0] = 0
				rErrors[0] = 0

				for j := 1; j <= gapGenome.Length(); j++ {
					lErrors[j] = lErrors[j-1]
					if (*read.Sequence)[gapRead.Start+j-1] != (*genomeIndex.Sequences[seqIndex])[gapGenome.Start+j-1] {
						lErrors[j]++
					}

					rErrors[j] = rErrors[j-1]
					if (*read.Sequence)[gapRead.End-j] != (*genomeIndex.Sequences[seqIndex])[gapGenome.End-j] {
						rErrors[j]++
					}

				}

				minMM := gapGenome.Length() + 1
				minSplit := -1

				j := gapGenome.Length()
				for k := 0; k <= gapGenome.Length(); k++ {
					split := lErrors[k] + rErrors[j]
					if split < minMM {
						minMM = split
						minSplit = k
					}
					j--

				}
				rSplit := gapGenome.Length() - minSplit
				// lDonorSeq := string((*read.Sequence)[gapRead.Start : gapRead.Start+minSplit])
				// rDonorSeq := string((*read.Sequence)[gapRead.End-rSplit : gapRead.End])

				// close gap in genome
				readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(gapGenomeStart, gapGenomeEnd)

				// add left part of read
				if minSplit > 0 {
					readMatchResult.MatchedRead.AddRegionNonOverlappingPanic(gapRead.Start, gapRead.Start+minSplit)
				}

				// add right part of read
				if rSplit > 0 {
					readMatchResult.MatchedRead.AddRegionNonOverlappingPanic(gapRead.End-rSplit, gapRead.End)
				}
				readMatchResult.NormalizeRegions()
				readMatchResult.IsGapFill = true
				readMatchResult.GapsFilled = append(readMatchResult.GapsFilled, gapRead.Length())
			}
		}
		indexRegionBeforeGap = readMatchResult.MatchedRead.GetGapIndexAfterPos(gapRead.End + 1)
	}
}

func determineBestSplit(
	genomeIndex *index.GenomeIndex,
	read *fastq.Read,
	seqIndex int,
	gapRead *regionvector.Region,
	gapGenome *regionvector.Region,
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

		spliceSitePenalty, _ := utils.ScoreSpliceSites(donorSiteSeq[0], donorSiteSeq[1],
			acceptorSiteSeq[0], acceptorSiteSeq[1], lookOnPlusStrand)
		numMismatches += spliceSitePenalty

		if numMismatches <= minErrors {
			minErrors = numMismatches
			minSplit = i
		}
	}

	return minSplit
}
