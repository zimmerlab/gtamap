package secondpass

import (
	"fmt"
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/core/mapper/thirdpass"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/sirupsen/logrus"
)

func SecondpassMappingWorker(secondPassChan *SecondPassChannel, wgIncompleteMapping *sync.WaitGroup, annotationChan <-chan map[int]*mapperutils.TargetAnnotation, thirdPassChan *thirdpass.ThirdPassChannel, genomeIndex *index.GenomeIndex) {
	// in here we receive all non confident readpairs. Some have multiple maps for fw and rv and some are
	// not completely mapped yet. Confident readpairs are contained in confidentChan.
	defer wgIncompleteMapping.Done()
	logrus.Info("Started second pass")

	// main seq id -> seq annot
	var annotation map[int]*mapperutils.TargetAnnotation
	for annot := range annotationChan {
		// is only one object
		annotation = annot
	}
	fmt.Println(annotation)

	var wgRemap sync.WaitGroup

	for {
		// logrus.Info("About to receive from secondPassChan...")
		task, ok := secondPassChan.Receive()
		// logrus.Infof("Received from secondPassChan, ok=%v", ok)
		if !ok {
			logrus.Info("secondPassChan is closed, breaking from loop")
			break
		}

		wgRemap.Add(1)
		go func(t *mapperutils.ReadPairMatchResults) {
			defer wgRemap.Done()

			remapReadPair(task, annotation, genomeIndex)

			logrus.Debugf("Secondpass map: %s", task.ReadPair.ReadR1.Header)

			thirdPassChan.Send(&thirdpass.ThirdPassTask{
				ReadPairId: task.ReadPair.ReadR1.Header,
				TargetInfo: task,
			})
		}(task)

		// TODO: REMOVE DEBUG
		// for _, mapping := range task.Fw {
		// 	if mapping.IncompleteMap {
		// 		continue
		// 	}
		// 	debugout.GenerateAlignmentView(genomeIndex, *mapping, task.ReadPair.ReadR1)
		// }
		// for _, mapping := range task.Rv {
		// 	if mapping.IncompleteMap {
		// 		continue
		// 	}
		// 	debugout.GenerateAlignmentView(genomeIndex, *mapping, task.ReadPair.ReadR2)
		// }

	}
	wgRemap.Wait()
	logrus.Info("Done with second pass")
}

func remapReadPair(readPairMapping *mapperutils.ReadPairMatchResults, annotationMap map[int]*mapperutils.TargetAnnotation, genomeIndex *index.GenomeIndex) {
	for _, mapping := range readPairMapping.Fw {
		// here we megre intervals in the genomic regions for easier handling
		mapping.MatchedGenome.MergeAlignmentBlocks()
		mainSeqId := mapping.SequenceIndex / 2
		remapRead(mapping, annotationMap[mainSeqId], readPairMapping.ReadPair.ReadR1, genomeIndex)
	}
	for _, mapping := range readPairMapping.Rv {
		// here we megre intervals in the genomic regions for easier handling
		mapping.MatchedGenome.MergeAlignmentBlocks()
		mainSeqId := mapping.SequenceIndex / 2
		remapRead(mapping, annotationMap[mainSeqId], readPairMapping.ReadPair.ReadR2, genomeIndex)
	}
}

func remapRead(readMapping *mapperutils.ReadMatchResult, annotation *mapperutils.TargetAnnotation, read *fastq.Read, genomeIndex *index.GenomeIndex) {
	if read.Header == "2543591" {
		fmt.Println("s")
	}
	// remap IncompleteMap
	if readMapping.IncompleteMap {
		logrus.Debug("Have to remap readMatchResult due to: IncompleteMap")
		return
	}
	anchorGuidedRemap(readMapping, annotation.Introns[readMapping.SequenceIndex], read, genomeIndex)
}

// rightRemapAlignmentBlockFromPos remaps a given readSequenceToRemap from a given remapStart position.
func rightRemapAlignmentBlockFromPos(remapStart int, anchorRegion *regionvector.Region, targetSeqIntronSet *regionvector.RegionSet, readSequenceToRemap []byte, refSeq *[]byte) ([]*regionvector.Region, int, []int) {
	remap := make([]*regionvector.Region, 0)

	// get next intron
	nextIntronFromAnchor := targetSeqIntronSet.GetNextIntron(anchorRegion.Start)
	refPos := remapStart
	score := 0
	mm := make([]int, 0)

	for j := 0; j < len(readSequenceToRemap)-1; j++ {

		// jump to next exon in this case
		//        |anchor                  | here we jump
		// (READ) +++++++--------------****####-------------->
		//  (REF) +++++++****--------------------------####++>
		//                  | boundary
		if nextIntronFromAnchor != nil && refPos == nextIntronFromAnchor.Start-1 {
			remap = append(remap, &regionvector.Region{Start: remapStart, End: refPos + 1})

			if readSequenceToRemap[j] == (*refSeq)[refPos] {
				score++
			} else {
				mm = append(mm, refPos)
			}

			refPos = nextIntronFromAnchor.End
			remapStart = nextIntronFromAnchor.End
			continue
		}

		if readSequenceToRemap[j] == (*refSeq)[refPos] {
			score++
		} else {
			mm = append(mm, refPos)
		}
		refPos++
	}

	remap = append(remap, &regionvector.Region{Start: remapStart, End: refPos + 1})

	return remap, score, mm
}

// leftRemapAlignmentBlockFromPos remaps a given readSequenceToRemap from a given remapStart position.
func leftRemapAlignmentBlockFromPos(remapStart int, anchorRegion *regionvector.Region, targetSeqIntronSet *regionvector.RegionSet, readSequenceToRemap []byte, refSeq *[]byte) ([]*regionvector.Region, int, []int) {
	remap := make([]*regionvector.Region, 0)

	// get next intron
	nextIntronFromAnchor := targetSeqIntronSet.GetPrevIntron(anchorRegion.End)
	refPos := remapStart
	mm := make([]int, 0)
	explainedBases := 0
	score := 0

	for j := len(readSequenceToRemap) - 1; j >= 0; j-- {

		// jump to next exon in this case
		//                      | jump here
		// (READ) ----------####****-------+++++++++++---->
		//  (REF) -####----------------****+++++++++++---->
		//
		if nextIntronFromAnchor != nil && refPos == nextIntronFromAnchor.End {
			remap = append(remap, &regionvector.Region{Start: refPos, End: remapStart + 1})

			if readSequenceToRemap[j] != (*refSeq)[refPos] {
				mm = append(mm, refPos)
			} else {
				score++
			}

			refPos = nextIntronFromAnchor.Start - 1
			remapStart = nextIntronFromAnchor.Start - 1
			continue
		}

		if readSequenceToRemap[j] != (*refSeq)[refPos] {
			mm = append(mm, refPos)
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

func anchorGuidedRemap(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, genomeIndex *index.GenomeIndex) {
	if readMatchResult.MatchedGenome.HasGaps() {
		// get larges anchor in map
		mainAnchor, mainAnchorRank := readMatchResult.MatchedGenome.GetLargestAnchor()

		// Check if there are gaps to the right direction of the anchor
		if mainAnchorRank < len(readMatchResult.MatchedGenome.Regions)-1 {
			for mainAnchorRank != len(readMatchResult.MatchedGenome.Regions)-1 {
				//////////////////////////////
				//////   RIGHT REMAP   ///////
				//////////////////////////////
				// mainAnchor is left anchor extend right as far as possible
				// there are gaps to the right
				rGap := readMatchResult.MatchedGenome.GetGap(mainAnchorRank)
				if rGap == nil {
					break
				}
				weakAnchor := readMatchResult.MatchedGenome.Regions[mainAnchorRank+1]

				overlappingIntrons := targetSeqIntronSet.GetIntersectingIntrons(rGap)

				if len(overlappingIntrons) == 0 {
					// Case I (Deletion)
					// (READ) +++++++-------------+++++++
					// (REF) ++++++++++++++++++++++++++++++++
					logrus.WithFields(logrus.Fields{
						"Implied gap in read": rGap,
						"Length of Deletion":  rGap.Length(),
						"Read Blocks":         readMatchResult.MatchedGenome.Regions,
					}).Debug("Found deletion in read.")
					// Just log / use for report later
					// check if left or right diag can be extended
					mainAnchor = weakAnchor
					mainAnchorRank++
					continue
				}

				// (READ)      +++-----------------+++++++
				// (REF)  ++++++++--------+++------+++++++
				//                lIntron    rIntron
				// if only one intron is overlapping rIntron == lIntron
				lIntron := overlappingIntrons[0]                         // left most intron -> end of alignment block
				rIntron := overlappingIntrons[len(overlappingIntrons)-1] // right most intron -> start of next alignment block

				if lIntron.Start == rGap.Start && rIntron.End == rGap.End {
					// if gap already complies with intron boundaries,don't do anything
					// and skip to next anchor
					mainAnchor = weakAnchor
					mainAnchorRank++
					continue
				}

				logrus.WithFields(logrus.Fields{
					"Implied gap in read": rGap,
					"Read Blocks":         readMatchResult.MatchedGenome.Regions,
					"Overlapping Introns": overlappingIntrons,
				}).Debug("Found inconsistent read junction not following inferred intron boundaries")

				correctedL := lIntron.Start - rGap.Start
				correctedR := rIntron.End - rGap.End

				// if correctedL != correctedR -> remapDiagonal
				// this would equal CASE II right remap
				// (READ) +++++--------------------++++++
				// (REF)  ++++++++-------+++++++++++++++++++++++
				// or
				// (READ) +++++---------------++++
				// (REF)  ++++++++------------------------++++++
				if correctedR != correctedL {
					// do right cascade remap
					mmOptions := make(map[int][]int)
					remapOptions := make(map[int][]*regionvector.Region)
					bestOption := [1]int{-1}
					bestScore := -1
					regionReadStart := regionvector.GenomicCoordToReadCoord(readMatchResult.MatchedRead.GetFirstRegion().Start, weakAnchor.Start, readMatchResult.MatchedGenome.Regions)
					rightCascadeRemap(readMatchResult, targetSeqIntronSet, read, weakAnchor, mainAnchor, genomeIndex, &remapOptions, &mmOptions, &bestOption, bestScore, regionReadStart, -1)
					// is best option better than already existing map

					// add best matching region back to result
					for _, region := range remapOptions[bestOption[0]] {
						readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(region.Start, region.End)
					}

					mainAnchor = weakAnchor
					mainAnchorRank++
					continue
				}

				// Case III
				// (READ) ++++++++++-------+++++++++++
				// (READ) ++++++-------+++++++++++++++ (or)
				// (REF)  ++++++++-------++++++++++++
				// to
				// (READ) ++++++++-------+++++++++++++
				// (READ) ++++++++-------+++++++++++++
				// (REF)  ++++++++-------++++++++++++
				mainAnchor.End = mainAnchor.End + correctedL
				weakAnchor.Start = weakAnchor.Start + correctedR
				// go to next anchor
				mainAnchor = weakAnchor
				mainAnchorRank++
			}
		} else {
			for mainAnchorRank > 0 {
				//////////////////////////////
				//////   LEFT REMAP    ///////
				//////////////////////////////
				// mainAnchor is right anchor -> extend as far left as possible
				// there are gaps to the left
				lGap := readMatchResult.MatchedGenome.GetGap(mainAnchorRank - 1)
				if lGap == nil {
					break
				}
				weakAnchor := readMatchResult.MatchedGenome.Regions[mainAnchorRank-1]

				overlappingIntrons := targetSeqIntronSet.GetIntersectingIntrons(lGap)

				if len(overlappingIntrons) == 0 {
					// Case I (Deletion)
					//                            |mainAnchor
					// (READ) ++++----------------+++++++
					// (REF) ++++++++++++++++++++++++++++++++
					logrus.WithFields(logrus.Fields{
						"Implied gap in read": lGap,
						"Length of Deletion":  lGap.Length(),
						"Read Blocks":         readMatchResult.MatchedGenome.Regions,
					}).Debug("Found deletion in read.")
					// Just log / use for report later
					// check if left or right diag can be extended
					mainAnchor = weakAnchor
					mainAnchorRank--
					continue
				}

				//             |weakAnchor         |mainAnchor
				// (READ)      +++---------------------+++++++
				// (REF)  ++++++++---------+++---------+++++++
				//                |lIntron|   |rIntron|
				// if only one intron is overlapping rIntron == lIntron
				lIntron := overlappingIntrons[0]
				rIntron := overlappingIntrons[len(overlappingIntrons)-1]

				if lIntron.Start == lGap.Start && rIntron.End == lGap.End {
					// if gap already complies with intron boundaries, don't do anything and skip to next anchor
					mainAnchor = weakAnchor
					mainAnchorRank--
					continue
				}

				logrus.WithFields(logrus.Fields{
					"Implied gap in read": lGap,
					"Read Blocks":         readMatchResult.MatchedGenome.Regions,
					"Overlapping Introns": overlappingIntrons,
				}).Debug("Found inconsistent read junction not following inferred intron boundaries")

				correctedL := lIntron.Start - lGap.Start
				correctedR := rIntron.End - lGap.End

				// if correctedL != correctedR -> remapDiagonal
				// this would equal CASE II left remap
				//                              |mainAnchor
				// (READ) +++++-----------------++++++++++++
				// (REF)  ++++++++-------+++++++++++++++++++++++
				// or                                    |mainAnchor
				// (READ)            +++++---------------++++++++++
				// (REF)  ++++++++-----------------------++++++++++
				if correctedR != correctedL {
					// do left cascade remap

					remapOptions := make(map[int]*RemapOption)
					regionReadEnd := regionvector.GenomicCoordToReadCoord(readMatchResult.MatchedRead.GetFirstRegion().Start, weakAnchor.End, readMatchResult.MatchedGenome.Regions)

					leftCascadeRemap(readMatchResult, targetSeqIntronSet, read, weakAnchor, mainAnchor, genomeIndex, &remapOptions, regionReadEnd, -1)
					// is best option better than already existing map?
					// TODO

					// add best matching region back to result
					bestOption := -1
					localMax := -1
					for start, remapOption := range remapOptions {
						if remapOption.Score > localMax {
							localMax = remapOption.Score
							bestOption = start
						}
					}

					for _, region := range remapOptions[bestOption].RemapVector {
						readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(region.Start, region.End)
					}

					mainAnchor = weakAnchor
					mainAnchorRank--
					continue
				}

				// Case III
				// (READ) ++++++++++-------+++++++++++
				// (READ) ++++++-------+++++++++++++++ (or)
				// (REF)  ++++++++-------++++++++++++
				// to
				// (READ) ++++++++-------+++++++++++++
				// (READ) ++++++++-------+++++++++++++
				// (REF)  ++++++++-------++++++++++++
				mainAnchor.Start = mainAnchor.Start + correctedL
				weakAnchor.End = weakAnchor.End + correctedR
				// go to next anchor
				mainAnchor = weakAnchor
				mainAnchorRank--
			}
		}
	} else { // here we hanlde reads with no gaps
		// if there are no gaps we need to remap overhangs
		region := readMatchResult.MatchedGenome.Regions[0]
		overlappingIntrons := targetSeqIntronSet.GetIntersectingIntrons(region)
		if len(overlappingIntrons) == 0 {
			// solid map, no remap possible
			return
		}
		mapStart := region.Start
		mapEnd := region.End
		// either
		// +++++++++ (READ)
		// ----+++++(REF)
		// or
		// +++++++++ (READ)
		// +++++---- (REF)
		// or  (unlikely i think)
		// ++++++++++++ (READ)
		// ----+++++---(REF)
		lIntron := overlappingIntrons[0]
		rIntron := overlappingIntrons[len(overlappingIntrons)-1]

		// check if read start needs left remap
		if lIntron.End > mapStart {
			// remap
			// region.Start = lIntron.End
			//
			// regionReadStart := regionvector.GenomicCoordToReadCoord(readMatchResult.MatchedRead.GetFirstRegion().Start, mapStart, readMatchResult.MatchedGenome.Regions)
			// regionReadEnd := regionvector.GenomicCoordToReadCoord(readMatchResult.MatchedRead.GetFirstRegion().Start, region.Start, readMatchResult.MatchedGenome.Regions)
			//
			// readMatchResult.MatchedGenome.RemoveRegion(region.Start, lIntron.End)
			//
			// readSequenceToRemap := (*read.Sequence)[regionReadStart:regionReadEnd]
			// // refSeq := genomeIndex.Sequences[readMatchResult.SequenceIndex]
			// fmt.Println(read.Header)
			//
			// fmt.Println(string(readSequenceToRemap))

			// remapOptions, bestOption, _ := leftRemapAlignmentBlockFromPos(region.Start, region, targetSeqIntronSet, readSequenceToRemap, refSeq)
			//
			// for _, remapRegion := range remapOptions[bestOption] {
			// 	readMatchResult.MatchedRead.AddRegionNonOverlappingPanic(remapRegion.S)
			// }
		}

		// check if read end needs right remap
		if rIntron.Start < mapEnd {
			// remap
		}
	}
}

func rightCascadeRemap(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, regionToRemap *regionvector.Region, anchorRegion *regionvector.Region, genomeIndex *index.GenomeIndex, remapOptions *map[int][]*regionvector.Region, mmOptions *map[int][]int, bestOption *[1]int, bestScore int, initialRegionReadStart int, cascadeStart int) {
	// remove already used region we want to remap from result
	regionReadEnd := initialRegionReadStart + regionToRemap.Length()

	readMatchResult.MatchedGenome.RemoveRegion(regionToRemap.Start, regionToRemap.End)

	readSequenceToRemap := (*read.Sequence)[initialRegionReadStart:regionReadEnd]
	refSeq := genomeIndex.Sequences[readMatchResult.SequenceIndex]

	// get next intron
	intronRight := targetSeqIntronSet.GetNextIntron(anchorRegion.Start)
	if intronRight != nil {
		// first check if end of anchorRegion is at intron boundary
		//        |anchor
		// (READ) +++++++--------------++---------->
		//  (REF) +++++++++++---------------------->
		//                   | boundary
		if anchorRegion.End < intronRight.Start {
			remapPos := anchorRegion.End
			regions, score, mm := rightRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
			// here I need to use the original index if we are in a cascade, we are in a cascade
			if cascadeStart == anchorRegion.Start {
				(*mmOptions)[anchorRegion.Start] = append((*mmOptions)[anchorRegion.Start], mm...)
				(*remapOptions)[anchorRegion.Start] = append((*remapOptions)[anchorRegion.Start], regions...)
				// TODO Halndel score and mm
			} else {
				(*mmOptions)[remapPos] = append((*mmOptions)[remapPos], mm...)
				(*remapOptions)[remapPos] = append((*remapOptions)[remapPos], regions...)
			}
			if bestScore < score {
				bestOption[0] = remapPos
				bestScore = score
			}
		}
		// check for potential overhangs into neighbor intron of anchor seed
		//(READ) ++++++++++*-----######---------->
		//(REF) +++++++++++---------------*######>
		if anchorRegion.End > intronRight.Start {
			padding := anchorRegion.End - intronRight.Start

			// add padding to region to remap
			// before:               |remap|
			// READ: ++++++++++*-----#######--------->
			//  REF: ++++++++++---------------*######>
			// after:          |   remap   |
			// READ: ++++++++++*-----#######--------->
			//  REF: ++++++++++---------------*######>
			initialRegionReadStart = initialRegionReadStart - padding
			readMatchResult.MatchedGenome.RemoveRegion(anchorRegion.End-padding, anchorRegion.End)
			anchorRegion.End = anchorRegion.End - padding
			readSequenceToRemap = (*read.Sequence)[initialRegionReadStart:regionReadEnd]
		}

		// remap for remaining exons (in right direction)
		for i := intronRight.Rank; i < len(targetSeqIntronSet.Regions); i++ {
			refStart := targetSeqIntronSet.Regions[i].End
			regions, score, mm := rightRemapAlignmentBlockFromPos(refStart, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
			(*mmOptions)[refStart] = append((*mmOptions)[refStart], mm...)
			(*remapOptions)[refStart] = append((*remapOptions)[refStart], regions...)
			if bestScore < score {
				bestOption[0] = refStart
				bestScore = score
			}
		}
	} else {
		remapPos := anchorRegion.End
		regions, score, mm := rightRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
		(*mmOptions)[remapPos] = append((*mmOptions)[remapPos], mm...)
		(*remapOptions)[remapPos] = append((*remapOptions)[remapPos], regions...)
		if bestScore < score {
			bestOption[0] = remapPos
			bestScore = score
		}
	}

	// now we check if the region we just remapped skips any other region downstream of the map
	for _, region := range readMatchResult.MatchedGenome.Regions {
		// skip region we just remapped
		if region.Start == (*remapOptions)[bestOption[0]][0].Start {
			continue
		}
		// skip all regions left to anchor region since we are in a right remap
		if region.Start <= anchorRegion.Start {
			continue
		}
		if region.Start < bestOption[0] {
			fmt.Println("SKIPPED REGION DURING RIGHT REMAP -> CASCADING RIGHT")
			// as soon as we skip any region while remapping to the right, we have to also remap that region into the same direction
			// lAnchor is the last mapped region to the right of our remap
			anchorRegion = (*remapOptions)[bestOption[0]][len((*remapOptions)[bestOption[0]])-1]
			rightCascadeRemap(readMatchResult, targetSeqIntronSet, read, region, anchorRegion, genomeIndex, remapOptions, mmOptions, bestOption, bestScore, regionReadEnd, anchorRegion.Start)
		}
	}
}

func leftCascadeRemap(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, regionToRemap *regionvector.Region, anchorRegion *regionvector.Region, genomeIndex *index.GenomeIndex, remapOptions *map[int]*RemapOption, initialRegionReadEnd int, cascadeStart int) {
	regionReadStart := initialRegionReadEnd - regionToRemap.Length()

	// remove already used region we want to remap from result
	readMatchResult.MatchedGenome.RemoveRegion(regionToRemap.Start, regionToRemap.End)

	readSequenceToRemap := (*read.Sequence)[regionReadStart:initialRegionReadEnd]
	refSeq := genomeIndex.Sequences[readMatchResult.SequenceIndex]

	// get next intron
	intronLeft := targetSeqIntronSet.GetPrevIntron(anchorRegion.Start)
	if intronLeft != nil {
		if anchorRegion.Start > intronLeft.End {
			// check if start of anchorRegion is at intron boundary
			//                                 |anchor
			// (READ) -----++------------------+++++++->
			//  (REF) -----------------+++++++++++++++->
			//           nextIntron    |boundary
			remapPos := anchorRegion.Start - 1
			regions, score, mm := leftRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
			if cascadeStart != -1 { // are we in a cascade?
				(*remapOptions)[cascadeStart].UpdateMM(mm)
				(*remapOptions)[cascadeStart].UpdateRegions(regions)
				(*remapOptions)[cascadeStart].UpdateScore(score)

			} else {
				(*remapOptions)[remapPos] = &RemapOption{
					RemapStart:  remapPos,
					RemapVector: regions,
					MisMatches:  mm,
					Score:       score,
				}
			}
		} else if anchorRegion.Start < intronLeft.Start {
			// check for potential overhangs into neighbor intron of anchor seed
			padding := intronLeft.Start - anchorRegion.Start
			// add padding to region to remap
			// before:     |remap|
			//(READ) -------####------*+++++---------->
			//(REF)  ####*-------------+++++---------->
			// after:       |  remap  |
			//(READ) -------####------*+++++---------->
			//(REF)  ####*-------------+++++---------->
			initialRegionReadEnd = initialRegionReadEnd + padding
			readMatchResult.MatchedGenome.RemoveRegion(anchorRegion.Start, anchorRegion.Start+padding)
			anchorRegion.Start = anchorRegion.Start + padding
			readSequenceToRemap = (*read.Sequence)[regionReadStart:initialRegionReadEnd]
		}

		// remap for remaining exons (in left direction)
		for i := intronLeft.Rank; i >= 0; i-- {
			refStart := targetSeqIntronSet.Regions[i].Start - 1
			regions, score, mm := leftRemapAlignmentBlockFromPos(refStart, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
			if cascadeStart != -1 { // are we in a cascade
				(*remapOptions)[cascadeStart].UpdateMM(mm)
				(*remapOptions)[cascadeStart].UpdateRegions(regions)
				(*remapOptions)[cascadeStart].UpdateScore(score)
			} else {
				(*remapOptions)[refStart] = &RemapOption{
					RemapStart:  refStart,
					RemapVector: regions,
					MisMatches:  mm,
					Score:       score,
				}
			}
		}
	} else {
		remapPos := anchorRegion.Start - 1
		regions, score, mm := leftRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)

		if cascadeStart != -1 {
			(*remapOptions)[cascadeStart].UpdateMM(mm)
			(*remapOptions)[cascadeStart].UpdateRegions(regions)
			(*remapOptions)[cascadeStart].UpdateScore(score)
		} else {
			(*remapOptions)[remapPos] = &RemapOption{
				RemapStart:  remapPos,
				RemapVector: regions,
				MisMatches:  mm,
				Score:       score,
			}
		}
	}

	localMax := -1
	bestOption := -1
	for start, remapOption := range *remapOptions {
		if remapOption.Score > localMax {
			localMax = remapOption.Score
			bestOption = start
		}
	}

	// now we check if the region we just remapped skips any other region downstream of the map
	// since this is a left remap we go throug hregions in reverse order
	for i := len(readMatchResult.MatchedGenome.Regions) - 1; i >= 0; i-- {
		region := readMatchResult.MatchedGenome.Regions[i]
		// skip region we just remapped
		if region.Start == (*remapOptions)[bestOption].RemapStart {
			continue
		}
		// skip all regions right to anchor region since we are in a left remap
		if region.Start >= anchorRegion.Start {
			continue
		}
		if region.Start > bestOption {
			fmt.Println("SKIPPED REGION DURING LEFT REMAP -> CASCADING LEFT")
			fmt.Println(read.Header)

			// as soon as we skip any region while remapping to the left, we have to also remap that region into the same direction
			// lAnchor is the last mapped region to the right of our remap
			anchorRegion = (*remapOptions)[bestOption].RemapVector[0]
			// start cascade
			if cascadeStart == -1 {
				leftCascadeRemap(readMatchResult, targetSeqIntronSet, read, region, anchorRegion, genomeIndex, remapOptions, regionReadStart, anchorRegion.End-1)
			} else { // continue cascade
				leftCascadeRemap(readMatchResult, targetSeqIntronSet, read, region, anchorRegion, genomeIndex, remapOptions, regionReadStart, cascadeStart)
			}
		}
	}
}

type RemapOption struct {
	RemapStart  int
	RemapVector []*regionvector.Region
	MisMatches  []int
	Score       int
}

func (r *RemapOption) UpdateMM(mm []int) {
	if r.MisMatches == nil {
		r.MisMatches = mm
	}
	r.MisMatches = append(r.MisMatches, mm...)
}

func (r *RemapOption) UpdateRegions(remap []*regionvector.Region) {
	r.RemapVector = append(r.RemapVector, remap...)
}

func (r *RemapOption) UpdateScore(score int) {
	r.Score += score
}
