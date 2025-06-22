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

			// remapReadPair(task, annotation, genomeIndex)

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
	if readPairMapping.ReadPair.ReadR1.Header == "2739268" {
		fmt.Println("s")
	}
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
					remapOptions := rightRemap(readMatchResult, targetSeqIntronSet, read, weakAnchor.Start, mainAnchor, genomeIndex)

					// add best matching region back to result
					// we want the remap with the higest score and if there are remaps with same score, we start with
					bestOption := -1
					localMax := -1
					distToAnchor := 100000 // init with large value
					for start, remapOption := range remapOptions {
						if remapOption.Score > localMax {
							localMax = remapOption.Score
							bestOption = start
						} else if remapOption.Score == localMax {
							dist := remapOption.Rank - mainAnchorRank
							if dist < distToAnchor {
								fmt.Println("same")
								distToAnchor = dist
								bestOption = start
							}
						}
					}

					for _, region := range remapOptions[bestOption].RemapVector {
						readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(region.Start, region.End)
					}

					// maybe break?
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
					remapOptions := leftRemap(readMatchResult, targetSeqIntronSet, read, weakAnchor.End, mainAnchor, genomeIndex)
					// is best option better than already existing map?

					// add best matching region back to result
					bestOption := -1
					localMax := -1
					distToAnchor := 100000
					for start, remapOption := range remapOptions {
						if remapOption.Score > localMax {
							localMax = remapOption.Score
							bestOption = start
						} else if remapOption.Score == localMax {
							dist := mainAnchorRank - remapOption.Rank
							if dist < distToAnchor {
								fmt.Println("same")

								distToAnchor = dist
								bestOption = start
							}
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

func rightRemap(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, remapStart int, anchorRegion *regionvector.Region, genomeIndex *index.GenomeIndex) map[int]*RemapOption {
	remapOptions := make(map[int]*RemapOption)

	regionReadStart := regionvector.GenomicCoordToReadCoord(readMatchResult.MatchedRead.GetFirstRegion().Start, remapStart, readMatchResult.MatchedGenome.Regions)
	regionReadEnd := len((*read.Sequence)) // remap from remap start to end of read

	remapEnd := readMatchResult.MatchedGenome.GetLastRegion().End
	// remove all regions after remap start
	readMatchResult.MatchedGenome.RemoveRegion(remapStart, remapEnd)

	readSequenceToRemap := (*read.Sequence)[regionReadStart:regionReadEnd]
	refSeq := genomeIndex.Sequences[readMatchResult.SequenceIndex]

	// get next intron
	intronRight := targetSeqIntronSet.GetNextIntron(anchorRegion.Start)
	if intronRight != nil {
		if anchorRegion.End < intronRight.Start {
			// first check if end of anchorRegion is at intron boundary
			// if now we need to also remap from anchorRegion.End
			//               |anchorRegion.End
			// (READ) +++++++--------------++---------->
			//  (REF) +++++++++++------------------++++>
			//                   | intronRight.Start
			remapPos := anchorRegion.End
			regions, score, mm := rightRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
			remapOptions[remapPos] = &RemapOption{
				RemapStart:  remapPos,
				RemapVector: regions,
				MisMatches:  mm,
				Score:       score,
				Rank:        intronRight.Rank,
			}
		}

		// check for padding
		if anchorRegion.End > intronRight.Start {
			// check for potential overhangs into neighbor intron of anchor seed
			// before:               |remap|
			// READ: ++++++++++*-----#######--------->
			//  REF: ++++++++++---------------*######>
			padding := anchorRegion.End - intronRight.Start
			// after:          |   remap   |
			// READ: ++++++++++*-----#######--------->
			//  REF: ++++++++++---------------*######>
			regionReadStart = regionReadStart - padding
			readMatchResult.MatchedGenome.RemoveRegion(anchorRegion.End-padding, anchorRegion.End)
			anchorRegion.End = anchorRegion.End - padding
			readSequenceToRemap = (*read.Sequence)[regionReadStart:regionReadEnd]
		}

		// remap for all remaining exons (in right direction)
		for i := intronRight.Rank; i < len(targetSeqIntronSet.Regions); i++ {
			remapPos := targetSeqIntronSet.Regions[i].End
			regions, score, mm := rightRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
			remapOptions[remapPos] = &RemapOption{
				RemapStart:  remapPos,
				RemapVector: regions,
				MisMatches:  mm,
				Score:       score,
				Rank:        intronRight.Rank,
			}
		}
	} else {
		// if there are no introns left in right direction, start remap from anchorregion end
		remapPos := anchorRegion.End
		regions, score, mm := rightRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
		remapOptions[remapPos] = &RemapOption{
			RemapStart:  remapPos,
			RemapVector: regions,
			MisMatches:  mm,
			Score:       score,
			Rank:        targetSeqIntronSet.Regions[len(targetSeqIntronSet.Regions)-1].Rank + 1, // remap is after last intron
		}
	}

	return remapOptions
}

func leftRemap(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, remapStart int, anchorRegion *regionvector.Region, genomeIndex *index.GenomeIndex) map[int]*RemapOption {
	remapOptions := make(map[int]*RemapOption)

	regionReadStart := 0 // remap entire read from o to remapPos
	regionReadEnd := regionvector.GenomicCoordToReadCoord(readMatchResult.MatchedRead.GetFirstRegion().Start, remapStart, readMatchResult.MatchedGenome.Regions)

	// remove already used region we want to remap from result
	remapEnd := readMatchResult.MatchedGenome.GetFirstRegion().Start
	readMatchResult.MatchedGenome.RemoveRegion(remapEnd, remapStart)

	readSequenceToRemap := (*read.Sequence)[regionReadStart:regionReadEnd]
	refSeq := genomeIndex.Sequences[readMatchResult.SequenceIndex]

	// get next intron
	intronLeft := targetSeqIntronSet.GetPrevIntron(anchorRegion.Start)
	if intronLeft != nil {
		if anchorRegion.Start > intronLeft.End {
			// check if start of anchorRegion is at intron boundary
			//                                 |anchor
			// (READ) -----++------------------+++++++->
			//  (REF) +++--------------+++++++++++++++->
			//                        |intronLeft.End
			remapPos := anchorRegion.Start - 1
			regions, score, mm := leftRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
			remapOptions[remapPos] = &RemapOption{
				RemapStart:  remapPos,
				RemapVector: regions,
				MisMatches:  mm,
				Score:       score,
				Rank:        intronLeft.Rank,
			}
		}

		// check for padding
		if anchorRegion.Start < intronLeft.Start {
			// OR
			// check for potential overhangs into neighbor intron of anchor seed
			// before:     |remap|
			//(READ) -------####------*+++++---------->
			//(REF)  ####*-------------+++++---------->
			padding := intronLeft.Start - anchorRegion.Start
			// after:       |  remap  |
			//(READ) -------####------*+++++---------->
			//(REF)  ####*-------------+++++---------->
			regionReadEnd = regionReadEnd + padding
			readMatchResult.MatchedGenome.RemoveRegion(anchorRegion.Start, anchorRegion.Start+padding)
			anchorRegion.Start = anchorRegion.Start + padding
			readSequenceToRemap = (*read.Sequence)[regionReadStart:regionReadEnd]
		}

		// remap for remaining exons (in left direction)
		for i := intronLeft.Rank; i >= 0; i-- {
			remapPos := targetSeqIntronSet.Regions[i].Start - 1
			regions, score, mm := leftRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
			remapOptions[remapPos] = &RemapOption{
				RemapStart:  remapPos,
				RemapVector: regions,
				MisMatches:  mm,
				Score:       score,
				Rank:        intronLeft.Rank,
			}
		}
	} else {
		// if there are no introns left in left direction, start remap from anchorregion start - 1
		remapPos := anchorRegion.Start - 1
		regions, score, mm := leftRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
		remapOptions[remapPos] = &RemapOption{
			RemapStart:  remapPos,
			RemapVector: regions,
			MisMatches:  mm,
			Score:       score,
			Rank:        -1, // remap is before first intron
		}
	}
	return remapOptions
}

type RemapOption struct {
	RemapStart  int
	RemapVector []*regionvector.Region
	MisMatches  []int
	Score       int
	Rank        int // this allows us to priotitize remaps if they have the same score. The rank is the intron rank of the corresponding intron
}

func (r *RemapOption) UpdateRemap(remap []*regionvector.Region, mm []int, score int) {
	r.RemapVector = append(r.RemapVector, remap...)
	r.MisMatches = append(r.MisMatches, mm...)
	r.Score += score
}
