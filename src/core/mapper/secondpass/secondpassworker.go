package secondpass

import (
	"sync"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/core/mapper/thirdpass"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
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

	var wgRemap sync.WaitGroup
	logrus.Info("Started second pass")

	for {
		// logrus.Info("About to receive from secondPassChan...")
		task, ok := secondPassChan.Receive()

		// logrus.Infof("Received from secondPassChan, ok=%v", ok)
		if !ok {
			logrus.Info("Done with second pass")
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
}

func remapReadPair(readPairMapping *mapperutils.ReadPairMatchResults, annotationMap map[int]*mapperutils.TargetAnnotation, genomeIndex *index.GenomeIndex) {
	for _, mapping := range readPairMapping.Fw {
		mainSeqId := mapping.SequenceIndex / 2
		remapRead(mapping, annotationMap[mainSeqId], readPairMapping.ReadPair.ReadR1, genomeIndex)
	}
	for _, mapping := range readPairMapping.Rv {
		mainSeqId := mapping.SequenceIndex / 2
		remapRead(mapping, annotationMap[mainSeqId], readPairMapping.ReadPair.ReadR2, genomeIndex)
	}

	// TODO: make mappings uniq
}

func remapRead(readMapping *mapperutils.ReadMatchResult, annotation *mapperutils.TargetAnnotation, read *fastq.Read, genomeIndex *index.GenomeIndex) {
	// remap IncompleteMap
	if len(annotation.Introns[0].Regions) == 0 {
		logrus.Debug("Second pass remap not possible because no introns were inferred.")
		return
	}
	if readMapping.IncompleteMap || readMapping.MatchedGenome.Length() != len(*read.Sequence) {
		logrus.Debug("Have to remap readMatchResult due to: IncompleteMap")
		incomplRemap(readMapping, annotation.Introns[readMapping.SequenceIndex], read, genomeIndex)
		return
	} else {
		anchorGuidedRemap(readMapping, annotation.Introns[readMapping.SequenceIndex], read, genomeIndex)
	}
}

// rightRemapAlignmentBlockFromPos remaps a given readSequenceToRemap from a given remapStart position.
func rightRemapAlignmentBlockFromPos(remapStart int, anchorRegion *regionvector.Region, targetSeqIntronSet *regionvector.RegionSet, readSequenceToRemap []byte, refSeq *[]byte, startInRead int) ([]*regionvector.Region, int, []int) {
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
			remap = append(remap, &regionvector.Region{Start: remapStart, End: refPos})

			//if readSequenceToRemap[j] == (*refSeq)[refPos] {
			//	score++
			//} else {
			//	mm = append(mm, startInRead+j)
			//}

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
func leftRemapAlignmentBlockFromPos(remapStart int, anchorRegion *regionvector.Region, targetSeqIntronSet *regionvector.RegionSet, readSequenceToRemap []byte, refSeq *[]byte) ([]*regionvector.Region, int, []int) {
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
			remap = append(remap, &regionvector.Region{Start: refPos, End: remapStart + 1})

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

func anchorGuidedRemap(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, genomeIndex *index.GenomeIndex) {
	readMatchResult.NormalizeRegions()
	beforeRemap := readMatchResult.Copy()

	if readMatchResult.MatchedGenome.HasGaps() {
		// get larges anchor in map for right remap
		mainAnchor, mainAnchorRank, mainAnchorIndex := readMatchResult.GetLargestAnchor(targetSeqIntronSet)

		// Check if there are gaps to the right direction of the anchor
		// if already corrected the mainAnchor to the right, we don't need to perform a right remap again
		if mainAnchorIndex < len(readMatchResult.MatchedGenome.Regions)-1 {
			for mainAnchorIndex != len(readMatchResult.MatchedGenome.Regions)-1 {
				//////////////////////////////
				//////   RIGHT REMAP   ///////
				//////////////////////////////
				// mainAnchor is left anchor extend right as far as possible
				// there are gaps to the right
				rGap := readMatchResult.MatchedGenome.GetGap(mainAnchorIndex)
				if rGap == nil {
					break
				}
				weakAnchor := readMatchResult.MatchedGenome.Regions[mainAnchorIndex+1]

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
					// startInRead := regionvector.GenomicCoordToReadCoord(readMatchResult.MatchedRead.GetFirstRegion().Start, weakAnchor.Start, readMatchResult.MatchedGenome.Regions)
					startInRead := readMatchResult.MatchedRead.Regions[mainAnchorIndex+1].Start
					remapOptions := rightRemap(readMatchResult, targetSeqIntronSet, read, weakAnchor, mainAnchor, genomeIndex, mainAnchorRank, startInRead)
					if remapOptions != nil {
						bestOption := chooseBestRemap(remapOptions)

						if len(remapOptions[bestOption].MismatchesRead) < len(readMatchResult.MismatchesRead) {
							// remove weakAnchor and everything right to it
							readMatchResult.MatchedGenome.RemoveRegion(weakAnchor.Start, readMatchResult.MatchedGenome.GetLastRegion().End)
							// update map
							for _, region := range remapOptions[bestOption].RemapVector {
								readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(region.Start, region.End)
							}
							// update mm
							readMatchResult.MismatchesRead = remapOptions[bestOption].MismatchesRead
						}
					}

					mainAnchor = weakAnchor
					mainAnchorIndex++
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
					mainAnchorIndex++
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

					// startInRead := regionvector.GenomicCoordToReadCoord(readMatchResult.MatchedRead.GetFirstRegion().Start, weakAnchor.Start, readMatchResult.MatchedGenome.Regions)
					startInRead := readMatchResult.MatchedRead.Regions[mainAnchorIndex+1].Start

					// we check for padding if part of the main anchor is reaching into intron
					if mainAnchor.End > lIntron.Start && mainAnchor.Start < lIntron.Start {
						// check for potential overhangs into neighbor intron of anchor seed
						// before:               |remap|
						// READ: ++++++++++*-----#######--------->
						//  REF: ++++++++++---------------*######>
						//                 |lIntron
						padding := mainAnchor.End - lIntron.Start
						// after:          |   remap   |
						// READ: ++++++++++*-----#######--------->
						//  REF: ++++++++++---------------*######>
						mainAnchor.End -= padding
						weakAnchor.Start -= padding
						startInRead -= padding
					}
					remapOptions := rightRemap(readMatchResult, targetSeqIntronSet, read, weakAnchor, mainAnchor, genomeIndex, mainAnchorRank, startInRead)

					if remapOptions != nil {
						// add best matching region back to result
						// we want the remap with the highest score and if there are remaps with same score, we start with
						bestOption := chooseBestRemap(remapOptions)
						if len(remapOptions[bestOption].MismatchesRead) < len(readMatchResult.MismatchesRead) {

							// remove weakAnchor and everything right to it
							readMatchResult.MatchedGenome.RemoveRegion(weakAnchor.Start, readMatchResult.MatchedGenome.GetLastRegion().End)

							for _, region := range remapOptions[bestOption].RemapVector {
								readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(region.Start, region.End)
							}

							// update mm in readMatchResult
							readMatchResult.MismatchesRead = remapOptions[bestOption].MismatchesRead
						}
					}
				} else {
					// Case III
					// (READ) ++++++++++-------+++++++++++
					// (READ) ++++++-------+++++++++++++++ (or)
					// (REF)  ++++++++-------++++++++++++
					// to
					// (READ) ++++++++-------+++++++++++++
					// (READ) ++++++++-------+++++++++++++
					// (REF)  ++++++++-------++++++++++++
					if mainAnchor.Length() > abs(correctedL) && weakAnchor.Length() > abs(correctedR) {
						mainAnchor.End = mainAnchor.End + correctedL
						weakAnchor.Start = weakAnchor.Start + correctedR
					}
				}
				// go to next anchor
				mainAnchor = weakAnchor
				mainAnchorIndex++
			}
		}

		// get larges anchor in map for left remap
		mainAnchor, mainAnchorRank, mainAnchorIndex = readMatchResult.GetLargestAnchor(targetSeqIntronSet)

		// if already corrected the mainAnchor to the left, we don't need to perform a left remap again
		if mainAnchorIndex > 0 {
			for mainAnchorIndex > 0 {
				//////////////////////////////
				//////   LEFT REMAP    ///////
				//////////////////////////////
				// mainAnchor is right anchor -> extend as far left as possible
				// there are gaps to the left
				lGap := readMatchResult.MatchedGenome.GetGap(mainAnchorIndex - 1)
				if lGap == nil {
					break
				}
				weakAnchor := readMatchResult.MatchedGenome.Regions[mainAnchorIndex-1]

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
					regionReadEnd, err := regionvector.GenomicCoordToReadCoord(0, weakAnchor.End, readMatchResult.MatchedGenome.Regions)
					if err != nil {
						logrus.Errorf("Error while converting genomic coord to read coord in read %s", read.Header)
						logrus.Fatal(err)
					}
					remapOptions := leftRemap(readMatchResult, targetSeqIntronSet, read, weakAnchor, mainAnchor, genomeIndex, mainAnchorRank, regionReadEnd)
					if remapOptions != nil {
						bestOption := chooseBestRemap(remapOptions)

						if len(remapOptions[bestOption].MismatchesRead) < len(readMatchResult.MismatchesRead) {
							// remove weakAnchor and everything left from it
							readMatchResult.MatchedGenome.RemoveRegion(readMatchResult.MatchedGenome.GetFirstRegion().Start, weakAnchor.End)
							for _, region := range remapOptions[bestOption].RemapVector {
								readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(region.Start, region.End)
							}

							readMatchResult.MismatchesRead = remapOptions[bestOption].MismatchesRead
						}
					}

					mainAnchor = weakAnchor
					mainAnchorIndex--
					continue
				}

				//             |weakAnchor             |mainAnchor
				// (READ)      +++---------------------+++++++
				// (REF)  ++++++++---------+++---------+++++++
				//                |lIntron|   |rIntron|
				// if only one intron is overlapping rIntron == lIntron
				lIntron := overlappingIntrons[0]
				rIntron := overlappingIntrons[len(overlappingIntrons)-1]

				if lIntron.Start == lGap.Start && rIntron.End == lGap.End {
					// if gap already complies with intron boundaries, don't do anything and skip to next anchor
					mainAnchor = weakAnchor
					mainAnchorIndex--
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
					// we check for padding if part of the main anchor is reaching into intron
					if mainAnchor.Start < rIntron.End && mainAnchor.End > rIntron.End {
						// before:      |remap|   mainAnchor
						//(READ) -------####------*+++++---------->
						//(REF)  ####*-------------+++++---------->
						//            |  rIntron  |
						padding := rIntron.End - mainAnchor.Start
						// after:       |  remap  |
						//(READ) -------####------*+++++---------->
						//(REF)  ####*-------------+++++---------->
						weakAnchor.End += padding
						mainAnchor.Start += padding
					}

					regionReadEnd, err := regionvector.GenomicCoordToReadCoord(0, weakAnchor.End, readMatchResult.MatchedGenome.Regions)
					if err != nil {
						logrus.Errorf("Error while converting genomic coord to read coord in read %s", read.Header)
						logrus.Fatal(err)
					}
					remapOptions := leftRemap(readMatchResult, targetSeqIntronSet, read, weakAnchor, mainAnchor, genomeIndex, mainAnchorRank, regionReadEnd)

					if remapOptions != nil {
						bestOption := chooseBestRemap(remapOptions)
						if len(remapOptions[bestOption].MismatchesRead) < len(readMatchResult.MismatchesRead) {
							// remove weakAnchor and everything left from it
							readMatchResult.MatchedGenome.RemoveRegion(readMatchResult.MatchedGenome.GetFirstRegion().Start, weakAnchor.End)

							// use remap
							for _, region := range remapOptions[bestOption].RemapVector {
								readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(region.Start, region.End)
							}
							// update mm in readMatchResult
							readMatchResult.MismatchesRead = remapOptions[bestOption].MismatchesRead
						}
					}
				} else {
					// Case III
					// (READ) ++++++++++-------+++++++++++
					// (READ) ++++++-------+++++++++++++++ (or)
					// (REF)  ++++++++-------++++++++++++
					// to
					// (READ) ++++++++-------+++++++++++++
					// (READ) ++++++++-------+++++++++++++
					// (REF)  ++++++++-------++++++++++++
					if mainAnchor.Length() > abs(correctedL) && weakAnchor.Length() > abs(correctedR) {
						mainAnchor.Start = mainAnchor.Start + correctedL
						weakAnchor.End = weakAnchor.End + correctedR
					}
				}
				// go to next anchor
				mainAnchor = weakAnchor
				mainAnchorIndex--
			}
		}

		////////////////////////////////////////////////////////
		///// CHECK IF WE NEED TO CORRECT START OR END REGION //
		////////////////////////////////////////////////////////
		correctOverhangs(readMatchResult, targetSeqIntronSet, read, genomeIndex, readMatchResult.MatchedGenome.GetFirstRegion())
		correctOverhangs(readMatchResult, targetSeqIntronSet, read, genomeIndex, readMatchResult.MatchedGenome.GetLastRegion())
	} else { // here we handle reads with no gaps
		////////////////////////////////////////////////////////
		///// IF READ HAS NO,JUNCTIONS CHECK FOR OVERHANGS /////
		////////////////////////////////////////////////////////
		// if there are no gaps we need to remap overhangs

		correctOverhangs(readMatchResult, targetSeqIntronSet, read, genomeIndex, readMatchResult.MatchedGenome.Regions[0])
	}

	// in any case we want to reprocess the remapped read again in order to annotate splicesites and new mm (caused by padding)
	mm := make([]int, 0)
	readMatchResult.SyncRegions()
	for i := 0; i < len(readMatchResult.MatchedRead.Regions); i++ {
		readRegion := readMatchResult.MatchedRead.Regions[i]
		genomeRegion := readMatchResult.MatchedGenome.Regions[i]
		for j := 0; j < readRegion.Length(); j++ {
			posRead := readRegion.Start + j
			posGenome := genomeRegion.Start + j
			if (*read.Sequence)[posRead] != (*genomeIndex.Sequences[readMatchResult.SequenceIndex])[posGenome] {
				mm = append(mm, readRegion.Start+j)
			}
		}
	}

	// TODO: BRANCH
	if len(mm) > len(beforeRemap.MismatchesRead) {
		*readMatchResult = *beforeRemap
	}
	// TODO: CAN BE REMOVED AFTER BRANCHING
	if readMatchResult.MatchedGenome.Length() != len(*read.Sequence) {
		*readMatchResult = *beforeRemap
	}

	// check if read could be mapped/remapped
	if uint8(float64(len(readMatchResult.MismatchesRead))*100/float64(len(*read.Sequence))) > config.MaxMismatchPercentage() {
		readMatchResult.IncompleteMap = true
	} else {
		readMatchResult.IncompleteMap = false
	}
}

func correctOverhangs(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, genomeIndex *index.GenomeIndex, region *regionvector.Region) {
	overlappingIntrons := targetSeqIntronSet.GetIntersectingIntrons(region)
	if len(overlappingIntrons) == 0 {
		// solid map, no remap possible
		return
	}
	weakAnchor := &regionvector.Region{}
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
	// or  (unlikely i think)
	// +++++++++++++++++++++++++++++ (READ)
	// +++++++++++++++----++++++++++(REF) (deletion) -> we do not want to correct here
	lIntron := overlappingIntrons[0]
	rIntron := overlappingIntrons[len(overlappingIntrons)-1]

	// handle deletion
	if lIntron.Start > mapStart && rIntron.End < mapEnd {
		return
	}

	// check if read start needs left remap
	if lIntron.End > mapStart && mapEnd > lIntron.End {
		//         |mapStart
		//         +++++++++ (READ)
		// +++---------+++++(REF)
		//            |lIntron.End

		// adjust anchors
		weakAnchor.Start = region.Start
		weakAnchor.End = lIntron.End
		region.Start = lIntron.End
		mainAnchorRank := lIntron.Rank + 1
		readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(weakAnchor.Start, weakAnchor.End)
		regionReadEnd, err := regionvector.GenomicCoordToReadCoord(0, weakAnchor.End, readMatchResult.MatchedGenome.Regions)
		if err != nil {
			logrus.Errorf("Error while converting genomic coord to read coord in read %s", read.Header)
			logrus.Fatal(err)
		}
		remapOptions := leftRemap(readMatchResult, targetSeqIntronSet, read, weakAnchor, region, genomeIndex, mainAnchorRank, regionReadEnd)

		if remapOptions != nil {
			bestOption := chooseBestRemap(remapOptions)
			if len(remapOptions[bestOption].MismatchesRead) < len(readMatchResult.MismatchesRead) {

				// remove weakAnchor
				readMatchResult.MatchedGenome.RemoveRegion(readMatchResult.MatchedGenome.GetFirstRegion().Start, weakAnchor.End)

				// use remap
				for _, r := range remapOptions[bestOption].RemapVector {
					readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(r.Start, r.End)
				}
				// update mm in readMatchResult
				readMatchResult.MismatchesRead = remapOptions[bestOption].MismatchesRead
			}
		}
	}

	// check if read end needs right remap
	if rIntron.Start < mapEnd && mapStart < rIntron.Start {
		// remap
		// +++++++++ (READ)
		// +++++---------++++++++ (REF)
		// adjust anchor
		weakAnchor.Start = rIntron.Start
		weakAnchor.End = mapEnd
		region.End = rIntron.Start
		mainAnchorRank := rIntron.Rank - 1

		startInRead, err := regionvector.GenomicCoordToReadCoord(readMatchResult.MatchedRead.GetFirstRegion().Start, weakAnchor.Start, readMatchResult.MatchedGenome.Regions)
		if err != nil {
			logrus.Errorf("Error while converting genomic coord to read coord")
			logrus.Fatal(err)
		}
		remapOptions := rightRemap(readMatchResult, targetSeqIntronSet, read, weakAnchor, region, genomeIndex, mainAnchorRank, startInRead)
		if remapOptions != nil {
			bestOption := chooseBestRemap(remapOptions)
			// TODO: branch if mm == remapMM
			if len(remapOptions[bestOption].MismatchesRead) < len(readMatchResult.MismatchesRead) {

				// remove weakAnchor
				readMatchResult.MatchedGenome.RemoveRegion(weakAnchor.Start, readMatchResult.MatchedGenome.GetLastRegion().End)

				// use remap
				for _, r := range remapOptions[bestOption].RemapVector {
					readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(r.Start, r.End)
				}
				// update mm in readMatchResult
				readMatchResult.MismatchesRead = remapOptions[bestOption].MismatchesRead
			}
		}
	}
}

func rightRemap(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, weakAnchor *regionvector.Region, anchorRegion *regionvector.Region, genomeIndex *index.GenomeIndex, anchorRank int, startInRead int) map[int]*RemapOption {
	refSeq := genomeIndex.Sequences[readMatchResult.SequenceIndex]
	if weakAnchor.End > len(*refSeq) {
		return nil
	}
	remapOptions := make(map[int]*RemapOption)

	regionReadEnd := len(*read.Sequence)

	readSequenceToRemap := (*read.Sequence)[startInRead:regionReadEnd]

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
			regions, score, mm := rightRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq, startInRead)
			mmCombined := alignMismatches(readMatchResult.MismatchesRead, mm, remapPos, true, startInRead)
			remapOptions[remapPos] = &RemapOption{
				RemapStart:       remapPos,
				RemapVector:      regions,
				MismatchesRead:   mmCombined,
				Score:            score,
				DistToMainAnchor: 0, // since remap starts within same exon as anchor
			}
		}

		// remap for all remaining exons (in right direction)
		for i := intronRight.Rank; i < len(targetSeqIntronSet.Regions); i++ {
			remapPos := targetSeqIntronSet.Regions[i].End
			regions, score, mm := rightRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq, startInRead)
			mmCombined := alignMismatches(readMatchResult.MismatchesRead, mm, startInRead, true, startInRead)
			remapOptions[remapPos] = &RemapOption{
				RemapStart:       remapPos,
				RemapVector:      regions,
				MismatchesRead:   mmCombined,
				Score:            score,
				DistToMainAnchor: abs(i - anchorRank), // since remap starts after intron
			}
		}
	} else {
		// if there are no introns left in right direction, start remap from anchorregion end
		remapPos := anchorRegion.End
		regions, score, mm := rightRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq, startInRead)
		mmCombined := alignMismatches(readMatchResult.MismatchesRead, mm, remapPos, true, startInRead)
		remapOptions[remapPos] = &RemapOption{
			RemapStart:       remapPos,
			RemapVector:      regions,
			MismatchesRead:   mmCombined,
			Score:            score,
			DistToMainAnchor: abs(targetSeqIntronSet.Regions[len(targetSeqIntronSet.Regions)-1].Rank - anchorRank), // remap is after last intron
		}
	}

	return remapOptions
}

func leftRemap(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, weakAnchor *regionvector.Region, anchorRegion *regionvector.Region, genomeIndex *index.GenomeIndex, anchorRank, regionReadEnd int) map[int]*RemapOption {
	if weakAnchor.Start < 0 {
		// if this is the case, theres not enough gene left to remap the read
		// GENOME [2,122]
		// READ   [30,150] -> 30 bases missing but cant do left remap since only 2 bases in gene are left
		// then weakAnchor == [-28, 2]
		// retrun and dont do anything
		return nil
	}
	remapOptions := make(map[int]*RemapOption)

	regionReadStart := 0 // remap entire read from 0 to end of weakAnchor

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
			mmCombined := alignMismatches(readMatchResult.MismatchesRead, mm, remapPos, false, regionReadEnd)

			remapOptions[remapPos] = &RemapOption{
				RemapStart:       remapPos,
				RemapVector:      regions,
				MismatchesRead:   mmCombined,
				Score:            score,
				DistToMainAnchor: 0, // since remap starts within same exon as anchor
			}
		}

		// remap for remaining exons (in left direction)
		for i := intronLeft.Rank; i >= 0; i-- {
			remapPos := targetSeqIntronSet.Regions[i].Start - 1
			regions, score, mm := leftRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
			mmCombined := alignMismatches(readMatchResult.MismatchesRead, mm, remapPos, false, regionReadEnd)
			remapOptions[remapPos] = &RemapOption{
				RemapStart:       remapPos,
				RemapVector:      regions,
				MismatchesRead:   mmCombined,
				Score:            score,
				DistToMainAnchor: abs(anchorRank-i) + 1,
			}
		}
	} else {
		// if there are no introns left in left direction, start remap from anchorregion start - 1
		remapPos := anchorRegion.Start - 1
		regions, score, mm := leftRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
		mmCombined := alignMismatches(readMatchResult.MismatchesRead, mm, remapPos, false, regionReadEnd)
		remapOptions[remapPos] = &RemapOption{
			RemapStart:       remapPos,
			RemapVector:      regions,
			MismatchesRead:   mmCombined,
			Score:            score,
			DistToMainAnchor: anchorRank, // remap is before first intron
		}
	}
	return remapOptions
}

type RemapOption struct {
	RemapStart       int
	RemapVector      []*regionvector.Region
	MismatchesRead   []int
	Score            int
	DistToMainAnchor int // this allows us to prioritize remaps if they have the same score
}

func (r *RemapOption) UpdateRemap(remap []*regionvector.Region, mm []int, score int) {
	r.RemapVector = append(r.RemapVector, remap...)
	r.MismatchesRead = append(r.MismatchesRead, mm...)
	r.Score += score
}

// chooseBestRemap returns the best possible remap from all remapOptions. If two options with same score are found, it currently
// returns the remap which is located nearest to the anchorRegion
func chooseBestRemap(remapOptions map[int]*RemapOption) int {
	// add best matching region back to result
	bestOption := -1
	bestScore := -1
	distToAnchor := 100000
	for start, remapOption := range remapOptions {
		if remapOption.Score > bestScore {
			bestScore = remapOption.Score
			bestOption = start
			distToAnchor = remapOption.DistToMainAnchor
		} else if remapOption.Score == bestScore {
			if remapOption.DistToMainAnchor < distToAnchor {
				distToAnchor = remapOption.DistToMainAnchor
				bestOption = start
			}
		}
	}
	return bestOption
}

// alignMismatches combines the mm of a readMatchResult with the mm of a remap
func alignMismatches(mmOld []int, remapMMs []int, bestOption int, isRightRemap bool, remapStart int) []int {
	remapMM := make([]int, 0)
	if isRightRemap {
		// update mm
		for _, mm := range mmOld {
			if mm < remapStart {
				remapMM = append(remapMM, mm)
			}
		}
		remapMM = append(remapMM, remapMMs...)
	} else {
		for i := len(mmOld) - 1; i >= 0; i-- {
			if mmOld[i] > remapStart {
				remapMM = append(remapMM, mmOld[i])
			}
		}
		remapMM = append(remapMM, remapMMs...)
	}

	return remapMM
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

func incomplRemap(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, genomeIndex *index.GenomeIndex) {
	// merge blocks s.t theres a 1 - 1 projection from genomic block to read block
	// then the index of main anchor references the adjacent block in read and mapStart and mapEnd are readRegions[mainAnchorRank].Start / .End
	readMatchResult.NormalizeRegions()

	// first get mainAnchor
	mainAnchor, mainAnchorRank, mainAnchorIndex := readMatchResult.GetLargestAnchor(targetSeqIntronSet)
	readMainAnchor := readMatchResult.MatchedRead.Regions[mainAnchorIndex]

	// and get overlapping introns
	overlappingIntrons := targetSeqIntronSet.GetIntersectingIntrons(mainAnchor)

	// check if main anchor is fully contained in intron -> stop remap and
	lPadding := 0
	rPadding := 0
	if len(overlappingIntrons) > 1 {
		lIntron := overlappingIntrons[0]
		rIntron := overlappingIntrons[len(overlappingIntrons)-1]

		if lIntron.Start > mainAnchor.Start && rIntron.End < mainAnchor.End {
			// READ ++++++++++++++++++++++++++++++++++++ -> not fine but shold never happen
			// REF  ++++--------+++++++---------++++++++
			return
		}
		// or
		// READ           +++++++++++ -> this is fine and can happen
		//                |         |
		// REF  ++++--------+++++++---------++++++++
		// we need padding for both start and end of anchor
		rPadding = mainAnchor.End - rIntron.Start
		mainAnchor.End -= rPadding
		lPadding = lIntron.End - mainAnchor.Start
		mainAnchor.Start += lPadding
		readMainAnchor.Start += lPadding
		readMainAnchor.End -= rPadding
	} else if len(overlappingIntrons) == 1 {
		// if len(overlappingIntrons) > 1 {
		// READ +++++++                           I -> padding
		// READ       |                +++++++    II -> padding
		// READ       |   +++++++      |          III -> abort
		// REF ++++++--------------------++++++++
		intron := overlappingIntrons[0]

		if intron.Start < mainAnchor.Start && intron.End > mainAnchor.End {
			// abort III
			// just in case
			readMatchResult.IncompleteMap = true
			return
		} else if mainAnchor.Start < intron.End && mainAnchor.End > intron.End {
			// padding II
			lPadding = intron.End - mainAnchor.Start
			mainAnchor.Start += lPadding
			readMainAnchor.Start += lPadding
		} else if mainAnchor.End > intron.Start && mainAnchor.Start < intron.Start {
			// padding I
			rPadding = mainAnchor.End - intron.Start
			mainAnchor.End -= rPadding
			readMainAnchor.End -= rPadding
		}
	}

	// -> check which parts of read are missing in mapping
	mappingStart := readMainAnchor.Start
	mappingEnd := readMainAnchor.End

	if mappingStart != 0 {
		// leftRemap
		regionToRemap := &regionvector.Region{Start: mainAnchor.Start - mappingStart, End: mainAnchor.Start}
		remapOptions := leftRemap(readMatchResult, targetSeqIntronSet, read, regionToRemap, mainAnchor, genomeIndex, mainAnchorRank, regionToRemap.Length())
		if remapOptions != nil && readMatchResult.MatchedRead.GetFirstRegion().Start != 0 {
			bestOption := chooseBestRemap(remapOptions)

			// remove everything left of main anchor
			readMatchResult.MatchedGenome.RemoveRegion(readMatchResult.MatchedGenome.GetFirstRegion().Start, mainAnchor.Start)
			readMatchResult.MatchedRead.RemoveRegion(0, mappingStart)

			for _, region := range remapOptions[bestOption].RemapVector {
				readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(region.Start, region.End)
			}

			readMatchResult.MatchedRead.AddRegionNonOverlappingPanic(0, mappingStart)

			// update mm in readMatchResult
			readMatchResult.MismatchesRead = remapOptions[bestOption].MismatchesRead
		}

	}

	if mappingEnd != len(*read.Sequence) {
		// right remap
		regionToRemap := &regionvector.Region{Start: mainAnchor.End, End: mainAnchor.End + (len(*read.Sequence) - mappingEnd)}

		remapOptions := rightRemap(readMatchResult, targetSeqIntronSet, read, regionToRemap, mainAnchor, genomeIndex, mainAnchorRank, mappingEnd)

		if remapOptions != nil && readMatchResult.MatchedRead.GetLastRegion().End != len(*read.Sequence) {
			bestOption := chooseBestRemap(remapOptions)

			// remove everything right of main anchor
			readMatchResult.MatchedGenome.RemoveRegion(mainAnchor.End, readMatchResult.MatchedGenome.GetLastRegion().End)
			readMatchResult.MatchedRead.RemoveRegion(mappingEnd, len(*read.Sequence))

			for _, region := range remapOptions[bestOption].RemapVector {
				readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(region.Start, region.End)
			}

			readMatchResult.MatchedRead.AddRegionNonOverlappingPanic(mappingEnd, len(*read.Sequence))

			// update mm in readMatchResult
			readMatchResult.MismatchesRead = remapOptions[bestOption].MismatchesRead
		}
	}

	// update bool to include res in sam
	if uint8(float64(len(readMatchResult.MismatchesRead))*100/float64(len(*read.Sequence))) > config.MaxMismatchPercentage() {
		readMatchResult.IncompleteMap = true
	} else {
		readMatchResult.IncompleteMap = false
	}
}
