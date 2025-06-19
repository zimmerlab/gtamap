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

	// var wgRemap sync.WaitGroup

	for {
		task, ok := secondPassChan.Receive()
		if !ok {
			break
		}

		// wgRemap.Add(1)
		// go func(t *mapperutils.ReadPairMatchResults) {
		// 	defer wgRemap.Done()

		remapReadPair(task, annotation, genomeIndex)

		logrus.Debugf("Secondpass map: %s", task.ReadPair.ReadR1.Header)

		thirdPassChan.Send(&thirdpass.ThirdPassTask{
			ReadPairId: task.ReadPair.ReadR1.Header,
			TargetInfo: task,
		})
		// }(task)

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
	// wgRemap.Wait()
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
	// remap IncompleteMap
	if readMapping.IncompleteMap {
		logrus.Debug("Have to remap readMatchResult due to: IncompleteMap")
		return
	}
	enforceIntrons(readMapping, annotation.Introns[readMapping.SequenceIndex], read, genomeIndex)
}

func enforceIntrons(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, genomeIndex *index.GenomeIndex) {
	gap := regionvector.Region{}

	// TODO: We should always start with the gap at the largest anchor
	if readMatchResult.MatchedGenome.HasGaps() {
		for i := 0; i < len(readMatchResult.MatchedGenome.Regions)-1; i++ {

			gapStart := readMatchResult.MatchedGenome.Regions[i].End
			gapStop := readMatchResult.MatchedGenome.Regions[i+1].Start
			gap.Start = gapStart
			gap.End = gapStop

			// we return a SLICE of overlapping introns because of cases like this:
			// (READ)      +++-----------------+++++++  -> one junction can span across several introns of inferred annotation
			// (REF)  ++++++++--------+++------+++++++
			overlappingIntrons := targetSeqIntronSet.GetIntersectingIntrons(&gap)

			if len(overlappingIntrons) == 0 {
				// Case I (Deletion)
				// (READ) +++++++-------------+++++++
				// (REF) ++++++++++++++++++++++++++++++++
				logrus.WithFields(logrus.Fields{
					"Implied gap in read": gap,
					"Length of Deletion":  gap.Length(),
					"Read Blocks":         readMatchResult.MatchedGenome.Regions,
				}).Debug("Found deletion in read.")
				// Just log / use for report later
				// check if left or right diag can be extended
				continue
			}

			// (READ)      +++-----------------+++++++
			// (REF)  ++++++++--------+++------+++++++
			//                lIntron    rIntron
			// if only one intron is overlapping rIntron == lIntron
			lIntron := overlappingIntrons[0]                         // left most intron -> end of alignment block
			rIntron := overlappingIntrons[len(overlappingIntrons)-1] // right most intron -> start of next alignment block

			lAnchor := readMatchResult.MatchedGenome.Regions[i]
			rAnchor := readMatchResult.MatchedGenome.Regions[i+1]

			if lIntron.Start == gapStart && rIntron.End == gapStop {
				// if gap already complies with intron boundaries,don't do anything
				continue
			}

			if rIntron.Start < rAnchor.Start && rIntron.End > rAnchor.End {
				// Case II right default
				// (READ) +++++++-------------------+--  -> small kmer gets mapped into some random part of the gene
				// (REF)  +++++++----+++---------------
				// -> we need to check after every inferred intron if the small kmer matches inside the exon

				logrus.WithFields(logrus.Fields{
					"Implied gap in read": gap,
					"Read Blocks":         readMatchResult.MatchedGenome.Regions,
					"Overlapping Introns": overlappingIntrons,
				}).Debug("Found inconsistent junction which has unstable (right) anchor mapped inside intron.")

				remapUsedDiagonal(readMatchResult, targetSeqIntronSet, read, rAnchor, false, lAnchor, genomeIndex)
				continue
			} else if lIntron.Start < lAnchor.Start && lIntron.End > lAnchor.End {
				// Case II left default
				// (READ) -------------+----+++++++++--  -> small kmer gets mapped into some random part of the gene
				// (REF)  +++++-------------+++++++++--
				// -> we need to check after every inferred intron if the small kmer matches inside the exon

				logrus.WithFields(logrus.Fields{
					"Implied gap in read": gap,
					"Read Blocks":         readMatchResult.MatchedGenome.Regions,
					"Overlapping Introns": overlappingIntrons,
				}).Debug("Found inconsistent junction which has unstable (left) anchor mapped inside intron.")

				remapUsedDiagonal(readMatchResult, targetSeqIntronSet, read, lAnchor, true, rAnchor, genomeIndex)
				continue
			}

			logrus.WithFields(logrus.Fields{
				"Implied gap in read": gap,
				"Read Blocks":         readMatchResult.MatchedGenome.Regions,
				"Overlapping Introns": overlappingIntrons,
			}).Debug("Found inconsistent read junction not following inferred intron boundaries")

			correctedL := lIntron.Start - gapStart
			correctedR := rIntron.End - gapStop

			// if correctedL != correctedR -> remapDiagonal
			// this would equal CASE II (Special case)
			// (READ) +++++--------------------++++++
			// (REF)  ++++++++-------+++++++++++++++++++++++
			if correctedR != correctedL {
				// we need to determine weak anchor and remap completely
				if lAnchor.Length() < rAnchor.Length() {
					remapUsedDiagonal(readMatchResult, targetSeqIntronSet, read, lAnchor, true, rAnchor, genomeIndex)
				} else {
					remapUsedDiagonal(readMatchResult, targetSeqIntronSet, read, rAnchor, false, lAnchor, genomeIndex)
				}
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
			lAnchor.End = lAnchor.End + correctedL
			rAnchor.Start = rAnchor.Start + correctedR
		}
	} else {
		// only check overhangs
	}
}

// remapUsedDiagonal gets a genomic region as input and remaps it by trying out every combination of introns
func remapUsedDiagonal(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, regionToRemap *regionvector.Region, isLeftRemap bool, anchorRegion *regionvector.Region, genomeIndex *index.GenomeIndex) {
	// remove already used region we want to remap from result
	regionReadStart := regionvector.GenomicCoordToReadCoord(readMatchResult.MatchedRead.GetFirstRegion().Start, regionToRemap.Start, readMatchResult.MatchedGenome.Regions)
	regionReadEnd := regionvector.GenomicCoordToReadCoord(readMatchResult.MatchedRead.GetFirstRegion().Start, regionToRemap.End, readMatchResult.MatchedGenome.Regions)

	readMatchResult.MatchedGenome.RemoveRegion(regionToRemap.Start, regionToRemap.End)
	readMatchResult.MatchedRead.RemoveRegion(regionReadStart, regionReadEnd)

	readSequenceToRemap := (*read.Sequence)[regionReadStart:regionReadEnd]
	refSeq := genomeIndex.Sequences[readMatchResult.SequenceIndex]

	bestRemap := -1
	bestScore := -1
	mmRemap := make(map[int][]int)
	regionsRemap := make(map[int][]*regionvector.Region)
	if isLeftRemap { // only look at exon ends (excluding last exon)

		// get next intron
		intronInfront := targetSeqIntronSet.GetPrevIntron(anchorRegion.Start)
		if intronInfront != nil {
			// check if start of anchorRegion is at intron boundary
			//                                 |anchor
			// (READ) -----++------------------+++++++->
			//  (REF) -----------------+++++++++++++++->
			//           nextIntron    |boundary
			if anchorRegion.Start > intronInfront.End {
				remapPos := anchorRegion.Start - 1
				regions, score, mm := leftRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
				mmRemap[remapPos] = mm
				regionsRemap[remapPos] = regions
				if bestScore < score {
					bestRemap = remapPos
					bestScore = score
				}
			}

			// check for potential overhangs into neighbor intron of anchor seed
			// (READ) ------####------*+++++---------->
			// (REF)  ####-------------+++++---------->
			if anchorRegion.Start < intronInfront.End {
				padding := intronInfront.End - anchorRegion.Start
				// add padding to region to remap
				// before:     |remap|
				//(READ) -------####------*+++++---------->
				//(REF)  ####*-------------+++++---------->
				// after:       |  remap  |
				//(READ) -------####------*+++++---------->
				//(REF)  ####*-------------+++++---------->
				regionReadEnd = regionReadEnd + padding

				readMatchResult.MatchedGenome.RemoveRegion(anchorRegion.Start, anchorRegion.Start+padding)
				readMatchResult.MatchedRead.RemoveRegion(regionReadEnd-padding, regionReadEnd)

				anchorRegion.Start = anchorRegion.Start + padding

				readSequenceToRemap = (*read.Sequence)[regionReadStart:regionReadEnd]
			}

			// remap for remaining exons (in left direction)
			for i := intronInfront.Rank; i >= 0; i-- {
				refStart := targetSeqIntronSet.Regions[i].Start - 1
				regions, score, mm := leftRemapAlignmentBlockFromPos(refStart, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
				mmRemap[refStart] = mm
				regionsRemap[refStart] = regions
				if bestScore < score {
					bestRemap = refStart
					bestScore = score
				}
			}
		} else {
			remapPos := anchorRegion.Start - 1
			regions, score, mm := leftRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
			mmRemap[remapPos] = mm
			regionsRemap[remapPos] = regions
			if bestScore < score {
				bestRemap = remapPos
				bestScore = score
			}
		}

		if read.Header == "2674361" {
			fmt.Println()
		}

		// add best matching region back to result
		for _, region := range regionsRemap[bestRemap] {
			readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(region.Start, region.End)
		}
		readMatchResult.MatchedRead.AddRegionNonOverlappingPanic(regionReadStart, regionReadEnd)

		// now we check if the region we just remapped skips any other region downstream of the map
		for _, region := range readMatchResult.MatchedGenome.Regions {
			// skip region we just remapped
			if region.Start == regionsRemap[bestRemap][0].Start {
				continue
			}
			// skip all regions before anchor region
			if region.Start <= anchorRegion.Start {
				continue
			}
			if region.Start < bestRemap {
				fmt.Println(read.Header)
				fmt.Println("SKIPPED REGION DURING LEFT REMAP -> CASCADING LEFT")

				// as soon as we skip any region while remapping to the left, we have to also remap that region into the same direction
				// rAnchor is the region we remapped
				rAnchor := regionsRemap[bestRemap][0]
				remapUsedDiagonal(readMatchResult, targetSeqIntronSet, read, region, true, rAnchor, genomeIndex)
			}
		}
	} else { // right remap
		// get next intron
		intronBehind := targetSeqIntronSet.GetNextIntron(anchorRegion.Start)
		if intronBehind != nil {
			// first check if end of anchorRegion is at intron boundary
			//        |anchor
			// (READ) +++++++--------------++---------->
			//  (REF) +++++++++++---------------------->
			//                   | boundary
			if anchorRegion.End < intronBehind.Start {
				remapPos := anchorRegion.End
				regions, score, mm := rightRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
				mmRemap[remapPos] = mm
				regionsRemap[remapPos] = regions
				if bestScore < score {
					bestRemap = remapPos
					bestScore = score
				}
			}

			// check for potential overhangs into neighbor intron of anchor seed
			//(READ) ++++++++++*-----######---------->
			//(REF) +++++++++++---------------*######>
			if anchorRegion.End > intronBehind.Start {
				padding := anchorRegion.End - intronBehind.Start

				// add padding to region to remap
				// before:               |remap|
				// READ: ++++++++++*-----#######--------->
				//  REF: ++++++++++---------------*######>
				// after:          |   remap   |
				// READ: ++++++++++*-----#######--------->
				//  REF: ++++++++++---------------*######>
				regionReadStart = regionReadStart - padding
				readMatchResult.MatchedGenome.RemoveRegion(anchorRegion.End-padding, anchorRegion.End)
				readMatchResult.MatchedRead.RemoveRegion(regionReadStart, regionReadStart+padding)

				anchorRegion.End = anchorRegion.End - padding

				readSequenceToRemap = (*read.Sequence)[regionReadStart:regionReadEnd]
			}

			// remap for remaining exons (in right direction)
			for i := intronBehind.Rank; i < len(targetSeqIntronSet.Regions); i++ {
				refStart := targetSeqIntronSet.Regions[i].End
				regions, score, mm := rightRemapAlignmentBlockFromPos(refStart, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
				mmRemap[refStart] = mm
				regionsRemap[refStart] = regions
				if bestScore < score {
					bestRemap = refStart
					bestScore = score
				}
			}
		} else {
			remapPos := anchorRegion.End
			regions, score, mm := rightRemapAlignmentBlockFromPos(remapPos, anchorRegion, targetSeqIntronSet, readSequenceToRemap, refSeq)
			mmRemap[remapPos] = mm
			regionsRemap[remapPos] = regions
			if bestScore < score {
				bestRemap = remapPos
				bestScore = score
			}
		}

		if read.Header == "2674361" {
			fmt.Println()
		}
		// add best matching region back to result
		for _, region := range regionsRemap[bestRemap] {
			readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(region.Start, region.End)
		}
		readMatchResult.MatchedRead.AddRegionNonOverlappingPanic(regionReadStart, regionReadEnd)

		// now we check if the region we just remapped skips any other region downstream of the map
		for _, region := range readMatchResult.MatchedGenome.Regions {
			// skip region we just remapped
			if region.Start == regionsRemap[bestRemap][0].Start {
				continue
			}
			// skip all regions before anchor region
			if region.Start <= anchorRegion.Start {
				continue
			}
			if region.Start < bestRemap {
				fmt.Println(read.Header)
				fmt.Println("SKIPPED REGION DURING RIGHT REMAP -> CASCADING RIGHT")
				// as soon as we skip any region while remapping to the right, we have to also remap that region into the same direction
				// lAnchor is the last mapped region to the right of our remap
				lAnchor := regionsRemap[bestRemap][len(regionsRemap[bestRemap])-1]
				remapUsedDiagonal(readMatchResult, targetSeqIntronSet, read, region, false, lAnchor, genomeIndex)
			}
		}
	}
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
	score := 0
	mm := make([]int, 0)
	explainedBases := 0

	for j := len(readSequenceToRemap) - 1; j >= 0; j-- {

		// jump to next exon in this case
		//                      | jump here
		// (READ) ----------####****-------+++++++++++---->
		//  (REF) -####----------------****+++++++++++---->
		//
		if nextIntronFromAnchor != nil && refPos == nextIntronFromAnchor.End {
			remap = append(remap, &regionvector.Region{Start: refPos, End: remapStart + 1})

			if readSequenceToRemap[j] == (*refSeq)[refPos] {
				score++
			} else {
				mm = append(mm, refPos)
			}

			refPos = nextIntronFromAnchor.Start - 1
			remapStart = nextIntronFromAnchor.Start - 1
			continue
		}

		if readSequenceToRemap[j] == (*refSeq)[refPos] {
			score++
		} else {
			mm = append(mm, refPos)
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
	mainAnchor, mainAnchorPos := readMatchResult.MatchedGenome.GetLargestAnchor()

	// Check if there are gaps to thr right direction of the anchor
}
