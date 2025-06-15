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

	for {
		task, ok := secondPassChan.Receive()
		if !ok {
			break
		}

		remapReadPair(task, annotation, genomeIndex)

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

		logrus.Debugf("Secondpass map: %s", task.ReadPair.ReadR1.Header)

		thirdPassChan.Send(&thirdpass.ThirdPassTask{
			ReadPairId: task.ReadPair.ReadR1.Header,
			TargetInfo: task,
		})
	}

	logrus.Info("Done with second pass")
}

func needsRemap(readMapping *mapperutils.ReadMatchResult, annotation *mapperutils.TargetAnnotation) bool {
	// remap IncompleteMap
	if readMapping.IncompleteMap {
		logrus.Debug("Have to remap readMatchResult due to: IncompleteMap")
		return true
	}
	if annotation.Introns[readMapping.SequenceIndex].IntersectsIntrons(readMapping.MatchedGenome.GetAlignmentBlocks()) {
		logrus.Debug("Have to remap readMatchResult due to: Intron boundaries not respected")
		return true
	} else if readMapping.MatchedGenome.HasGaps() {
		// if the read doesn't intersect with introns but has at least one gap
		return true
	}
	return false
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
}

func remapRead(readMapping *mapperutils.ReadMatchResult, annotation *mapperutils.TargetAnnotation, read *fastq.Read, genomeIndex *index.GenomeIndex) {
	// here we megre intervals in the genomic regions for easier handling
	readMapping.MatchedGenome.MergeAlignmentBlocks()

	// remap IncompleteMap
	if readMapping.IncompleteMap {
		logrus.Debug("Have to remap readMatchResult due to: IncompleteMap")
	}
	// remap partially overlapping introns and annotate deletions
	if annotation.Introns[readMapping.SequenceIndex].IntersectsIntrons(readMapping.MatchedGenome.GetAlignmentBlocks()) {
		logrus.Debug("Have to remap readMatchResult due to: Intron boundaries not respected")
		enforceIntrons(readMapping, annotation.Introns[readMapping.SequenceIndex], read, genomeIndex)
	} else if readMapping.MatchedGenome.HasGaps() {
		// if the read doesn't intersect with introns but has at least one gap
	}
}

func enforceIntrons(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, genomeIndex *index.GenomeIndex) {
	alignedGenomeBlocks := readMatchResult.MatchedGenome.Regions

	var gap *regionvector.Region = &regionvector.Region{}

	// here we iterate over all junctions of a read
	// and handle cases I-VI accordingly
	if readMatchResult.MatchedGenome.HasGaps() {
		for i := 0; i < len(alignedGenomeBlocks)-1; i++ {

			gapStart := alignedGenomeBlocks[i].End
			gapStop := alignedGenomeBlocks[i+1].Start
			gap.Start = gapStart
			gap.End = gapStop

			// we return a SLICE of overlapping introns because of cases like this:
			// (READ)      +++-----------------+++++++  -> one junction can span across several introns of inferred annotation
			// (REF)  ++++++++--------+++------+++++++
			overlappingIntrons := targetSeqIntronSet.GetIntersectingIntrons(gap)

			if len(overlappingIntrons) == 0 {
				// Case V (Deletion)
				// (READ) +++++++-------------+++++++
				// (REF) ++++++++++++++++++++++++++++++++
				logrus.WithFields(logrus.Fields{
					"Implied gap in read": gap,
					"Length of Deletion":  gap.Length(),
					"Read Blocks":         alignedGenomeBlocks,
				}).Debug("Found deletion in read.")
				// Just log / use for report later
				continue

			}

			//                lIntron    rIntron
			// (REF)  ++++++++--------+++------+++++++
			// (READ)      +++-----------------+++++++
			lIntron := overlappingIntrons[0]                         // left most intron -> end of alignment block
			rIntron := overlappingIntrons[len(overlappingIntrons)-1] // right most intron -> start of next alignment block

			if rIntron.Start < alignedGenomeBlocks[i+1].Start && rIntron.End > alignedGenomeBlocks[i+1].End {
				// Case III right
				// read +++++++-------------------+--  -> small kmer gets mapped into some random part of the gene
				// ++++++++++++----+++--------------- (REF)
				// -> we need to check after every inferred intron if the small kmer matches inside the exon

				logrus.WithFields(logrus.Fields{
					"Implied gap in read": gap,
					"Read Blocks":         alignedGenomeBlocks,
					"Overlapping Introns": overlappingIntrons,
				}).Debug("Found inconsistent junction which has unstable (right) anchor mapped inside intron")

				//remapUsedDiagonal(readMatchResult, targetSeqIntronSet, read, alignedGenomeBlocks[i+1], false, alignedGenomeBlocks[i], genomeIndex, lIntron.End)
				continue
			} else if lIntron.Start < alignedGenomeBlocks[i].Start && lIntron.End > alignedGenomeBlocks[i].End {
				// Case III left
				// read ------+-----------+++++++++--  -> small kmer gets mapped into some random part of the gene
				// -----------------------+++++++++-- (REF)
				// -> we need to check after every inferred intron if the small kmer matches inside the exon

				logrus.WithFields(logrus.Fields{
					"Implied gap in read": gap,
					"Read Blocks":         alignedGenomeBlocks,
					"Overlapping Introns": overlappingIntrons,
				}).Debug("Found inconsistent junction which has unstable (left) anchor mapped inside intron")

				//remapUsedDiagonal(readMatchResult, targetSeqIntronSet, read, alignedGenomeBlocks[i], true, alignedGenomeBlocks[i+1], genomeIndex, rIntron.Start)
				continue
			}

			logrus.WithFields(logrus.Fields{
				"Implied gap in read": gap,
				"Read Blocks":         alignedGenomeBlocks,
				"Overlapping Introns": overlappingIntrons,
			}).Debug("Found inconsistent read junction not following inferred intron boundaries")

			// Case VI
			correctedL := lIntron.Start - gapStart
			correctedR := rIntron.End - gapStop

			alignedGenomeBlocks[i].End = alignedGenomeBlocks[i].End + correctedL
			alignedGenomeBlocks[i+1].Start = alignedGenomeBlocks[i+1].Start + correctedR
		}
	}
}

// remapUsedDiagonal gets a genomic region as input and remaps it by trying out every combination of introns
func remapUsedDiagonal(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read, regionToRemap *regionvector.Region, isLeftRemap bool, anchorRegion *regionvector.Region, genomeIndex *index.GenomeIndex, intronBoundary int) {
	// remove already used region we want to remap from result
	regionReadStart := regionvector.GenomicCoordToReadCoord(readMatchResult.MatchedRead.GetFirstRegion().Start, regionToRemap.Start, readMatchResult.MatchedGenome.Regions)
	regionReadEnd := regionvector.GenomicCoordToReadCoord(readMatchResult.MatchedRead.GetFirstRegion().Start, regionToRemap.End, readMatchResult.MatchedGenome.Regions)

	readMatchResult.MatchedGenome.RemoveRegion(regionToRemap.Start, regionToRemap.End)
	readMatchResult.MatchedRead.RemoveRegion(regionReadStart, regionReadEnd)

	readSequenceToRemap := (*read.Sequence)[regionReadStart:regionReadEnd]
	refSeq := genomeIndex.Sequences[readMatchResult.SequenceIndex]

	if isLeftRemap {
		// only look at exon ends (excluding last exon)
		nextIntronFromAnchor := targetSeqIntronSet.GetPrevIntron(anchorRegion.End)
		scores := make(map[int]int)
		bestRemap := -1

		// check if start of anchorRegion is at intron boundary
		//                                 |anchor
		// (READ) -----++------------------+++++++->
		//  (REF) -----------------+++++++++++++++->
		//                         | boundary
		if anchorRegion.Start > intronBoundary {
			score := 0
			refStart := anchorRegion.Start - 1
			refPos := refStart
			for j := len(readSequenceToRemap) - 1; j >= 0; j-- {

				// jump to next exon in this case
				//                                  | anchor
				// (READ) -----------####****-------+++++++---------->
				//  (REF) ####------------------****+++++++---------->
				//           | jump to here
				if refPos == nextIntronFromAnchor.End {
					refPos = nextIntronFromAnchor.Start - 1
				}

				if readSequenceToRemap[j] == (*refSeq)[refPos] {
					score++
				}
				refPos++
			}
			if score >= bestRemap {
				bestRemap = refStart
			}
			scores[refStart] = score
		}

		// for all remaining introns with rank <= nextIntronFromAnchor, check matches manually
		for i := nextIntronFromAnchor.Rank; i >= 0; i-- {
			score := 0
			refPos := targetSeqIntronSet.Regions[i].Start - 1
			for j := len(readSequenceToRemap) - 1; j >= 0; j-- {
				if readSequenceToRemap[j] == (*refSeq)[refPos] {
					score++
				}
				refPos--
			}
			if score >= bestRemap {
				bestRemap = i
			}
			scores[i] = score
		}
	} else {
		scores := make(map[int]int)
		bestRemap := -1

		// get next intron
		nextIntronFromAnchor := targetSeqIntronSet.GetNextIntron(anchorRegion.Start)

		// first check if end of anchorRegion is at intron boundary
		//        |anchor
		// (READ) +++++++--------------++---------->
		//  (REF) +++++++++++---------------------->
		//                  | boundary
		if anchorRegion.End < intronBoundary {
			score := 0
			refStart := anchorRegion.End
			refPos := refStart
			for j := 0; j < len(readSequenceToRemap); j++ {

				// jump to next exon in this case
				//        |anchor                  | here we jump
				// (READ) +++++++--------------****####-------------->
				//  (REF) +++++++****--------------------------####++>
				//                  | boundary
				if refPos == nextIntronFromAnchor.Start {
					refPos = nextIntronFromAnchor.End
				}

				if readSequenceToRemap[j] == (*refSeq)[refPos] {
					score++
				}
				refPos++
			}
			if score >= bestRemap {
				bestRemap = refStart
			}
			scores[refStart] = score
		}

		for i := nextIntronFromAnchor.Rank; i < len(targetSeqIntronSet.Regions); i++ {
			score := 0
			refStart := targetSeqIntronSet.Regions[i].Start - 1
			refPos := refStart
			for j := 0; j < len(readSequenceToRemap); j++ {
				if readSequenceToRemap[j] == (*refSeq)[refPos] {
					score++
				}
				refPos++
			}
			if score >= bestRemap {
				bestRemap = refStart
			}
			scores[refStart] = score
		}

		// add best matching region back to result
		readMatchResult.MatchedGenome.AddRegionNonOverlappingPanic(bestRemap, bestRemap+len(readSequenceToRemap))
		readMatchResult.MatchedRead.AddRegionNonOverlappingPanic(regionReadStart, regionReadEnd)
	}
}
