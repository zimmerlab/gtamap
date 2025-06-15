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

		remapReadPair(task, annotation)

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

func remapReadPair(readPairMapping *mapperutils.ReadPairMatchResults, annotationMap map[int]*mapperutils.TargetAnnotation) {
	for _, mapping := range readPairMapping.Fw {
		mainSeqId := mapping.SequenceIndex / 2
		if readPairMapping.ReadPair.ReadR1.Header == "436071" {
			fmt.Println("s")
			x := annotationMap[mapping.SequenceIndex/2].Introns[mapping.SequenceIndex].GetIntersectingIntrons(mapping.MatchedGenome.GetFirstGap())
			fmt.Println(x)
		}

		if needsRemap(mapping, annotationMap[mainSeqId]) {
			remapRead(mapping, annotationMap[mainSeqId], readPairMapping.ReadPair.ReadR1)
		}
	}
	for _, mapping := range readPairMapping.Rv {
		mainSeqId := mapping.SequenceIndex / 2
		if needsRemap(mapping, annotationMap[mainSeqId]) {
			remapRead(mapping, annotationMap[mainSeqId], readPairMapping.ReadPair.ReadR2)
		}
	}
}

func remapRead(readMatchResult *mapperutils.ReadMatchResult, targetAnnotation *mapperutils.TargetAnnotation, read *fastq.Read) {
	enforceIntrons(readMatchResult, targetAnnotation.Introns[readMatchResult.SequenceIndex], read)
}

func enforceIntrons(readMatchResult *mapperutils.ReadMatchResult, targetSeqIntronSet *regionvector.RegionSet, read *fastq.Read) {
	// here we megre intervals in the genomic regions for easier handling
	readMatchResult.MatchedGenome.MergeAlignmentBlocks()

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

			if rIntron.Start < alignedGenomeBlocks[i+1].Start && rIntron.End > alignedGenomeBlocks[i+1].End || lIntron.Start < alignedGenomeBlocks[i].Start && lIntron.End > alignedGenomeBlocks[i].End {
				// Case III
				// read +++++++-------------------+--  -> small kmer gets mapped into some random part of the gene
				// ++++++++++++----+++--------------- (REF)
				// -> we need to check after every inferred intron if the small kmer matches inside the exon

				logrus.WithFields(logrus.Fields{
					"Implied gap in read": gap,
					"Read Blocks":         alignedGenomeBlocks,
					"Overlapping Introns": overlappingIntrons,
				}).Debug("Found inconsistent junction which has instable anchors mapped inside intron")

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
