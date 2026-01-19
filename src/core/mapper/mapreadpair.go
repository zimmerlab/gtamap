package mapper

import (
	"github.com/KleinSamuel/gtamap/src/core/mapper/confidentmappingpass"
	"github.com/KleinSamuel/gtamap/src/core/mapper/events"
	"github.com/KleinSamuel/gtamap/src/core/mapper/secondpass"

	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/core/timer"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/sirupsen/logrus"
)

// var considered int = 0
// var keepAfterFilter int = 0
func MapReadPair(
	readPair *fastq.ReadPair,
	genomeIndex *index.GenomeIndex,
	secondpassChan *secondpass.SecondPassChannel,
	confidentMatchesChan *confidentmappingpass.ConfidentPassChan,
	timerChannel chan<- *timer.Timer,
	progressChan chan<- events.Event,
	progressStats *ProgressStats,
	// mu *sync.Mutex,
) {
	keepFw := GlobalFilter(readPair.ReadR1.Sequence, genomeIndex)

	keepRw := true

	if readPair.ReadR2 != nil {
		keepRw = GlobalFilter(readPair.ReadR2.Sequence, genomeIndex)
	}

	// keepFw := BinnedFilter(readPair.ReadR1.Sequence, genomeIndex)
	// keepRw := BinnedFilter(readPair.ReadR2.Sequence, genomeIndex)

	// logrus.WithFields(logrus.Fields{
	// 	"keepFw": keepFw,
	// 	"keepRv": keepRw,
	// }).Debug("Filter results")

	// TODO: allow single-end mapping?
	if !keepFw || !keepRw {
		return
	}

	// mu.Lock()
	progressStats.ReadsAfterFiltering++
	// mu.Unlock()
	if progressStats.ReadsAfterFiltering >= 100 {
		progressChan <- events.Event{
			Type: events.EventTypeReadsAfterFiltering,
			Data: progressStats.ReadsAfterFiltering,
		}
		progressStats.ReadsAfterFiltering = 0
	}

	resultFw, isMappableFw := MapRead(readPair.ReadR1, genomeIndex, false)

	// TODO: allow mapping of only one read in pair?
	if !isMappableFw || len(resultFw) == 0 {
		return
	}

	resultRv := []*mapperutils.ReadMatchResult{}
	isMappableRv := false

	if readPair.ReadR2 != nil {
		resultRv, isMappableRv = MapRead(readPair.ReadR2, genomeIndex, false)

		// TODO: REMOVE DEBUG
		// if isMappableFw {
		// 	debugout.GenerateAlignmentView(genomeIndex, resultFw[0], readPair.ReadR1)
		// }
		// if isMappableRv {
		// 	debugout.GenerateAlignmentView(genomeIndex, resultRv[0], readPair.ReadR2)
		// }

		if !isMappableRv || len(resultRv) == 0 {
			// logrus.WithFields(logrus.Fields{
			// 	"isMappableFw":  isMappableFw,
			// 	"isMappableRv":  isMappableRv,
			// 	"num resultsFw": len(resultFw),
			// 	"num resultsRv": len(resultRv),
			// }).Debug("readpair not mappable")
			return
		}
	}

	// mu.Lock()
	progressStats.ReadsMapped++
	progressStats.NumMappingLocations += uint64(len(resultFw) + len(resultRv))
	// mu.Unlock()

	if progressStats.ReadsMapped >= 100 {
		progressChan <- events.Event{
			Type: events.EventTypeReadsMapped,
			Data: progressStats.ReadsMapped,
		}
		progressStats.ReadsMapped = 0
	}
	if progressStats.NumMappingLocations >= 1000 {
		progressChan <- events.Event{
			Type: events.EventTypeNumMappingLocations,
			Data: progressStats.NumMappingLocations,
		}
		progressStats.NumMappingLocations = 0
	}

	// postprocess every potential match
	// commented out for now
	// IMPORTANT: Currently leads to bug because read rv is not ajdusted when genome rv is adjsuted.
	// beofre postprocessReadMatch, every region in MatchedGenome and MatchedRead match 1:1 but
	// with the current implementation, this 1:1 mapping is destroyed.
	// for _, resFw := range resultFw {
	// 	postprocessReadMatch(genomeIndex, readPair.ReadR1, resFw)
	// }
	// for _, resRv := range resultRv {
	// 	postprocessReadMatch(genomeIndex, readPair.ReadR2, resRv)
	// }

	secondpassChan.Send(&mapperutils.ReadPairMatchResults{
		ReadPair: readPair,
		Fw:       resultFw,
		Rv:       resultRv,
	})

	// check if confident map
	// hasConfdentConfiguration allows to discover conf maps even if len(fw) > 1 && len(rv) > 1
	possibleConfMap, confMap := hasConfdentConfiguration(resultFw, resultRv, readPair)
	if possibleConfMap {
		confidentMatchesChan.Send(confMap)

		// mu.Lock()
		progressStats.NumConfidentMappings++
		// mu.Unlock()
		if progressStats.NumConfidentMappings >= 10 {
			progressChan <- events.Event{
				Type: events.EventTypeNumConfidentMappings,
				Data: progressStats.NumConfidentMappings,
			}
			progressStats.NumConfidentMappings = 0
		}
	}
	// alternatively use more strict way of determining conf map
	// if isStrictConfidentMap(resultFw, resultRv) {
	// 	confidentMatchesChan.Send(&confidentmappingpass.ConfidentTask{
	// 		ReadPair: readPair,
	// 		ResultFw: resultFw[0], // there should only exist fw[0] and rv[0] in a confident match
	// 		ResultRv: resultRv[0],
	// 	})
	// }

	// readPairMapping := &mapperutils.ReadPairMatchResults{
	// 	ReadPair: readPair,
	// 	Fw:       resultFw,
	// 	Rv:       resultRv,
	// }
	// here it is okay to also pass the pointers of resultFw and resultRv since paralogMappingChan is readOnly
	// paralogMappingChan <- readPairMapping
}

func isStrictConfidentMap(resultFw []*mapperutils.ReadMatchResult, resultRv []*mapperutils.ReadMatchResult) bool {
	return len(resultFw) == 1 && len(resultRv) == 1 && len(resultFw[0].MismatchesRead)+len(resultRv[0].MismatchesRead) < 6 && (resultRv[0].SequenceIndex-1 == resultFw[0].SequenceIndex || resultRv[0].SequenceIndex == resultFw[0].SequenceIndex-1)
}

func hasConfdentConfiguration(resultFw []*mapperutils.ReadMatchResult, resultRv []*mapperutils.ReadMatchResult, rp *fastq.ReadPair) (bool, *confidentmappingpass.ConfidentTask) {
	fwMapPerSeqIndex, rvMapPerSeqIndex, mappedIds := mapperutils.AssignReadMatchResults(resultFw, resultRv)

	for targetId := range mappedIds {
		fwMapsOfTargetId := fwMapPerSeqIndex[targetId]
		rvMapsOfTargetId := rvMapPerSeqIndex[targetId]
		validCombination := mapperutils.GetBestPossibleMappingCombination(fwMapsOfTargetId, rvMapsOfTargetId)
		if validCombination != nil {
			return true, &confidentmappingpass.ConfidentTask{
				ReadPair: rp,
				ResultFw: validCombination.Fw,
				ResultRv: validCombination.Rv,
			}
		}
	}
	return false, nil
}

func postprocessReadMatch(genomeIndex *index.GenomeIndex, read *fastq.Read, result *mapperutils.ReadMatchResult) {
	applyLeftNormalization(genomeIndex, read, result)
}

// applyLeftNormalization shifts the mapping in genome to the left if there are gaps in the genome
// that have ambiguous bases (N) before and after the gap. If this is the case, the mapping can not
// be determined, and it is up to the mapper to decide where to place the read.
// This function shifts the gap to the leftmost position.
func applyLeftNormalization(genomeIndex *index.GenomeIndex, read *fastq.Read, result *mapperutils.ReadMatchResult) {
	if len(result.MatchedGenome.Regions) < 2 {
		logrus.WithFields(logrus.Fields{
			"read":      read.Header,
			"isForward": genomeIndex.IsSequenceForward(result.SequenceIndex),
		}).Debug("No left normalization required")
		return
	}

	logrus.WithFields(logrus.Fields{
		"read":      read.Header,
		"isForward": genomeIndex.IsSequenceForward(result.SequenceIndex),
	}).Debug("apply left normalization")

	logrus.WithFields(logrus.Fields{
		"genome": result.MatchedGenome,
	}).Debug("match before left normalization")

	if genomeIndex.IsSequenceForward(result.SequenceIndex) {
		determineLeftNormalizationShiftFw(genomeIndex, read, result)
	} else {
		determineLeftNormalizationShiftRv(genomeIndex, read, result)
	}

	logrus.WithFields(logrus.Fields{
		"genome": result.MatchedGenome,
	}).Debug("match after left normalization")
}

// determineLeftNormalizationShiftFw determines the left normalization shift for forward reads.
// It finds gaps in the genome and shifts the mapping in genome to the left based on the
// value of shift.
// The shift is the number of equal bases between the suffix of the read sequence before the gap
// and the suffix of the genome sequence within the gap.
// The maximum value of the shift is the size if the gap.
func determineLeftNormalizationShiftFw(
	genomeIndex *index.GenomeIndex,
	read *fastq.Read,
	result *mapperutils.ReadMatchResult,
) {
	logrus.Debug("determine left normalization for forward seq")

	regionIndexBeforeGap := result.MatchedGenome.GetGapIndexAfterPos(0)
	gapRank := 0 // stores the gap rank (first gap = 0, second gap = 1 ...)

	// find gaps in genome (bases in reference that are not present in read)
	for regionIndexBeforeGap > -1 {

		gapGenome, _ := result.MatchedGenome.GetGapAfterRegionIndex(regionIndexBeforeGap)

		// skip gaps which have known splice sites
		if result.SpliceSitesInfo[gapRank] > 0 {
			gapRank++
			regionIndexBeforeGap = result.MatchedGenome.GetGapIndexAfterPos(gapGenome.End + 1)
			continue
		}

		// only handle gaps in genome that have no mutual gap in read
		// (only introns or deletions)
		gapRead, ok := result.MatchedRead.GetGapAfterRegionIndex(regionIndexBeforeGap)
		if !ok {
			logrus.WithFields(logrus.Fields{
				"gapRead":   gapRead,
				"gapGenome": gapGenome,
			}).Debug("skip gap in genome that has gap in read")
			regionIndexBeforeGap = result.MatchedGenome.GetGapIndexAfterPos(gapGenome.End + 1)
			continue
		}

		// the rightmost position of the subsequence to be tested in the read
		// (last position in the read before the gap in forward direction)
		rRead, sizeErr := result.MatchedGenome.GetSizeLeftIncluding(regionIndexBeforeGap)
		if sizeErr != nil {
			logrus.Fatal("Error getting size left including", sizeErr)
		}

		// the rightmost position of the subsequence to be tested in the genome
		// (end of the gap in the genome)
		rGenome := gapGenome.End

		// the amount of bases that the gap must be shifted to be left normalized
		shift := 0
		// the maximum amount that the gap can be shifted is until the beginning of the read
		maxShift := rRead

		for i := 0; i < maxShift; i++ {

			charRead := (*read.Sequence)[rRead-1-i]

			charGenome := (*genomeIndex.Sequences[result.SequenceIndex])[rGenome-1-i]

			if charRead != charGenome {
				break
			}

			shift++
		}

		logrus.WithFields(logrus.Fields{
			"shift": shift,
		}).Debug("determined shift")

		// shift the mapping in genome to left based on value of shift
		result.MatchedGenome.Regions[regionIndexBeforeGap].End -= shift
		result.MatchedGenome.Regions[regionIndexBeforeGap+1].Start -= shift

		regionIndexBeforeGap = result.MatchedGenome.GetGapIndexAfterPos(gapGenome.End + 1)
	}
}

// determineLeftNormalizationShiftRv determines the left normalization shift for reverse reads.
// The procedure is equivalent to the one for forward reads, but the direction of the
// comparison is reversed.
// The shift is now determined by the number of equal bases between the prefix of the read sequence
// after the gap and the prefix of the genome sequence within the gap.
func determineLeftNormalizationShiftRv(
	genomeIndex *index.GenomeIndex,
	read *fastq.Read,
	result *mapperutils.ReadMatchResult,
) {
	logrus.Debug("determine left normalization for reverse seq")

	regionIndexBeforeGap := result.MatchedGenome.GetGapIndexAfterPos(0)
	gapRank := 0 // stores the gap rank (first gap = 0, second gap = 1 ...)

	// find gaps in genome (bases in reference that are not present in read)
	for regionIndexBeforeGap > -1 {

		gapGenome, _ := result.MatchedGenome.GetGapAfterRegionIndex(regionIndexBeforeGap)

		// skip gaps which have known splice sites
		if result.SpliceSitesInfo[gapRank] > 0 {
			gapRank++
			regionIndexBeforeGap = result.MatchedGenome.GetGapIndexAfterPos(gapGenome.End + 1)
			continue
		}

		// only handle gaps in genome that have no mutual gap in read
		// (only introns or deletions)
		gapRead, ok := result.MatchedRead.GetGapAfterRegionIndex(regionIndexBeforeGap)
		if !ok {
			logrus.WithFields(logrus.Fields{
				"gapRead":   gapRead,
				"gapGenome": gapGenome,
			}).Debug("skip gap in genome that has gap in read")
			regionIndexBeforeGap = result.MatchedGenome.GetGapIndexAfterPos(gapGenome.End + 1)
			continue
		}

		// the rightmost position of the subsequence to be tested in the read
		// (last position in the read before the gap in forward direction)
		// actually this is the first position in the read after the gap (because rv)
		rRead, sizeErr := result.MatchedGenome.GetSizeLeftIncluding(regionIndexBeforeGap)
		if sizeErr != nil {
			logrus.Fatal("Error getting size left including", sizeErr)
		}

		// the rightmost position of the subsequence to be tested in the genome
		// (end of the gap in the genome in forward direction)
		// actually this is the first position in the gap (because rv)
		rGenome := gapGenome.Start

		// the amount of bases that the gap must be shifted to be left normalized
		shift := 0
		// the maximum amount that the gap can be shifted is until the end of the read
		maxShift := result.MatchedGenome.Length() - rRead

		for i := 0; i < maxShift; i++ {

			charRead := (*read.Sequence)[rRead+i]

			charGenome := (*genomeIndex.Sequences[result.SequenceIndex])[rGenome+i]

			if charRead != charGenome {
				break
			}

			shift++
		}

		logrus.WithFields(logrus.Fields{
			"shift": shift,
		}).Debug("determined shift")

		// shift the mapping in genome to right based on value of shift
		// result.MatchedGenome.Regions[i].End += shift
		// result.MatchedGenome.Regions[i+1].Start += shift

		result.MatchedGenome.Regions[regionIndexBeforeGap].End += shift
		result.MatchedGenome.Regions[regionIndexBeforeGap+1].Start += shift

		regionIndexBeforeGap = result.MatchedGenome.GetGapIndexAfterPos(gapGenome.End + 1)
	}
}
