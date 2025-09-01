package mapper

import (
	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/index"
)

// Filter the read
// TODO: finalize the filter step as this is the most crucial step in terms of runtime
// TODO: this has to be adjusted in order to allow for multiple targets in index
// it should then return a list of potential sequences (n,n+1 grouped)
func GlobalFilter(readSequence *[]byte, genomeIndex *index.GenomeIndex) bool {
	// indexHasBlacklist := len(genomeIndex.RepeatRegions) > 0
	numMatching := 0

	numMatchingFw := 0
	numMatchingRv := 0

	for i := 0; i <= len(*readSequence)-(int(config.KmerLength())); i += int(config.KmerLength()) {

		kmer := (*readSequence)[i : i+int(config.KmerLength())]

		// matches := genomeIndex.KeywordTree.FindKeyword(&kmer, i)
		matches := genomeIndex.GetKeywordFromMap(*(*[10]byte)(kmer))

		okFw := false
		okRv := false
		for _, m := range matches {
			//if indexHasBlacklist {
			//	repeatIndices := genomeIndex.Blacklist[int(m.SequenceIndex)].Including(float64(m.Position))
			//	if len(repeatIndices) > 1 {
			//		break
			//	}
			//}
			if m.SequenceIndex == 0 {
				okFw = true
			} else {
				okRv = true
			}
			if okRv && okFw {
				break // early break
			}
		}

		if okFw {
			numMatchingFw++
		}
		if okRv {
			numMatchingRv++
		}

		if matches != nil {
			numMatching++
		}
	}

	// TODO: Needs to be dynamic based on read length
	// INFO: DNA RNA MODE
	if config.IsOriginRNA {
		return numMatchingFw >= 6 || numMatchingRv >= 6
	} else {
		return numMatchingFw >= 7 || numMatchingRv >= 7
	}
}

func BinnedFilter(readSequence *[]byte, genomeIndex *index.GenomeIndex) bool {
	if config.IsOriginRNA {
		// TODO: Needs to be dynamic based on read length
		// return FilterWithBins(readSequence, genomeIndex, 80, 3)
		return FilterWithOverlappingBins(readSequence, genomeIndex, 40, 3)
	} else {
		// TODO: Needs to be dynamic based on read length
		// return FilterWithBins(readSequence, genomeIndex, 150, 7)
		return FilterWithOverlappingBins(readSequence, genomeIndex, 30, 2)
	}
}

func FilterWithBins(readSequence *[]byte, genomeIndex *index.GenomeIndex, binSize uint32, kmerMatchThreshold int) bool {
	kmerLen := int(config.KmerLength())
	readLen := len(*readSequence)

	if readLen < kmerLen || kmerLen != 10 {
		return false
	}

	binCounts := make(map[uint32]int, 16)
	readBytes := *readSequence

	for i := 0; i <= readLen-kmerLen; i += kmerLen {
		kmerKey := [10]byte{
			readBytes[i], readBytes[i+1], readBytes[i+2], readBytes[i+3], readBytes[i+4],
			readBytes[i+5], readBytes[i+6], readBytes[i+7], readBytes[i+8], readBytes[i+9],
		}

		matches := genomeIndex.GetKeywordFromMap(kmerKey)
		for _, m := range matches {

			binID := m.Position / binSize
			binCounts[binID]++
			if binCounts[binID] >= kmerMatchThreshold {
				return true
			}
		}
	}

	return false
}

func FilterWithOverlappingBins(readSequence *[]byte, genomeIndex *index.GenomeIndex, binSize uint32, kmerMatchThreshold int) bool {
	kmerLen := int(config.KmerLength())
	readLen := len(*readSequence)
	if readLen < kmerLen || kmerLen != 10 {
		return false
	}

	binCounts := make(map[uint32]int, 32)
	readBytes := *readSequence

	for i := 0; i <= readLen-kmerLen; i += kmerLen {
		kmerKey := [10]byte{
			readBytes[i], readBytes[i+1], readBytes[i+2], readBytes[i+3], readBytes[i+4],
			readBytes[i+5], readBytes[i+6], readBytes[i+7], readBytes[i+8], readBytes[i+9],
		}

		matches := genomeIndex.GetKeywordFromMap(kmerKey)
		for _, m := range matches {
			binID1 := m.Position / binSize
			binID2 := (m.Position + binSize/2) / binSize

			binCounts[binID1]++
			if binCounts[binID1] >= kmerMatchThreshold {
				return true
			}

			// second bin set with offset marker to distinguish
			offsetBinID := binID2 | 0x80000000 // high bit to mark offset bins
			binCounts[offsetBinID]++
			if binCounts[offsetBinID] >= kmerMatchThreshold {
				return true
			}
		}
	}
	return false
}
