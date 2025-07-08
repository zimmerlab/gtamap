package mapper

import (
	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/index"
)

// Filter the read
// TODO: finalize the filter step as this is the most crucial step in terms of runtime
func Filter(readSequence *[]byte, genomeIndex *index.GenomeIndex) bool {
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

	//if numMatchingFw > 0 || numMatchingRw > 0 {
	//	logrus.WithFields(logrus.Fields{
	//		"numMatchingFw": numMatchingFw,
	//		"numMatchingRv": numMatchingRw,
	//	}).Info("Filtering result")
	//}

	// fmt.Println("numMatchingFw: ", numMatchingFw)
	// fmt.Println("numMatchingRw: ", numMatchingRw)

	return numMatching >= 6
}
