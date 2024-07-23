package config

var env string = "development"

const toolVersion string = "0.1"

var kmerLength int = 8

var includeReadsImproperlyPaired bool = false
var includeReadsAmbiguouslyMapped bool = false
var includeReadsUnmapped bool = false

var outputIncludeReadSequence bool = true
var outputIncludeReadQuality bool = true
var outputIncludeDetailedCigar bool = false

// properPairLevel sets the level of reads that are considered properly paired
// 0 = both reads are mapped to the same reference on different strands
// 1 = 0 and the distance between the reads is consistent with the fragment length
var properPairLevel uint8 = 0

var fragmentLength int = 0
var fragmentLengthDeviation int = 0

func Env() string {
	return env
}

func ToolVersion() string {
	return toolVersion
}

func KmerLength() int {
	return kmerLength
}

func IncludeReadsImproperlyPaired() bool {
	return includeReadsImproperlyPaired
}

func IncludeReadsAmbiguouslyMapped() bool {
	return includeReadsAmbiguouslyMapped
}

func IncludeReadsUnmapped() bool {
	return includeReadsUnmapped
}

func SetOutputIncludeReadSequence(value bool) {
	outputIncludeReadSequence = value
}
func OutputIncludeReadSequence() bool {
	return outputIncludeReadSequence
}

func SetOutputIncludeReadQuality(value bool) {
	outputIncludeReadQuality = value
}
func OutputIncludeReadQuality() bool {
	return outputIncludeReadQuality
}

func SetOutputIncludeDetailedCigar(value bool) {
	outputIncludeDetailedCigar = value
}
func OutputIncludeDetailedCigar() bool {
	return outputIncludeDetailedCigar
}

func ProperPairLevel() uint8 {
	return properPairLevel
}
