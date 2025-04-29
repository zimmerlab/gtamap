package config

var env string = "development"

const toolVersion string = "0.2"

var kmerLength uint8 = 10

// the maximum percentage of mismatches allowed in a read (50 = 50% percent mismatches allowed)
var maxMismatchPercentage uint8 = 50

// the minimum length of an intron (in base pairs)
// used to decide whether a gap is a deletion (if below this length) or an intron
var intronLengthMin = 20

// the maximum error rate allowed per read
var errorRate float64 = 0.05

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

func ErrorRate() float64 {
	return errorRate
}

func KmerLength() uint8 {
	return kmerLength
}

func MaxMismatchPercentage() uint8 {
	return maxMismatchPercentage
}

func IntronLengthMin() int {
	return intronLengthMin
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
