package config

import (
	"os"
	"path/filepath"

	"github.com/sirupsen/logrus"
)

var env string = "development"

const toolVersion string = "0.3.1"

var kmerLength uint8 = 10

// the maximum percentage of mismatches allowed in a read (50 = 50% percent mismatches allowed)
var maxMismatchPercentage uint8 = 10

// the minimum length of an intron (in base pairs)
// used to decide whether a gap is a deletion (if below this length) or an intron
var intronLengthMin = 20

// the maximum error rate allowed per read
var errorRate float64 = 0.05

var outputDirectory string = "~/gtamap-output"

var (
	includeReadsImproperlyPaired  bool = false
	includeReadsAmbiguouslyMapped bool = false
	includeReadsUnmapped          bool = false
)

var (
	outputIncludeReadSequence  bool = true
	outputIncludeReadQuality   bool = true
	outputIncludeDetailedCigar bool = false
)

// properPairLevel sets the level of reads that are considered properly paired
// 0 = both reads are mapped to the same reference on different strands
// 1 = 0 and the distance between the reads is consistent with the fragment length
var properPairLevel uint8 = 0

var (
	fragmentLength          int = 0
	fragmentLengthDeviation int = 0
)

// Mapping params
var (
	MaxBranchPoints int = 10 // how often if applyPossibleDiagonals allowed to branch
)

// ConfidentWorkerParams
var (
	MaxConfMm                 int = 6   // how many mm is a conf map allowed to have
	MinConfAnchorLengthRNA    int = 20  // how long does each ali block in a conf map have to be to be considered conf in RNA
	MinConfAnchorLengthDNA    int = 50  // how long does each ali block in a conf map have to be to be considered conf in DNA
	IntronClusterDelta        int = 100 // by default, an intron cluster only absorbes an incoming gap (extending its reach) if the delta of gap.start/stop and cluster.start/stop is less than 100. This allows overlapping introns but also resolves intron coord confilct within close proximity
	IntronClusterRepairWindow int = 5
)

// SAM options
var (
	IncludeMMinSAM     bool = true  // if set to true, CIGAR will include "=" and "X" runes instead of only "M"
	IncludeAllPairings bool = false // if set to true, do all vs all in SAM out (we can later implement a method which does that by also looling at tlen, etc)
)

// RNA/DNA Flag
var IsOriginRNA bool = true

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

func OutputDirectory() string {
	// check if directory exists and create it if not
	if outputDirectory == "" {
		outputDirectory = "~/gtamap-output"
	}
	if outputDirectory[0] == '~' {
		homeDir, err := os.UserHomeDir()
		if err != nil {
			logrus.Fatal("Error getting home directory", err)
		}
		outputDirectory = filepath.Join(homeDir, outputDirectory[1:])
	}
	if _, err := os.Stat(outputDirectory); os.IsNotExist(err) {
		err := os.MkdirAll(outputDirectory, os.ModePerm)
		if err != nil {
			logrus.Fatal("Error creating output directory", err)
		}
	}
	logrus.WithFields(logrus.Fields{
		"outputDirectory": outputDirectory,
	}).Info("Using output directory")
	return outputDirectory
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
