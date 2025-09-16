package config

var env string = "development"

const toolVersion string = "0.4.0"

var kmerLength uint8 = 10

// SAM options
var (
	IncludeMMinSAM     bool = true  // if set to true, CIGAR will include "=" and "X" runes instead of only "M"
	IncludeAllPairings bool = false // if set to true, do all vs all in SAM out (we can later implement a method which does that by also looling at tlen, etc)
)

// RNA/DNA Flag
var IsOriginRNA bool = true

var LogOut string = "gta_log.tsv"

func Env() string {
	return env
}

func ToolVersion() string {
	return toolVersion
}

func KmerLength() uint8 {
	return kmerLength
}
