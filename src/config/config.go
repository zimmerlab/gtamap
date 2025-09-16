package config

var env string = "development"

const toolVersion string = "0.4.0"

var kmerLength uint8 = 10

// Mapping params
var (
	MaxBranchPoints int = 30 // how often if applyPossibleDiagonals allowed to branch
)

// ConfidentWorkerParams
var (
	MaxConfMm                 int = 6   // how many mm is a conf map allowed to have
	MinConfAnchorLengthRNA    int = 20  // how long does each ali block in a conf map have to be considered conf in RNA
	MinConfAnchorLengthDNA    int = 50  // how long does each ali block in a conf map have to be considered conf in DNA
	IntronClusterDelta        int = 100 // by default, an intron cluster only absorbes an incoming gap (extending its reach) if the delta of gap.start/stop and cluster.start/stop is less than 100. This allows overlapping introns but also resolves intron coord confilct within close proximity
	IntronClusterRepairWindow int = 10
)

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
