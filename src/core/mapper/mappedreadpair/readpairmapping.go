package mappedreadpair

import (
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
)

type ReadPairMatchResult struct {
	ReadPair *fastq.ReadPair
	Fw       *mapperutils.ReadMatchResult
	Rv       *mapperutils.ReadMatchResult
	Index    *index.GenomeIndex
}
