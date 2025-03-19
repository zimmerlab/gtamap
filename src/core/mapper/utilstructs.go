package mapper

import (
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
)

type MappingTask struct {
	ID       int
	ReadPair *fastq.ReadPair
}
