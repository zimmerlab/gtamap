package confidentmappingpass

import (
	"strconv"
	"strings"

	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"

	"github.com/KleinSamuel/gtamap/src/formats/fastq"
)

type ConfidentMappingTask struct {
	ReadPair *fastq.ReadPair
	ResultFw *mapperutils.ReadMatchResult
	ResultRv *mapperutils.ReadMatchResult
	Index    *index.GenomeIndex
}

func (i ConfidentMappingTask) String() string {
	var builder strings.Builder
	builder.Write([]byte("ReadPairR1 Header: "))
	builder.Write([]byte(i.ReadPair.ReadR1.Header))
	builder.Write([]byte("\n"))
	builder.Write([]byte("  <== FW MAPPING ==>"))
	builder.Write([]byte("\n"))
	builder.Write([]byte("\t SeqIndex: "))
	mapping := i.ResultFw
	seqIndex := strconv.Itoa(mapping.SequenceIndex)
	builder.WriteString(seqIndex)
	builder.WriteString("\n")
	builder.WriteString("\t GENOME -> ")
	builder.WriteString(mapping.MatchedGenome.String())
	builder.WriteString("\n")
	builder.WriteString("\t READ   -> ")
	builder.WriteString(mapping.MatchedRead.String())
	builder.WriteString("\n")
	builder.WriteString("\t MISMAT -> ")
	ints := mapping.MismatchesRead
	strs := make([]string, len(ints))
	for i, v := range ints {
		strs[i] = strconv.Itoa(v)
	}
	builder.WriteString(strings.Join(strs, ","))
	builder.WriteString("\n")
	builder.Write([]byte("  <== RV MAPPING ==>"))
	builder.Write([]byte("\n"))
	builder.Write([]byte("\t SeqIndex: "))
	mapping = i.ResultRv
	seqIndex = strconv.Itoa(mapping.SequenceIndex)
	builder.WriteString(seqIndex)
	builder.WriteString("\n")
	builder.WriteString("\t GENOME -> ")
	builder.WriteString(mapping.MatchedGenome.String())
	builder.WriteString("\n")
	builder.WriteString("\t READ   -> ")
	builder.WriteString(mapping.MatchedRead.String())
	builder.WriteString("\n")
	builder.WriteString("\t MISMAT -> ")
	ints = mapping.MismatchesRead
	strs = make([]string, len(ints))
	for i, v := range ints {
		strs[i] = strconv.Itoa(v)
	}
	builder.WriteString(strings.Join(strs, ","))
	builder.WriteString("\n")
	return builder.String()
}
