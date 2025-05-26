package confidentmappingpass

import (
	"strconv"
	"strings"

	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"

	"github.com/KleinSamuel/gtamap/src/formats/fastq"
)

type ConfidentMappingTask struct {
	ReadPair *fastq.ReadPair
	ResultFw *mapperutils.ReadMatchResult
	ResultRv *mapperutils.ReadMatchResult
}

func (i ConfidentMappingTask) String() string {
	var builder strings.Builder
	builder.WriteString("ReadPairR1 Header: ")
	builder.WriteString(i.ReadPair.ReadR1.Header)
	builder.WriteString("\n")
	builder.WriteString("  <== FW MAPPING ==>")
	builder.WriteString("\n")
	builder.WriteString("\t SeqIndex: ")
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
	builder.WriteString("  <== RV MAPPING ==>")
	builder.WriteString("\n")
	builder.WriteString("\t SeqIndex: ")
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
