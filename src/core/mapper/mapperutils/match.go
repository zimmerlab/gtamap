package mapperutils

import (
	"fmt"
	"strconv"
	"strings"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/sirupsen/logrus"
)

type Match struct {
	SequenceIndex int // the index of the sequence in the genome
	FromGenome    int // the start position of the match in the genome
	ToGenome      int // the end position of the match in the genome
	FromRead      int // the start position of the match in the read
	ToRead        int // the end position of the match in the read
	StartGenome   int // the start position of the match in the genome (diagonal)
	Used          bool
}

func (m *Match) String() string {
	return fmt.Sprintf("[%d, %d, %d, %d, %t]", m.FromRead, m.ToRead, m.FromGenome, m.ToGenome, m.Used)
}

type GlobalMatchResult struct {
	MatchesPerSequence []*SequenceMatchResult
}

type SequenceMatchResult struct {
	MatchesPerDiagonal map[int][]*Match
}

type ReadMatchResult struct {
	SequenceIndex   int                        // the index of the sequence in the genome
	MatchedRead     *regionvector.RegionVector // region vector containing the matched positions in the read
	MatchedGenome   *regionvector.RegionVector // region vector containing the matched positions in the genome
	MismatchesRead  []int                      // the positions of the mismatches in the read
	diagonalHandler *DiagonalHandler
	IncompleteMap   bool
	NeedRemap       bool // true if this result must undergo a remap
}

func (m ReadMatchResult) Copy() *ReadMatchResult {
	var dhCopy *DiagonalHandler
	if m.diagonalHandler != nil {
		dhCopy = m.diagonalHandler.Copy()
	}

	return &ReadMatchResult{
		SequenceIndex:   m.SequenceIndex,
		MatchedRead:     m.MatchedRead.Copy(),
		MatchedGenome:   m.MatchedGenome.Copy(),
		MismatchesRead:  append([]int{}, m.MismatchesRead...),
		NeedRemap:       m.NeedRemap,
		diagonalHandler: dhCopy,
	}
}

func (m ReadMatchResult) GetCigar() (string, error) {
	var builder strings.Builder

	isForwardStrand := m.SequenceIndex == 0

	if isForwardStrand {

		// number of cumulative aligned positions (matches or mismatches) between the read and the genome
		numMatchesSum := 0

		// the regions in MatchedGenome and MatchedRead have the same dimensions
		for i := 0; i < len(m.MatchedGenome.Regions); i++ {

			// each pair of regions represent a match
			numMatchesSum += m.MatchedGenome.Regions[i].Length()

			// if this is the last match, add the number of matches
			if i == len(m.MatchedGenome.Regions)-1 {
				builder.WriteString(strconv.Itoa(numMatchesSum))
				builder.WriteString("M")
				break
			}

			gapInGenome := m.MatchedGenome.Regions[i].End < m.MatchedGenome.Regions[i+1].Start
			gapInRead := m.MatchedRead.Regions[i].End < m.MatchedRead.Regions[i+1].Start

			if gapInGenome && gapInRead {
				logrus.WithFields(logrus.Fields{
					"read":   m.MatchedRead,
					"genome": m.MatchedGenome,
				}).Warn("Gap in genome and read at the same time")

				return "", fmt.Errorf("gap in genome and read at the same time")
			}

			if !gapInGenome && !gapInRead {
				continue
			}

			builder.WriteString(strconv.Itoa(numMatchesSum))
			builder.WriteString("M")
			numMatchesSum = 0

			// intron or deletion
			if gapInGenome {
				// number of skipped bases in the reference
				numSkipped := m.MatchedGenome.Regions[i+1].Start - m.MatchedGenome.Regions[i].End

				builder.WriteString(strconv.Itoa(numSkipped))

				// determine whether the gap is an intron or a deletion
				if numSkipped < config.IntronLengthMin() {
					builder.WriteString("D")
				} else {
					builder.WriteString("N")
				}
			}

			// insertion
			if gapInRead {
				// number of skipped bases in the read
				numSkipped := m.MatchedRead.Regions[i+1].Start - m.MatchedRead.Regions[i].End

				builder.WriteString(strconv.Itoa(numSkipped))
				builder.WriteString("I")
			}
		}

	} else {

		// reverse the order of the regions in MatchedGenome
		// this is because the read is on the reverse strand

		numMatchesSum := 0

		// the regions in MatchedGenome and MatchedRead have the same dimensions
		for i := len(m.MatchedGenome.Regions) - 1; i >= 0; i-- {

			numMatchesSum += m.MatchedGenome.Regions[i].Length()

			// if this is the last match, add the number of matches
			if i == 0 {
				builder.WriteString(strconv.Itoa(numMatchesSum))
				builder.WriteString("M")
				break
			}

			gapInGenome := m.MatchedGenome.Regions[i].Start > m.MatchedGenome.Regions[i-1].End
			gapInRead := m.MatchedRead.Regions[i].Start > m.MatchedRead.Regions[i-1].End

			if gapInGenome && gapInRead {
				logrus.WithFields(logrus.Fields{
					"read":   m.MatchedRead,
					"genome": m.MatchedGenome,
				}).Warn("Gap in genome and read at the same time")

				return "", fmt.Errorf("gap in genome and read at the same time")
			}

			if !gapInGenome && !gapInRead {
				continue
			}

			builder.WriteString(strconv.Itoa(numMatchesSum))
			builder.WriteString("M")
			numMatchesSum = 0

			// intron or deletion
			if gapInGenome {
				// number of skipped bases in the reference
				numSkipped := m.MatchedGenome.Regions[i].Start - m.MatchedGenome.Regions[i-1].End

				builder.WriteString(strconv.Itoa(numSkipped))

				// determine whether the gap is an intron or a deletion
				if numSkipped < config.IntronLengthMin() {
					builder.WriteString("D")
				} else {
					builder.WriteString("N")
				}
			}

			// insertion
			if gapInRead {
				// number of skipped bases in the read
				numSkipped := m.MatchedRead.Regions[i].Start - m.MatchedRead.Regions[i-1].End

				builder.WriteString(strconv.Itoa(numSkipped))
				builder.WriteString("I")
			}

		}
	}

	return builder.String(), nil
}

// holds all relevant mapping information of potentially several hits per readpair
// bundled into one type
type ReadPairMatchResults struct {
	ReadPair *fastq.ReadPair
	Fw       []*ReadMatchResult
	Rv       []*ReadMatchResult
}

func (i ReadPairMatchResults) String() string {
	var builder strings.Builder
	builder.Write([]byte("ReadPairR1 Header: "))
	builder.Write([]byte(i.ReadPair.ReadR1.Header))
	builder.Write([]byte("\n"))
	builder.Write([]byte("  <== FW MAPPINGS ==>"))
	builder.Write([]byte("\n"))
	for _, mapping := range i.Fw {
		builder.Write([]byte("\t SeqIndex: "))
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
	}
	builder.Write([]byte("  <== RV MAPPINGS ==>"))
	builder.Write([]byte("\n"))
	for _, mapping := range i.Rv {
		builder.Write([]byte("\t SeqIndex: "))
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
	}
	return builder.String()
}
