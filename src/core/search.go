package core

import (
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
	"sort"
)

// TODO: (performance) remove pattern as it is probably already known and does not have to be stored twice
type PatternSearchResult struct {
	Pattern *string
	Matches []PatternMatch
}

type PatternMatch struct {
	// the index of the sequence within the trees sequences
	SequenceIndex int
	// the 0-based start position of the pattern within that sequence
	From int
	// the end-exclusive position of the pattern within that sequence
	To int
}

type SequenceMatch struct {
	SequenceIndex       int   // the index of the sequence within the trees sequences
	FromSource          int   // the 0-based start position of the match within the source sequence (read)
	ToSource            int   // the end-exclusive position of the match within the source sequence (read)
	FromTarget          int   // the 0-based start position of the match within the target sequence (transcript)
	ToTarget            int   // the end-exclusive position of the match within the target sequence (transcript)
	Mismatches          []int // the locations of the mismatches in the source sequence (read)
	EquivalenceClassIds []uint32
}

type MatchesSortableByFromSource []*SequenceMatch

func (m MatchesSortableByFromSource) Len() int {
	return len(m)
}
func (m MatchesSortableByFromSource) Less(i, j int) bool {
	return m[i].FromSource < m[j].FromSource
}
func (m MatchesSortableByFromSource) Swap(i, j int) {
	m[i], m[j] = m[j], m[i]
}
func (m MatchesSortableByFromSource) Sort() {
	sort.Sort(m)
}
func (m MatchesSortableByFromSource) Last() (SequenceMatch, bool) {
	if len(m) == 0 {
		return SequenceMatch{}, false
	}
	return *m[len(m)-1], true
}
func (m MatchesSortableByFromSource) Add(match *SequenceMatch) MatchesSortableByFromSource {
	return append(m, match)
}

type ExactMatchResult struct {
	Matches []*SequenceMatch
}

type DiscardStepMatchInformation struct {
	NumMismatches int
	Matches       MatchesSortableByFromSource
}

type SequenceMatches struct {
	SequenceIndex int                         // the index of the sequence within the target sequences
	NumMatches    int                         // the number of unique k-mer matchutils to this sequence
	Matches       MatchesSortableByFromSource // the matchutils to this sequence
}

type ReadMapResult struct {
	SequenceMatches map[int]SequenceMapResult // key = sequence index in suffix tree
}

type SequenceMapResult struct {
	SequenceContextMatches map[int]SequenceContextMatch // key = start index of read match on target sequence
}

type SequenceContextMatch struct {
	// the index of the sequence within the target sequences (contained in the index tree)
	SequenceIndex int
	// the 0-based start position of the match within the target sequence
	TargetIndex int
	Matches     MatchesSortableByFromSource
}

func (c SequenceContextMatch) GetUnmatchedRegions(readLength int) []SequenceMatch {

	unmatchedRegions := make([]SequenceMatch, 0)

	for i, match := range c.Matches {

		if i == 0 {
			if match.FromSource > 0 {
				unmatchedRegions = append(unmatchedRegions, SequenceMatch{
					SequenceIndex: c.SequenceIndex,
					FromSource:    0,
					ToSource:      match.FromSource,
					FromTarget:    c.TargetIndex,
					ToTarget:      c.TargetIndex + match.FromSource,
				})
			}

		} else {

			prevMatch := c.Matches[i-1]
			unmatchedRegions = append(unmatchedRegions, SequenceMatch{
				SequenceIndex: c.SequenceIndex,
				FromSource:    prevMatch.ToSource,
				ToSource:      match.FromSource,
				FromTarget:    c.TargetIndex + prevMatch.ToSource,
				ToTarget:      c.TargetIndex + match.FromSource,
			})
		}

		// the last match
		if i == len(c.Matches)-1 {
			// the end of the read is not matched yet
			if match.ToSource < readLength {
				unmatchedRegions = append(unmatchedRegions, SequenceMatch{
					SequenceIndex: c.SequenceIndex,
					FromSource:    match.ToSource,
					ToSource:      readLength,
					FromTarget:    c.TargetIndex + match.ToSource,
					ToTarget:      c.TargetIndex + readLength,
				})
			}
		}
	}

	return unmatchedRegions
}

type InexactMatchResult struct {
	SequenceIndex       int      // the id of the sequence within the suffix tree (transcript)
	FromTarget          int      // the 0-based start position of the match within the target sequence (transcript)
	ToTarget            int      // the end-exclusive position of the match within the target sequence (transcript)
	Mismatches          []int    // the locations of the mismatches in the source sequence (read)
	EquivalenceClassIds []uint32 // the ids of the equivalence classes that the match overlaps with
	FromGene            uint32   // the 0-based start position of the match within the gene
}

type ProperPairCandidate struct {
	ReferenceIndex int  // the index of the reference (transcript) determined by the sequence indices of R1 and R2
	FragmentLength int  // the fragment length of the pair (distance between the 5' ends of the reads)
	R1isForward    bool // true if R1 is on the forward strand
	ResultR1       *InexactMatchResult
	ResultR2       *InexactMatchResult
}

// ReadMappingPreResult is the result of mapping step where the reads are mapped to the suffix tree
// and each read can have multiple locations on multiple references.
type ReadMappingPreResult struct {
	ReadR1    *fastq.Read
	ReadR2    *fastq.Read
	ResultsR1 *map[int][]*InexactMatchResult
	ResultsR2 *map[int][]*InexactMatchResult
}

// ReadMappingResult is the result of the mapping step where the positions of the reads are determined
// using the relevant parameters (proper pair level, fragment length, etc.) and SAM record pairs are created.
type ReadMappingResult struct {
	ReadR1         *fastq.Read
	ReadR2         *fastq.Read
	SamRecordPairs []*sam.RecordPair
}
