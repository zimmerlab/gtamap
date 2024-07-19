package core

import (
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
	// the index of the sequence within the trees sequences
	SequenceIndex int
	// the 0-based start position of the match within the source sequence (read)
	FromSource int
	// the end-exclusive position of the match within the source sequence (read)
	ToSource int
	// the 0-based start position of the match within the target sequence (transcript)
	FromTarget int
	// the end-exclusive position of the match within the target sequence (transcript)
	ToTarget int
	// the locations of the mismatches in the source sequence (read)
	Mismatches []int
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
	NumMatches    int                         // the number of unique k-mer matches to this sequence
	Matches       MatchesSortableByFromSource // the matches to this sequence
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

type Interval struct {
	Start int
	End   int
}

type Intervals []Interval

func (intervals Intervals) Len() int {
	return len(intervals)
}
func (intervals Intervals) Less(i, j int) bool {
	return intervals[i].Start < intervals[j].Start
}
func (intervals Intervals) Swap(i, j int) {
	intervals[i], intervals[j] = intervals[j], intervals[i]
}

type InexactMatchResult struct {
	SequenceIndex int   // the id of the sequence within the suffix tree (transcript)
	FromTarget    int   // the 0-based start position of the match within the target sequence
	ToTarget      int   // the end-exclusive position of the match within the target sequence
	Mismatches    []int // the locations of the mismatches in the source sequence (read)
}
