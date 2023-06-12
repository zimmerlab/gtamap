package core

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
