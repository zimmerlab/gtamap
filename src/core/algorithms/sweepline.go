package algorithms

import (
	"sort"
)

type Event struct {
	Position int
	Delta    int
}

type ByPosition []Event

func (a ByPosition) Len() int      { return len(a) }
func (a ByPosition) Swap(i, j int) { a[i], a[j] = a[j], a[i] }
func (a ByPosition) Less(i, j int) bool {
	if a[i].Position == a[j].Position {
		return a[i].Delta > a[j].Delta
	}
	return a[i].Position < a[j].Position
}

func FindHighCoverageRegions(intervals [][2]int) [][2]int {
	var events []Event
	for _, interval := range intervals {
		start, end := interval[0], interval[1]
		events = append(events, Event{Position: start, Delta: +1})
		events = append(events, Event{Position: end, Delta: -1})
	}

	sort.Sort(ByPosition(events))

	var result [][2]int
	var regionStart *int
	coverage := 0
	coverageSlice := make([]int, 0)
	totalCoverage := 0

	for _, e := range events {
		coverage += e.Delta
		totalCoverage += coverage
		coverageSlice = append(coverageSlice, coverage)
	}

	// TODO: add multiplyer as param
	threshold := totalCoverage / len(coverageSlice) * 2
	//threshold := 70
	for _, e := range events {
		coverage += e.Delta
		if coverage >= threshold && regionStart == nil {
			regionStart = &e.Position
		}
		if coverage < threshold && regionStart != nil {
			if e.Position > *regionStart {
				result = append(result, [2]int{*regionStart, e.Position})
			}
			regionStart = nil
		}
	}

	return result
}
