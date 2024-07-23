package interval

import (
	"sort"
)

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

func SortIntervals(intervals []*Interval) {
	sort.Slice(intervals, func(i, j int) bool {
		return intervals[i].Start < intervals[j].Start
	})
}

func MergeIntervals(intervals []*Interval) []*Interval {

	if len(intervals) <= 1 {
		return intervals
	}

	SortIntervals(intervals)

	merged := make([]*Interval, 1)
	merged[0] = &Interval{Start: intervals[0].Start, End: intervals[0].End}

	for _, nextInterval := range intervals[1:] {

		if nextInterval.Start <= merged[len(merged)-1].End {
			if nextInterval.End > merged[len(merged)-1].End {
				merged[len(merged)-1].End = nextInterval.End
			}
		} else {
			merged = append(merged, &Interval{Start: nextInterval.Start, End: nextInterval.End})
		}
	}

	return merged
}
