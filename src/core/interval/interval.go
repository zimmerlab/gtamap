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

// GenerateEquivalenceClasses generates equivalence classes from a list of intervals.
// The given list of intervals can be unsorted and overlapping.
// The function returns a list of non-overlapping intervals that represent the equivalence classes.
// An equivalence class is a distinct interval which is defined by start and end positions of the given intervals.
// Idea:
// 1. Extract a list of positions where the given intervals start
// 2. Extract a list of positions where the given intervals end
// 3. Sort both lists
// 4. Iterate over both lists and consider the smallest value of both lists
// 5. If the smallest value is an opening, create a new interval and increment the number of open intervals
// 6. If the smallest value is a closing, create a new interval and decrement the number of open intervals
// 7. If the smallest value is an opening and closing at the same time, create a new interval and do not increment
// the number of open intervals
func GenerateEquivalenceClasses(intervals []*Interval) []*Interval {

	openings := make([]int, 0)
	closings := make([]int, 0)

	for _, i := range intervals {
		openings = append(openings, i.Start)
		closings = append(closings, i.End)
	}

	sort.Ints(openings)
	sort.Ints(closings)

	ecs := make([]*Interval, 0)
	numOpen := 0
	indexOpen := 0
	indexClose := 0
	lastPos := 0

	for indexOpen < len(openings) || indexClose < len(closings) {

		// only closes left
		if indexOpen >= len(openings) {

			if closings[indexClose] == lastPos {
				indexClose++
				continue
			}

			ecs = append(ecs, &Interval{
				Start: lastPos,
				End:   closings[indexClose],
			})
			lastPos = closings[indexClose]
			indexClose++
			numOpen--
			continue
		}

		if openings[indexOpen] == closings[indexClose] {
			// open and close are the same
			// create a new interval and update the last position
			// increment both indices
			// but do not increment the number of open intervals as one is opening and closing at the same time

			// if the opening is at the last position, skip
			if openings[indexOpen] == lastPos {
				indexOpen++
				indexClose++
				continue
			}

			ecs = append(ecs, &Interval{
				Start: lastPos,
				End:   openings[indexOpen],
			})

			lastPos = openings[indexOpen]
			indexOpen++
			indexClose++

		} else if openings[indexOpen] < closings[indexClose] {

			if openings[indexOpen] == lastPos {
				indexOpen++
				continue
			}

			if numOpen > 0 {
				// already open -> create new interval
				ecs = append(ecs, &Interval{
					Start: lastPos,
					End:   openings[indexOpen],
				})
			}

			lastPos = openings[indexOpen]
			indexOpen++
			numOpen++

		} else {

			if closings[indexClose] == lastPos {
				indexClose++
				continue
			}

			ecs = append(ecs, &Interval{
				Start: lastPos,
				End:   closings[indexClose],
			})

			lastPos = closings[indexClose]
			indexClose++
			numOpen--
		}
	}

	return ecs
}
