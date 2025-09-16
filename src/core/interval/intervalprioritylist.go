package interval

import (
	"fmt"
	"sort"

	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
)

type PriorityListItem struct {
	Start     int
	End       int
	Priority  int
	Threshold int
	Name      string
}

type PriorityList struct {
	Items       []*PriorityListItem
	MaxPriority int
}

func NewList() *PriorityList {
	return &PriorityList{
		Items:       make([]*PriorityListItem, 0),
		MaxPriority: -1,
	}
}

func NewListWithData(
	// priorities map[string]int,
	entries []*PriorityListItem,
) *PriorityList {

	pList := NewList()

	// sort bed file entries by priorities descending
	sort.Slice(entries, func(i int, j int) bool {
		return entries[i].Priority > entries[j].Priority
	})

	for _, interval := range entries {
		pList.InsertInterval(
			interval.Start,
			interval.End,
			interval.Priority,
			interval.Name,
			interval.Threshold,
		)
	}

	return pList
}

func (l *PriorityList) ContainedIn(position int) int {

	left := 0
	right := len(l.Items) - 1
	mid := 0

	for {

		if left > right || left < 0 || right >= len(l.Items) {
			return -1
		}

		mid = left + ((right - left) / 2)

		fmt.Println("check", mid, "with", l.Items[mid])

		if position < l.Items[mid].Start {
			right = mid - 1
			fmt.Println("go left")
		} else if position > l.Items[mid].End {
			left = mid + 1
			fmt.Println("go right")
		} else {
			// right = mid - 1
			// fmt.Println("found but go left")
			return mid
		}
	}
}

func (l *PriorityList) ApplyMaskToRegion(start int, end int) {

	regions := make([]regionvector.Region, 0)
	labels := make([]string, 0)

	index, offset := l.GetInsertionPos(start)

	fmt.Println(index, offset)

	// the given region is completely unmasked when:
	// - there are no items in the mask list
	// - the insertion index is after the last item
	// - the given region is between list items (offset -1 means not contained
	//   in item at index and should be inserted before)
	if len(l.Items) == 0 ||
		index >= len(l.Items) ||
		(offset == -1 && end <= l.Items[index].Start) {
		regions = append(regions, regionvector.Region{Start: start, End: end})
		labels = append(labels, "unmasked")
		return
	}

	item := l.Items[index]

	// add unmasked if before first item
	if offset == -1 {
		// there is an overlap with item or otherwise the function would have
		// returned above
		regions = append(regions, regionvector.Region{
			Start: start,
			End:   item.Start,
		})
		labels = append(labels, "unmasked")
		index++
	}

	for {

		regions = append(regions, regionvector.Region{
			Start: max(l.Items[index].Start, start),
			End:   min(l.Items[index].End, end),
		})
		labels = append(labels, l.Items[index].Name)

		index++

		if index >= len(l.Items) || l.Items[index].Start >= end {
			break
		}
	}

	fmt.Println("after regionmask is applied")

	for i := 0; i < len(regions); i++ {
		fmt.Println(regions[i], labels[i])
	}

	fmt.Println("-")
}

func (l *PriorityList) GetItemAtPosition(pos int) (bool, string, int) {

	index, offset := l.GetInsertionPos(pos)

	if len(l.Items) == 0 ||
		offset == -1 ||
		index >= len(l.Items) {
		return false, "", 0
	}

	return true, l.Items[index].Name, l.Items[index].Priority
}

func (l *PriorityList) GetMostImportantItemThatOverlaps(
	start int,
	end int,
) (bool, string, int) {

	index, offset := l.GetInsertionPos(start)

	if (offset != -1 && index == 0) || index >= len(l.Items) {
		return false, "", 0
	}

	item := l.Items[index]

	// fmt.Println("region\t", start, end)
	// fmt.Println("item\t", item.Start, item.End)

	prio := -1
	size := 0
	if offset > -1 {

		// fmt.Println(min(item.End, end))
		// fmt.Println("item.Start", item.Start, "offset", offset)

		prio = item.Priority
		size = min(item.End, end) - (item.Start + offset)
	}

	for {
		index++

		if index >= len(l.Items) || l.Items[index].Start >= end {
			break
		}

		// skip if item does not overlap region
		if l.Items[index].End <= start {
			continue
		}

		if prio == -1 || l.Items[index].Priority > prio {

			prio = l.Items[index].Priority
			size = min(l.Items[index].End, end) - max(l.Items[index].Start, start)
			item = l.Items[index]

			// if item.Priority == l.MaxPriority {
			// 	// premature exit if max priority reached
			// 	return item
			// }

		} else if l.Items[index].Priority == prio {
			size += min(l.Items[index].End, end) - max(l.Items[index].Start, start)
		}
	}

	if size == 0 {
		return false, "", 0
	}

	return true, item.Name, size
}

func (l *PriorityList) GetInsertionPos(position int) (int, int) {

	left := 0
	right := len(l.Items) - 1
	mid := 0

	for {

		if left > right || left < 0 || right >= len(l.Items) {
			dir := 0
			if mid < len(l.Items) && l.Items[mid].Start < position {
				dir = 1
			}
			return mid + dir, -1
		}

		mid = left + ((right - left) / 2)

		// fmt.Println("check", mid, "with", l.Items[mid])

		if position < l.Items[mid].Start {
			right = mid - 1
			// fmt.Println("go left")
		} else if position >= l.Items[mid].End {
			left = mid + 1
			// fmt.Println("go right")
		} else {
			return mid, position - l.Items[mid].Start
		}
	}
}

func (l *PriorityList) InsertInterval(
	start int,
	end int,
	priority int,
	name string,
	threshold int,
) {

	// fmt.Println("insert ", start, "-", end)

	if start >= end {
		l.MaxPriority = max(l.MaxPriority, priority)
		return
	}

	index, offset := l.GetInsertionPos(start)

	// fmt.Println(index, offset)

	if offset == -1 {
		// start in region with no interval

		// fmt.Println("start not contained. index", index)

		newEnd := end
		if index < len(l.Items) && l.Items[index].Start < end {
			newEnd = l.Items[index].Start
		}

		l.InsertAtIndex(index, start, newEnd, priority, name, threshold)

		l.InsertInterval(newEnd, end, priority, name, threshold)

	} else {
		// interval overlaps with existing interval

		oldEnd := l.Items[index].End
		oldPriority := l.Items[index].Priority
		oldName := l.Items[index].Name
		oldThreshold := l.Items[index].Threshold

		newStart := l.Items[index].Start + offset
		newEnd := min(oldEnd, end)

		// fmt.Println("overlap with", l.Items[index], "newStart", newStart, "newEnd", newEnd)

		if l.Items[index].Priority >= priority {
			l.InsertInterval(newEnd, end, priority, name, threshold)
			return
		}

		// split existing interval
		if offset == 0 {
			// convert item to new interval
			l.Items[index].End = newEnd
			l.Items[index].Priority = priority
			l.Items[index].Name = name
			l.Items[index].Threshold = threshold
		} else {
			// keep item but split
			l.Items[index].End = newStart

			index = index + 1

			l.InsertAtIndex(index, newStart, newEnd, priority, name, threshold)
		}

		if newEnd < oldEnd {
			l.InsertAtIndex(index+1, newEnd, oldEnd, oldPriority, oldName, oldThreshold)
		}

		l.InsertInterval(newEnd, end, priority, name, threshold)
	}
}

func (l *PriorityList) InsertAtIndex(
	index int,
	start int,
	end int,
	priority int,
	name string,
	threshold int,
) {
	if index < 0 || index > len(l.Items) {
		return
	}
	if index == len(l.Items) {
		l.Items = append(l.Items, &PriorityListItem{
			Start:     start,
			End:       end,
			Priority:  priority,
			Name:      name,
			Threshold: threshold,
		})
	} else {
		l.Items = append(l.Items[:index+1], l.Items[index:]...)
		l.Items[index] = &PriorityListItem{
			Start:     start,
			End:       end,
			Priority:  priority,
			Name:      name,
			Threshold: threshold,
		}
	}
}

func (l *PriorityList) Print() {
	fmt.Printf("\n%5s ", "index")
	for i, _ := range l.Items {
		fmt.Printf("%3d ", i)
	}
	fmt.Printf("\n%5s ", "start")
	for _, item := range l.Items {
		fmt.Printf("%3d ", item.Start)
	}
	fmt.Printf("\n%5s ", "end")
	for _, item := range l.Items {
		fmt.Printf("%3d ", item.End)
	}
	fmt.Printf("\n%5s ", "type")
	for _, item := range l.Items {
		fmt.Printf("%3s ", item.Name)
	}
	fmt.Printf("\n%5s ", "prio")
	for _, item := range l.Items {
		fmt.Printf("%3d ", item.Priority)
	}
	fmt.Printf("\n%5s ", "thres")
	for _, item := range l.Items {
		fmt.Printf("%3d ", item.Threshold)
	}
	fmt.Printf("\n\n")

}
