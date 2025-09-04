package interval

import (
	"fmt"
	"sort"
)

type PriorityListItem struct {
	Start    int
	End      int
	Priority int
	Name     string
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

func NewListWithData(priorities map[string]int, entries []*PriorityListItem) *PriorityList {

	pList := NewList()

	// sort bed file entries by priorities descending
	sort.Slice(entries, func(i int, j int) bool {
		pI := 0
		if p, exists := priorities[entries[i].Name]; exists {
			pI = p
		}
		pJ := 0
		if p, exists := priorities[entries[j].Name]; exists {
			pJ = p
		}
		return pI > pJ
	})

	for _, interval := range entries {
		pList.InsertInterval(interval.Start, interval.End, priorities[interval.Name], interval.Name)
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

func (l *PriorityList) GetMostImportantItemThatOverlaps(start int, end int) *PriorityListItem {

	index, offset := l.GetInsertionPos(start)

	if offset != -1 && index == 0 {
		return nil
	}

	item := l.Items[index]

	for {
		index++

		if index >= len(l.Items) || l.Items[index].Start >= end {
			break
		}

		if l.Items[index].Priority > item.Priority {
			item = l.Items[index]

			if item.Priority == l.MaxPriority {
				// premature exit if max priority reached
				return item
			}
		}
	}

	return item
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

func (l *PriorityList) InsertInterval(start int, end int, priority int, name string) {

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

		l.InsertAtIndex(index, start, newEnd, priority, name)

		l.InsertInterval(newEnd, end, priority, name)

	} else {
		// interval overlaps with existing interval

		oldEnd := l.Items[index].End
		oldPriority := l.Items[index].Priority
		oldName := l.Items[index].Name

		newStart := l.Items[index].Start + offset
		newEnd := min(oldEnd, end)

		// fmt.Println("overlap with", l.Items[index], "newStart", newStart, "newEnd", newEnd)

		if l.Items[index].Priority >= priority {
			l.InsertInterval(newEnd, end, priority, name)
			return
		}

		// split existing interval
		if offset == 0 {
			// convert item to new interval
			l.Items[index].End = newEnd
			l.Items[index].Priority = priority
			l.Items[index].Name = name
		} else {
			// keep item but split
			l.Items[index].End = newStart

			index = index + 1

			l.InsertAtIndex(index, newStart, newEnd, priority, name)
		}

		if newEnd < oldEnd {
			l.InsertAtIndex(index+1, newEnd, oldEnd, oldPriority, oldName)
		}

		l.InsertInterval(newEnd, end, priority, name)
	}
}

func (l *PriorityList) InsertAtIndex(index int, start int, end int, priority int, name string) {
	if index < 0 || index > len(l.Items) {
		return
	}
	if index == len(l.Items) {
		l.Items = append(l.Items, &PriorityListItem{
			Start:    start,
			End:      end,
			Priority: priority,
			Name:     name,
		})
	} else {
		l.Items = append(l.Items[:index+1], l.Items[index:]...)
		l.Items[index] = &PriorityListItem{
			Start:    start,
			End:      end,
			Priority: priority,
			Name:     name,
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
	fmt.Printf("\n\n")

}
