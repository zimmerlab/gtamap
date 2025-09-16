package interval

import (
	"container/list"
	"fmt"
	"sort"
)

type List struct {
	InsertList *list.List
	Elements   []*Interval
}

func NewIntervalList() *List {
	return &List{
		InsertList: list.New(),
	}
}

func NewIntervalListFromData(intervals []*Interval) *List {
	l := NewIntervalList()
	l.AddAll(intervals)
	l.ToSlice()
	l.SortSlice()
	return l
}

func (l *List) Add(i *Interval) {
	l.InsertList.PushBack(i)
}

func (l *List) AddAll(intervals []*Interval) {
	for _, i := range intervals {
		l.Add(i)
	}
}

func (l *List) ToSlice() {
	l.Elements = make([]*Interval, 0, l.InsertList.Len())
	for e := l.InsertList.Front(); e != nil; e = e.Next() {
		if interval, ok := e.Value.(*Interval); ok {
			l.Elements = append(l.Elements, interval)
		}
	}

	l.InsertList = nil
}

func (l *List) SortSlice() {
	sort.Slice(l.Elements, func(i, j int) bool {
		if l.Elements[i].Start == l.Elements[j].Start {
			return l.Elements[i].End < l.Elements[j].End
		}
		return l.Elements[i].Start < l.Elements[j].Start
	})
}

func (l *List) GetFirstIndexContaining(position int) int {

	left := 0
	right := len(l.Elements) - 1
	mid := 0

	for {

		fmt.Printf("\nleft: %d, right: %d, mid: %d\n", left, right, mid)

		if left > right || left < 0 || right >= len(l.Elements) {
			// dir := 0
			// if mid < len(l.Elements) && l.Elements[mid].Start < position {
			// 	dir = 1
			// }
			// return mid + dir
			return -1
		}

		mid = left + ((right - left) / 2)

		fmt.Printf("checking mid %d (%d-%d)\n", mid, l.Elements[mid].Start, l.Elements[mid].End)

		if position < l.Elements[mid].Start {
			right = mid - 1
			fmt.Println("go left")
		} else if position >= l.Elements[mid].End {
			left = mid + 1
			fmt.Println("go right")
		} else {
			// return mid, position - l.Elements[mid].Start
			fmt.Println("found")
			return mid
		}
	}
}

func (l *List) Print() {
	fmt.Printf("\n%5s ", "index")
	for i, _ := range l.Elements {
		fmt.Printf("%3d ", i)
	}
	fmt.Printf("\n%5s ", "start")
	for _, item := range l.Elements {
		fmt.Printf("%3d ", item.Start)
	}
	fmt.Printf("\n%5s ", "end")
	for _, item := range l.Elements {
		fmt.Printf("%3d ", item.End)
	}
	fmt.Printf("\n\n")
}
