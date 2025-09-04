package main

import (
	"fmt"

	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/interval"
)

func iTest() {

	l := interval.NewList()

	l.InsertInterval(1, 3, 1, "a")
	l.InsertInterval(6, 10, 1, "a")
	l.InsertInterval(10, 11, 1, "a")
	l.InsertInterval(11, 20, 1, "a")
	l.InsertInterval(30, 40, 1, "a")

	l.Print()

	// os.Exit(1)

	// p := 11

	// index := l.ContainedIn(p)
	// fmt.Println("result", index)

	// i, offset := l.GetInsertionPos(p)
	// fmt.Println("insert at ", i, "with offset", offset)

	l.InsertInterval(4, 7, 2, "b")

	fmt.Println("after")
	l.Print()

	i := l.GetMostImportantItemThatOverlaps(5, 28)
	fmt.Println("most important", i)
}

func main() {

	// bedFilePath := "/home/sam/Data/gtamap/regionmask/mask.bed"
	//
	// bedFile, err := os.Open(bedFilePath)
	// if err != nil {
	// 	logrus.Fatal(err)
	// }
	//
	// b := bed.New(bedFile)
	//
	// prioFilePath := "/home/sam/Data/gtamap/regionmask/mask.priority"
	//
	// p, errP := bed.ReadPriorities(prioFilePath)
	// if errP != nil {
	// 	logrus.Fatal(errP)
	// }
	//
	// fmt.Println(len(b.Entries))
	// fmt.Println(len(p))
	// fmt.Println(p)
	//
	// interval.NewListWithData(p, b.NameMap["X"])

	data := []*interval.Interval{
		&interval.Interval{Start: 1, End: 3},
		&interval.Interval{Start: 4, End: 5},
		&interval.Interval{Start: 4, End: 7},
		&interval.Interval{Start: 7, End: 10},
		&interval.Interval{Start: 7, End: 10},
		&interval.Interval{Start: 11, End: 20},
		&interval.Interval{Start: 13, End: 15},
		&interval.Interval{Start: 13, End: 17},
		&interval.Interval{Start: 20, End: 31},
	}

	// l := interval.NewIntervalListFromData(data)
	//
	// l.Print()
	//
	// fmt.Println(l.GetFirstIndexContaining(5))

	rv := regionvector.NewRegionVector()

	for _, d := range data {
		rv.AddRegionAndMerge(d.Start, d.End)
	}

	fmt.Println(rv.StringTable())

	fmt.Println(rv.OverlapsAny(10, 12))

}
