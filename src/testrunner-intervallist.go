package main

import (
	"fmt"

	"github.com/KleinSamuel/gtamap/src/core/interval"
)

func main() {

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
