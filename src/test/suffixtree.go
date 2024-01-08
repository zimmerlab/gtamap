package test

import (
	"bufio"
	"container/list"
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core/datastructure"
	"log"
	"os"
	"slices"
	"strconv"
	"strings"
)

type EdgeListItem struct {
	From  int
	To    int
	Label string
}

func ReadEdgeList(path string) (string, []EdgeListItem) {
	items := make([]EdgeListItem, 0)

	file, err := os.Open(path)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	pattern := ""

	for scanner.Scan() {
		lineArray := strings.Fields(scanner.Text())

		if lineArray[0] == "#" {
			pattern = lineArray[1]
			continue
		}

		from, errFrom := strconv.Atoi(lineArray[0])
		to, errTo := strconv.Atoi(lineArray[1])
		label := lineArray[2]

		if errFrom != nil || errTo != nil {
			log.Fatal(err)
		}
		items = append(items, EdgeListItem{from, to, label})
	}

	return pattern, items
}

func TestTree(edgelistFile string) bool {

	pattern, _ := ReadEdgeList(edgelistFile)

	tree := datastructure.CreateNewTree()
	tree.AddSequence(pattern, 0)

	treeEdgeList := make([]EdgeListItem, 0)

	queue := list.New()
	queue.PushBack(tree.Root)
	visited := make(map[int]bool)

	for queue.Len() > 0 {

		item := queue.Front()
		queue.Remove(item)
		node := item.Value.(*datastructure.Node)
		visited[node.Id] = true

		for _, edge := range node.Edges {

			treeEdgeList = append(treeEdgeList, EdgeListItem{node.Id, edge.To.Id, tree.GetEdgeSequence(edge)})

			if !visited[edge.To.Id] {
				queue.PushBack(edge.To)
			}
		}

		if node.Link != nil {
			treeEdgeList = append(treeEdgeList, EdgeListItem{node.Id, node.Link.Id, "link"})
		}
	}

	passed := true

	fmt.Println("Tree edges:")
	for _, edge := range treeEdgeList {
		fmt.Println(edge)

		if slices.Contains() {
			
		}
	}

	return passed
}
