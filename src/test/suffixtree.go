package test

import (
	"bufio"
	"fmt"
	"github.com/sirupsen/logrus"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

type Edge struct {
	From  int
	To    int
	Label string
}

type SortableEdges []Edge

func (s SortableEdges) Len() int {
	return len(s)
}

func (s SortableEdges) Less(i, j int) bool {
	if s[i].From < s[j].From {
		return true
	}
	if s[i].From > s[j].From {
		return false
	}
	if s[i].To < s[j].To {
		return true
	}
	if s[i].To > s[j].To {
		return false
	}
	return s[i].Label < s[j].Label
}

func (s SortableEdges) Swap(i, j int) {
	s[i], s[j] = s[j], s[i]
}

type Position struct {
	Index int
	Start int
}

type SortablePositions []Position

func (s SortablePositions) Len() int {
	return len(s)
}

func (s SortablePositions) Less(i, j int) bool {
	if s[i].Index < s[j].Index {
		return true
	}
	if s[i].Index > s[j].Index {
		return false
	}
	return s[i].Start < s[j].Start
}

func (s SortablePositions) Swap(i, j int) {
	s[i], s[j] = s[j], s[i]
}

type Testcase struct {
	Sequences []string
	Edges     []Edge
	Positions map[int][]Position
}

func (t *Testcase) Print() {
	fmt.Println("Sequences:")
	for i, sequence := range t.Sequences {
		fmt.Println(i, sequence)
	}
	fmt.Println("Edges:")
	for _, edge := range t.Edges {
		fmt.Println(edge)
	}
	fmt.Println("Positions:")
	for key, positions := range t.Positions {
		fmt.Println(key, positions)
	}
}

func (t *Testcase) IsEqual(o *Testcase) bool {

	logrus.Debug("Compare testcases")

	// compare sequences
	if !compareSequences(t.Sequences, o.Sequences) {
		return false
	}

	// compare edges
	if !compareEdges(t.Edges, o.Edges) {
		return false
	}

	// compare positions
	if !comparePositions(t.Positions, o.Positions) {
		return false
	}

	logrus.Debug("Testcases are equal")

	return true
}

func compareSequences(truth []string, computed []string) bool {

	if len(truth) != len(computed) {
		logrus.WithFields(logrus.Fields{
			"truth":    len(truth),
			"computed": len(computed),
		}).Warn("different number of sequences!")
		return false
	}

	sort.Strings(truth)
	sort.Strings(computed)

	for i, sequence := range truth {
		if sequence != truth[i] {
			logrus.WithFields(logrus.Fields{
				"truth":    sequence,
				"computed": computed[i],
			}).Warn("sequences are not equal!")
			return false
		}
	}
	return true
}

func compareEdges(truth []Edge, computed []Edge) bool {

	if len(truth) != len(computed) {
		logrus.WithFields(logrus.Fields{
			"truth":    len(truth),
			"computed": len(computed),
		}).Warn("different number of edges!")
		return false
	}

	sort.Sort(SortableEdges(truth))
	sort.Sort(SortableEdges(computed))

	for i, edge := range truth {
		if edge != computed[i] {
			logrus.WithFields(logrus.Fields{
				"edge truth":    edge,
				"edge computed": computed[i],
			}).Warn("edges are not equal!")
			return false
		}
	}
	return true
}

func comparePositions(truth map[int][]Position, computed map[int][]Position) bool {

	if len(truth) != len(computed) {
		logrus.WithFields(logrus.Fields{
			"len(truth)":    len(truth),
			"len(computed)": len(computed),
		}).Warn("different number of positions!")
		return false
	}

	for key, positions := range truth {
		if len(positions) != len(computed[key]) {
			logrus.WithFields(logrus.Fields{
				"nodeId":             key,
				"len(truth[key])":    len(positions),
				"len(computed[key])": len(computed[key]),
			}).Warn("different number of positions!")
			return false
		}

		sort.Sort(SortablePositions(positions))
		sort.Sort(SortablePositions(computed[key]))

		for i, position := range positions {
			if position != computed[key][i] {
				logrus.WithFields(logrus.Fields{
					"position truth":    position,
					"position computed": computed[key][i],
				}).Warn("positions are not equal!")
				return false
			}
		}
	}
	return true
}

func ReadTestcase(path string) *Testcase {

	file, errFile := os.Open(path)
	if errFile != nil {
		logrus.Fatal("Error reading testcase file", errFile)
		return nil
	}

	testCase := &Testcase{
		Sequences: make([]string, 0),
		Edges:     make([]Edge, 0),
		Positions: make(map[int][]Position),
	}

	scanner := bufio.NewScanner(file)

	annotationType := 0

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		if strings.Index(line, "sequences:") == 0 {
			annotationType = 1
			continue
		}
		if strings.Index(line, "edge list:") == 0 {
			annotationType = 2
			continue
		}
		if strings.Index(line, "positions:") == 0 {
			annotationType = 3
			continue
		}

		if annotationType == 1 {
			lineArray := strings.Split(line, " ")
			//sequenceIndex := strconv.Atoi(lineArray[0])
			sequence := lineArray[1]

			testCase.Sequences = append(testCase.Sequences, sequence)
		}
		if annotationType == 2 {
			lineArray := strings.Split(line, " ")
			fromId, _ := strconv.Atoi(lineArray[0])
			toId, _ := strconv.Atoi(lineArray[1])
			label := lineArray[2]

			edge := Edge{
				From:  fromId,
				To:    toId,
				Label: label,
			}
			testCase.Edges = append(testCase.Edges, edge)
		}
		if annotationType == 3 {
			lineArray := strings.Split(line, ": ")
			nodeId, _ := strconv.Atoi(lineArray[0])

			if len(lineArray) == 2 {
				positionsArray := strings.Split(lineArray[1], ",")

				testCase.Positions[nodeId] = make([]Position, len(positionsArray))

				for i, position := range positionsArray {
					positionArray := strings.Split(position, "-")
					index, _ := strconv.Atoi(positionArray[0])
					start, _ := strconv.Atoi(positionArray[1])

					testCase.Positions[nodeId][i] = Position{
						Index: index,
						Start: start,
					}
				}
			} else {
				testCase.Positions[nodeId] = make([]Position, 0)
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return nil
	}

	return testCase
}

func collectFilePaths(folderPath string) ([]string, error) {
	var filePaths []string

	err := filepath.WalkDir(folderPath, func(fp string, fi os.DirEntry, err error) error {
		if err != nil {
			return err
		}
		if !fi.IsDir() {
			filePaths = append(filePaths, fp)
		}
		return nil
	})

	if err != nil {
		return nil, err
	}

	return filePaths, nil
}

func GetTestcaseFiles() []string {
	folderPath := "../resources/test/"

	files, _ := collectFilePaths(folderPath)

	return files
}
