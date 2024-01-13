package test

import (
	"bufio"
	"fmt"
	"github.com/sirupsen/logrus"
	"os"
	"strings"
)

type Edge struct {
	Id    int
	From  int
	To    int
	Label string
}

type Position struct {
	Index int
	Start int
}

type Testcase struct {
	Sequences []string
	Nodes     []int
	Edges     map[int]Edge
	Positions map[int]Position
}

func ReadTestcase(path string) *Testcase {

	file, errFile := os.Open(path)
	if errFile != nil {
		logrus.Fatal("Error reading testcase file", errFile)
		return nil
	}

	testCase := &Testcase{
		Sequences: make([]string, 0),
		Nodes:     make([]int, 0),
		Edges:     make(map[int]Edge),
		Positions: make(map[int]Position),
	}

	scanner := bufio.NewScanner(file)

	annotationType := 0

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		fmt.Print(line)

		if strings.Index(line, "sequences:") == 0 {
			fmt.Println("got sequences line")
			annotationType = 1
		}
		if strings.Index(line, "edge list:") == 0 {
			fmt.Println("got edge list line")
			annotationType = 2
		}
		if strings.Index(line, "positions:") == 0 {
			fmt.Println("got positions line")
			annotationType = 3
		}

		if annotationType == 1 {
			lineArray := strings.Split(line, " ")
			//sequenceIndex := strconv.Atoi(lineArray[0])
			sequence := lineArray[1]

			testCase.Sequences = append(testCase.Sequences, sequence)
		}

	}

	if err := scanner.Err(); err != nil {
		return nil
	}

	return testCase
}
