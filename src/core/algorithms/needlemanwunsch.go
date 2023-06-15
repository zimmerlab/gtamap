package algorithms

import (
	"fmt"
)

type Cell struct {
	score   int
	pointer string
}

const (
	GapPenalty      = -2
	MismatchPenalty = -1
	MatchReward     = 1
)

func NeedlemanWunsch(reference, read string) (int, []rune, string, string) {
	rows := len(reference) + 1
	cols := len(read) + 1

	// Initialize the matrix
	matrix := make([][]Cell, rows)
	for i := 0; i < rows; i++ {
		matrix[i] = make([]Cell, cols)
	}

	// Initialize the first row and column with gap penalties
	for i := 0; i < rows; i++ {
		matrix[i][0].score = i * GapPenalty
		matrix[i][0].pointer = "U"
	}
	for j := 0; j < cols; j++ {
		matrix[0][j].score = j * GapPenalty
		matrix[0][j].pointer = "L"
	}

	// Fill the matrix
	for i := 1; i < rows; i++ {
		for j := 1; j < cols; j++ {
			scoreMatch := matrix[i-1][j-1].score + matchScore(reference[i-1], read[j-1])
			scoreDelete := matrix[i-1][j].score + GapPenalty
			scoreInsert := matrix[i][j-1].score + GapPenalty

			// Choose the maximum score and set the corresponding pointer
			maxScore := max(scoreMatch, scoreDelete, scoreInsert)
			switch maxScore {
			case scoreMatch:
				matrix[i][j].score = scoreMatch
				matrix[i][j].pointer = "D"
			case scoreDelete:
				matrix[i][j].score = scoreDelete
				matrix[i][j].pointer = "U"
			case scoreInsert:
				matrix[i][j].score = scoreInsert
				matrix[i][j].pointer = "L"
			}
		}
	}

	// Traceback to find the aligned sequences
	align1, align2 := "", ""

	cigarList := make([]rune, 0)

	i, j := rows-1, cols-1
	for i > 0 || j > 0 {
		switch matrix[i][j].pointer {
		case "D":
			align1 = string(reference[i-1]) + align1
			align2 = string(read[j-1]) + align2
			i--
			j--
			if reference[i] == read[j] {
				cigarList = append(cigarList, 'M')
			} else {
				cigarList = append(cigarList, 'X')
			}
		case "U":
			align1 = string(reference[i-1]) + align1
			align2 = "-" + align2
			i--
			cigarList = append(cigarList, 'D')
		case "L":
			align1 = "-" + align1
			align2 = string(read[j-1]) + align2
			j--
			cigarList = append(cigarList, 'I')
		}
	}

	return matrix[rows-1][cols-1].score, cigarList, align1, align2
}

func matchScore(c1, c2 byte) int {
	if c1 == c2 {
		return MatchReward
	}
	return MismatchPenalty
}

func max(a, b, c int) int {
	if a >= b && a >= c {
		return a
	} else if b >= a && b >= c {
		return b
	}
	return c
}

func Main() {
	seq1 := "AGTACGCA"
	seq2 := "TATGC"

	score, cigarList, align1, align2 := NeedlemanWunsch(seq1, seq2)

	fmt.Println("Alignment score:", score)
	fmt.Println("CIGAR:", string(cigarList))
	fmt.Println("Aligned sequence 1:", align1)
	fmt.Println("Aligned sequence 2:", align2)
}
