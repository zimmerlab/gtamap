package utils

import "strings"

func Arange(start, end, size int) []int {
	step := (end - start) / (size - 1)
	values := make([]int, size)

	for i := 0; i < size; i++ {
		values[i] = start + i*step
	}

	return values
}

var complementMap = map[string]string{
	"A": "T",
	"T": "A",
	"G": "C",
	"C": "G",
}

func ReverseComplementDNA(dna string) string {
	return ReverseString(ComplementDNA(dna))
}

func ComplementDNA(dna string) string {
	var complement strings.Builder

	for _, nucleotide := range dna {
		complement.WriteString(complementMap[string(nucleotide)])
	}

	return complement.String()
}

func ReverseString(s string) string {
	runes := []rune(s)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	return string(runes)
}
