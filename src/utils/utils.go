package utils

import (
	"fmt"
	"strings"
	"time"
)

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

var complementMapBytes = map[byte]byte{
	'A': 'T',
	'T': 'A',
	'G': 'C',
	'C': 'G',
	'N': 'N',
}

func ReverseComplementDnaBytes(dna []byte) ([]byte, error) {
	complement, err := ComplementDnaBytes(dna)
	if err != nil {
		return nil, err
	}
	return ReverseBytes(complement), nil
}

func ComplementDnaBytes(dna []byte) ([]byte, error) {
	complement := make([]byte, len(dna))

	for i, nucleotide := range dna {

		if _, ok := complementMapBytes[nucleotide]; !ok {
			return nil, fmt.Errorf("unknown nucleotide: %s", string(nucleotide))
		}

		complement[i] = complementMapBytes[nucleotide]
	}

	return complement, nil
}

func ReverseBytes(s []byte) []byte {
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
	return s
}

func FormatDuration(d time.Duration) string {
	if d.Hours() < 1 {
		if d.Minutes() < 1 {
			if d.Seconds() < 1 {
				return fmt.Sprintf("%dms", int(d.Milliseconds()))
			} else {
				return fmt.Sprintf("%ds", int(d.Seconds()))
			}
		} else {
			return fmt.Sprintf("%dm %ds", int(d.Minutes()), int(d.Seconds())%60)
		}
	} else {
		return fmt.Sprintf("%dh %dm %ds", int(d.Hours()), int(d.Minutes())%60, int(d.Seconds())%60)
	}
}
