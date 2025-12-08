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

func Abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

func ScoreSpliceSites(
	donorFirstBase byte,
	donorSecondBase byte,
	acceptorFirstBase byte,
	acceptorSecondBase byte,
	isForwardStrand bool,
) (int, bool) {
	if isForwardStrand {
		if donorFirstBase == byte('G') && donorSecondBase == byte('T') &&
			acceptorFirstBase == byte('A') && acceptorSecondBase == byte('G') {
			// canonical splice site GT/AG
			return 0, true
		} else if donorFirstBase == byte('G') && donorSecondBase == byte('C') &&
			acceptorFirstBase == byte('A') && acceptorSecondBase == byte('G') {
			// non-canonical splice site GC/AG
			return 1, true
		} else if donorFirstBase == byte('A') && donorSecondBase == byte('T') &&
			acceptorFirstBase == byte('A') && acceptorSecondBase == byte('C') {
			// non-canonical splice site AT/AC
			return 1, true
		}
	} else {
		if donorFirstBase == byte('C') && donorSecondBase == byte('T') &&
			acceptorFirstBase == byte('A') && acceptorSecondBase == byte('C') {
			// canonical splice site GT/AG on rev strand
			return 0, true
		} else if donorFirstBase == byte('C') && donorSecondBase == byte('T') &&
			acceptorFirstBase == byte('G') && acceptorSecondBase == byte('C') {
			// non-canonical splice site GC/AG on rev strand
			return 1, true
		} else if donorFirstBase == byte('G') && donorSecondBase == byte('T') &&
			acceptorFirstBase == byte('A') && acceptorSecondBase == byte('T') {
			// non-canonical splice site AT/AC on rev strand
			return 1, true
		}
	}

	// all other non-canonical splice sites
	return 2, false
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

func ReverseBytesInplace(s []byte) {
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
}

func ReverseBytes(s []byte) []byte {
	newBytes := make([]byte, len(s))
	for i := range s {
		newBytes[i] = s[len(s)-1-i]
	}
	return newBytes
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

func ArrayStringToString(arr []string, sep string) string {
	if len(arr) == 0 {
		return ""
	}
	var builder strings.Builder
	builder.WriteString(arr[0])
	for _, item := range arr[1:] {
		builder.WriteString(sep)
		builder.WriteString(item)
	}
	return builder.String()
}

func ArrayIntToString(arr []int, sep string) string {
	if len(arr) == 0 {
		return ""
	}
	var builder strings.Builder
	builder.WriteString(fmt.Sprintf("%d", arr[0]))
	for _, item := range arr[1:] {
		builder.WriteString(sep)
		builder.WriteString(fmt.Sprintf("%d", item))
	}
	return builder.String()
}
