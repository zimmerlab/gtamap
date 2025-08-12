package cigar

import (
	"fmt"
	"strconv"
	"strings"
)

// cigar elements
var ElementTypes = map[byte]bool{
	'M': true, // match or mismatch
	'X': true, // mismatch
	'=': true, // match
	'I': true, // insertion
	'D': true, // deletion
	'N': true, // skipped region from the reference
	'S': true, // soft clipping
	'H': true, // hard clipping
	'P': true, // padding
}

type Element struct {
	Length int    // length of the cigar element
	Type   string // type of the cigar element (M, I, D, etc.)
}

func (ce *Element) IsMatch() bool {
	return ce.Type == "M" || ce.Type == "X" || ce.Type == "="
}

func (ce *Element) ConsumesRead() bool {
	return ce.Type == "M" || ce.Type == "I" || ce.Type == "S" || ce.Type == "=" || ce.Type == "X"
}

func (ce *Element) ConsumesReference() bool {
	return ce.Type == "M" || ce.Type == "D" || ce.Type == "N" || ce.Type == "=" || ce.Type == "X"
}

type Object struct {
	Elements []Element // list of cigar elements
}

func NewCigarFromString(cigarString string) (*Object, error) {

	l := 0

	elements := make([]Element, 0)

	for i := 0; i < len(cigarString); i++ {

		// check if char at i is a cigar element char
		if _, ok := ElementTypes[cigarString[i]]; !ok {
			continue
		}

		num, err := strconv.ParseInt(string(cigarString[l:i]), 10, 64)
		if err != nil {
			return nil, fmt.Errorf("error parsing cigar string: %w", err)
		}

		elements = append(elements, Element{
			Length: int(num),
			Type:   string(cigarString[i]),
		})

		l = i + 1
	}

	cigar := &Object{
		Elements: elements,
	}

	return cigar, nil
}

func (c *Object) String() string {
	var sb strings.Builder
	for _, elem := range c.Elements {
		sb.WriteString(fmt.Sprintf("%d%s", elem.Length, elem.Type))
	}
	return sb.String()
}

func (c *Object) StringUniform() string {
	var sb strings.Builder
	currentLen := 0
	for _, elem := range c.Elements {
		if elem.Type == "M" || elem.Type == "X" || elem.Type == "=" {
			currentLen += elem.Length
			continue
		}
		if currentLen > 0 {
			sb.WriteString(fmt.Sprintf("%d%s", currentLen, "M"))
			currentLen = 0
		}
		sb.WriteString(fmt.Sprintf("%d%s", elem.Length, elem.Type))
	}
	if currentLen > 0 {
		sb.WriteString(fmt.Sprintf("%d%s", currentLen, "M"))
	}

	return sb.String()
}

func (c *Object) GetLengthOfMatched() int {
	length := 0

	for _, elem := range c.Elements {
		if elem.Type == "M" || elem.Type == "X" || elem.Type == "=" {
			length += elem.Length
		}
	}

	return length
}

func (c *Object) GetNumMismatches() int {
	num := 0
	for _, elem := range c.Elements {
		if elem.Type == "X" {
			num += elem.Length
		}
	}
	return num
}

func (c *Object) GetGapLength() int {
	num := 0
	for _, elem := range c.Elements {
		if elem.Type == "D" || elem.Type == "N" {
			num += elem.Length
		}
	}
	return num
}

func (c *Object) GetNumGaps() int {
	num := 0
	for _, elem := range c.Elements {
		// deletions (D) and skipped regions (N) should not be consecutive
		if elem.Type == "D" || elem.Type == "N" {
			num++
		}
	}
	return num
}
