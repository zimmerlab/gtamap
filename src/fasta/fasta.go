package fasta

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

type Index struct {
	Entries map[string]*IndexEntry
}

type IndexEntry struct {
	Reference string
	Length    uint32
	Offset    int64
	Linebases uint16
	Linewidth uint16
}

func ReadFastaIndex(filePath string) *Index {

	fastaIndex := Index{
		Entries: make(map[string]*IndexEntry),
	}

	file, err := os.Open(filePath)
	if err != nil {
		return nil
	}
	defer file.Close()

	// buffered reader
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := strings.Split(scanner.Text(), "\t")

		reference := line[0]
		length64, errLrngth := strconv.ParseUint(line[1], 10, 16)
		offset64, errOffset := strconv.ParseUint(line[2], 10, 16)
		linebases64, errLinebases := strconv.ParseUint(line[3], 10, 16)
		linewidth64, errLinewidth := strconv.ParseUint(line[4], 10, 16)

		if errLrngth != nil || errOffset != nil || errLinebases != nil || errLinewidth != nil {
			fmt.Println("Error while reading fasta index")
			return nil
		}

		entry := IndexEntry{
			Reference: reference,
			Length:    uint32(length64),
			Offset:    int64(offset64),
			Linebases: uint16(linebases64),
			Linewidth: uint16(linewidth64),
		}

		fastaIndex.Entries[reference] = &entry
	}

	if err := scanner.Err(); err != nil {
		return nil
	}

	return &fastaIndex
}
