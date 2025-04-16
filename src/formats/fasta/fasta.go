package fasta

import (
	"bufio"
	"github.com/sirupsen/logrus"
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

func ReadFastaIndexUsingPath(fastaIndexPath string) *Index {
	fastaIndexFile, err := os.Open(fastaIndexPath)
	if err != nil {
		logrus.Error(err)
		return nil
	}
	defer fastaIndexFile.Close()

	return ReadFastaIndex(fastaIndexFile)
}

func ReadFastaIndex(fastaIndexFile *os.File) *Index {

	fastaIndex := Index{
		Entries: make(map[string]*IndexEntry),
	}

	// buffered reader
	scanner := bufio.NewScanner(fastaIndexFile)

	for scanner.Scan() {
		line := strings.Split(scanner.Text(), "\t")

		reference := line[0]
		length64, errLrngth := strconv.ParseUint(line[1], 10, 32)
		offset64, errOffset := strconv.ParseUint(line[2], 10, 32)
		linebases64, errLinebases := strconv.ParseUint(line[3], 10, 16)
		linewidth64, errLinewidth := strconv.ParseUint(line[4], 10, 16)

		if errLrngth != nil {
			logrus.Error("Could not parse fasta entry length: ", errLrngth)
			return nil
		}
		if errOffset != nil {
			logrus.Error("Could not parse fasta entry offset: ", errOffset)
			return nil
		}
		if errLinebases != nil {
			logrus.Error("Could not parse fasta entry linebases: ", errLinebases)
			return nil
		}
		if errLinewidth != nil {
			logrus.Error("Could not parse fasta entry linewidth: ", errLinewidth)
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
