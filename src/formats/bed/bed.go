package bed

import (
	"bufio"
	"os"
	"strconv"
	"strings"

	"github.com/sirupsen/logrus"
)

type Entry struct {
	Contig string
	Start  int
	End    int
	Name   string
}

type File struct {
	Entries   []*Entry
	ContigMap map[string][]*Entry
	NameMap   map[string][]*Entry
}

func NewFromPath(filePath string) (*File, error) {
	file, err := os.Open(filePath)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	return New(file), nil
}

func New(file *os.File) *File {

	entries := make([]*Entry, 0)
	nameMap := make(map[string][]*Entry)
	contigMap := make(map[string][]*Entry)

	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := scanner.Text()

		lineParts := strings.Split(line, "\t")

		if len(lineParts) < 4 {
			logrus.Error("Invalid BED line: ", line)
			continue
		}

		start, errStart := strconv.Atoi(lineParts[1])
		if errStart != nil {
			logrus.Error("Error parsing start position: ", errStart)
			continue
		}
		end, errEnd := strconv.Atoi(lineParts[2])
		if errEnd != nil {
			logrus.Error("Error parsing end position: ", errEnd)
			continue
		}

		entry := &Entry{
			Contig: lineParts[0],
			Start:  start,
			End:    end,
			Name:   lineParts[3],
		}

		entries = append(entries, entry)

		if _, exists := nameMap[entry.Name]; !exists {
			nameMap[entry.Name] = make([]*Entry, 0)
		}
		nameMap[entry.Name] = append(nameMap[entry.Name], entry)

		if _, exists := contigMap[entry.Contig]; !exists {
			contigMap[entry.Contig] = make([]*Entry, 0)
		}
		contigMap[entry.Contig] = append(contigMap[entry.Contig], entry)
	}

	return &File{
		Entries:   entries,
		ContigMap: contigMap,
		NameMap:   nameMap,
	}
}
