package index

import (
	"bufio"
	"container/list"
	"os"
	"strconv"
	"strings"

	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/interval"
	"github.com/KleinSamuel/gtamap/src/formats/bed"
	"github.com/sirupsen/logrus"
)

type RegionMask struct {
	ContigMasks map[string]*interval.PriorityList
	Priorities  map[string]int
}

func NewEmptyRegionMask() *RegionMask {
	return &RegionMask{
		ContigMasks: make(map[string]*interval.PriorityList),
		Priorities:  make(map[string]int),
	}
}

// func NewRegionMask(priorities map[string]int, bedFile *bed.File) *RegionMask {
//
// 	rm := &RegionMask{
// 		ContigMasks: make(map[string]*interval.PriorityList),
// 		Priorities:  priorities,
// 	}
//
// 	for contig, entries := range bedFile.NameMap {
// 		plist := interval.NewListWithData(priorities, entries)
// 		rm.ContigMasks[contig] = plist
// 	}
//
// 	return rm
// }

func NewRegionMaskFromPaths(
	prioritiesPath string,
	bedFilePath string,
	contigToTargetRegions map[string]*regionvector.RegionVector,
) (*RegionMask, error) {

	prioritiesFile, errP := os.Open(prioritiesPath)
	if errP != nil {
		return nil, errP
	}

	bedFile, errB := os.Open(bedFilePath)
	if errB != nil {
		return nil, errB
	}

	return NewRegionMask(prioritiesFile, bedFile, contigToTargetRegions)
}

func NewRegionMask(
	prioritiesFile *os.File,
	bedFile *os.File,
	contigToTargetRegions map[string]*regionvector.RegionVector,
) (*RegionMask, error) {

	priorities, errP := bed.ReadPriorities(prioritiesFile)
	if errP != nil {
		return nil, errP
	}

	mask := NewEmptyRegionMask()

	scanner := bufio.NewScanner(bedFile)

	contigListTmp := make(map[string]*list.List)

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

		contig := lineParts[0]

		// the contig is not part of the target regions
		if _, exists := contigToTargetRegions[contig]; !exists {
			continue
		}
		// the region does not overlap any target region
		if !contigToTargetRegions[contig].OverlapsAny(start, end) {
			continue
		}

		name := lineParts[3]

		if _, exists := contigListTmp[contig]; !exists {
			contigListTmp[contig] = list.New()
		}

		contigListTmp[contig].PushBack(&interval.PriorityListItem{
			Start:    start,
			End:      end,
			Name:     name,
			Priority: priorities[name],
		})
	}

	// convert tmp lists to slice and sort
	for contig, lst := range contigListTmp {
		regionSlice := make([]*interval.PriorityListItem, 0, lst.Len())
		for e := lst.Front(); e != nil; e = e.Next() {
			regionSlice = append(regionSlice, e.Value.(*interval.PriorityListItem))
		}
		mask.ContigMasks[contig] = interval.NewListWithData(priorities, regionSlice)
	}

	return mask, nil
}
