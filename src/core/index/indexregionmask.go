package index

import (
	"bufio"
	"container/list"
	"os"
	"strconv"
	"strings"

	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/interval"
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

func NewRegionMaskFromPaths(
	// prioritiesPath string,
	regionmaskFilePath string,
	contigToTargetRegions map[string]*regionvector.RegionVector,
) (*RegionMask, error) {

	// prioritiesFile, errP := os.Open(prioritiesPath)
	// if errP != nil {
	// 	return nil, errP
	// }

	regionmaskFile, err := os.Open(regionmaskFilePath)
	if err != nil {
		return nil, err
	}

	return NewRegionMask(regionmaskFile, contigToTargetRegions)
}

// NewRegionMask creates a region mask from a BED file and a map of contig
// to target regions. These target regions are used to filter out regions in
// the BED. The priority file maps the name of a bed entry to an integer.
// The priority is used to construct the priority list for regions.
func NewRegionMask(
	// prioritiesFile *os.File,
	regionmaskFile *os.File,
	contigToTargetRegions map[string]*regionvector.RegionVector,
) (*RegionMask, error) {

	// priorities, errP := bed.ReadPriorities(prioritiesFile)
	// if errP != nil {
	// 	return nil, errP
	// }

	mask := NewEmptyRegionMask()

	scanner := bufio.NewScanner(regionmaskFile)

	contigListTmp := make(map[string]*list.List)

	for scanner.Scan() {
		line := scanner.Text()

		lineParts := strings.Split(line, "\t")

		if len(lineParts) < 6 {
			logrus.Error("Invalid BED line (needs 6 columns): ", line)
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

		priority, errPriority := strconv.Atoi(lineParts[4])
		if errPriority != nil {
			logrus.Error("Error parsing priority: ", errPriority)
			continue
		}

		maxMismatches, errMaxMismatches := strconv.Atoi(lineParts[5])
		if errMaxMismatches != nil {
			logrus.Error("Error parsing max mismatches: ", errMaxMismatches)
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
			Start:     start,
			End:       end,
			Name:      name,
			Priority:  priority,
			Threshold: maxMismatches,
		})
	}

	// convert tmp lists to slice and sort
	for contig, lst := range contigListTmp {
		regionSlice := make([]*interval.PriorityListItem, 0, lst.Len())
		for e := lst.Front(); e != nil; e = e.Next() {
			regionSlice = append(regionSlice, e.Value.(*interval.PriorityListItem))
		}
		mask.ContigMasks[contig] = interval.NewListWithData(regionSlice)
	}

	return mask, nil
}
