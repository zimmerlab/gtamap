package index

import (
	"bufio"
	"container/list"
	"os"
	"strconv"
	"strings"

	"github.com/KleinSamuel/gtamap/src/core/interval"
	"github.com/KleinSamuel/gtamap/src/formats/gtf"
	"github.com/sirupsen/logrus"
)

type RegionMask struct {
	TargetMasks []*interval.PriorityList
	ContigMasks map[string]*interval.PriorityList
	Priorities  map[string]int
}

func NewEmptyRegionMask(numTargets int) *RegionMask {
	return &RegionMask{
		TargetMasks: make([]*interval.PriorityList, numTargets),
		ContigMasks: make(map[string]*interval.PriorityList),
		Priorities:  make(map[string]int),
	}
}

func NewRegionMaskFromPaths(
	// prioritiesPath string,
	regionmaskFilePath string,
	// contigToTargetRegions map[string]*regionvector.RegionVector,
	targets []*gtf.GeneBasic,
) (*RegionMask, error) {
	// prioritiesFile, errP := os.Open(prioritiesPath)
	// if errP != nil {
	// 	return nil, errP
	// }

	regionmaskFile, err := os.Open(regionmaskFilePath)
	if err != nil {
		return nil, err
	}

	return NewRegionMask(
		regionmaskFile,
		// contigToTargetRegions,
		targets,
	)
}

// NewRegionMask creates a region mask from a BED file and a map of contig
// to target regions. These target regions are used to filter out regions in
// the BED. The priority file maps the name of a bed entry to an integer.
// The priority is used to construct the priority list for regions.
func NewRegionMask(
	regionmaskFile *os.File,
	// contigToTargetRegions map[string]*regionvector.RegionVector,
	targets []*gtf.GeneBasic,
) (*RegionMask, error) {
	mask := NewEmptyRegionMask(len(targets) * 2)

	if regionmaskFile == nil {

		for i := range mask.TargetMasks {
			mask.TargetMasks[i] = interval.NewList()
		}

		return mask, nil
	}

	scanner := bufio.NewScanner(regionmaskFile)

	// contigListTmp := make(map[string]*list.List)

	targetListTmp := make([]*list.List, len(targets)*2)
	for i := range targetListTmp {
		targetListTmp[i] = list.New()
	}

	for scanner.Scan() {
		line := scanner.Text()

		lineParts := strings.Split(line, "\t")

		if len(lineParts) < 6 {
			logrus.Error("Invalid BED line (needs 6 columns): ", line)
			continue
		}

		contig := lineParts[0]

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

		name := lineParts[3]

		priority, errPriority := strconv.Atoi(lineParts[4])
		if errPriority != nil {
			logrus.Error("Error parsing priority: ", errPriority)
			continue
		}

		maxMismatches, errMaxMismatches := strconv.Atoi(lineParts[5])
		if errMaxMismatches != nil {
			logrus.Error("Error parsing max mismatches: ", errMaxMismatches)
		}

		// find all target regions which overlap
		for i, target := range targets {

			if contig != target.Contig || start > int(target.EndGenomic) || end < int(target.StartGenomic) {
				continue
			}

			// add the 5->3 region

			startRelFw := start - int(target.StartGenomic)
			endRelFw := end - int(target.StartGenomic)

			targetListTmp[i].PushBack(
				&interval.PriorityListItem{
					Start:     startRelFw,
					End:       endRelFw,
					Name:      name,
					Priority:  priority,
					Threshold: maxMismatches,
				})

			// add the 3->5 region

			startRelRv := int(target.EndGenomic) - end
			endRelRv := int(target.EndGenomic) - start

			targetListTmp[i+1].PushBack(
				&interval.PriorityListItem{
					Start:     startRelRv,
					End:       endRelRv,
					Name:      name,
					Priority:  priority,
					Threshold: maxMismatches,
				})
		}

		// // the contig is not part of the target regions
		// if _, exists := contigToTargetRegions[contig]; !exists {
		// 	continue
		// }
		// // the region does not overlap any target region
		// if !contigToTargetRegions[contig].OverlapsAny(start, end) {
		// 	continue
		// }

		// if _, exists := contigListTmp[contig]; !exists {
		// 	contigListTmp[contig] = list.New()
		// }
		//
		// contigListTmp[contig].PushBack(&interval.PriorityListItem{
		// 	Start:     start,
		// 	End:       end,
		// 	Name:      name,
		// 	Priority:  priority,
		// 	Threshold: maxMismatches,
		// })
	}

	// convert the target lists to slice and sort
	for i, lst := range targetListTmp {
		regionSlice := make([]*interval.PriorityListItem, 0, lst.Len())
		for e := lst.Front(); e != nil; e = e.Next() {
			regionSlice = append(regionSlice, e.Value.(*interval.PriorityListItem))
		}
		mask.TargetMasks[i] = interval.NewListWithData(regionSlice)
	}

	// // convert tmp lists to slice and sort
	// for contig, lst := range contigListTmp {
	// 	regionSlice := make([]*interval.PriorityListItem, 0, lst.Len())
	// 	for e := lst.Front(); e != nil; e = e.Next() {
	// 		regionSlice = append(regionSlice, e.Value.(*interval.PriorityListItem))
	// 	}
	// 	mask.ContigMasks[contig] = interval.NewListWithData(regionSlice)
	// }

	return mask, nil
}
