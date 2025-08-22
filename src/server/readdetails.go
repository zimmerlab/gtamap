package server

import (
	"fmt"
	"slices"
	"strings"

	"github.com/KleinSamuel/gtamap/src/dataloader"
)

type ReadDetailsDto struct {
	IsAcceptedR1 bool                `json:"isAcceptedR1"`
	IsAcceptedR2 bool                `json:"isAcceptedR2"`
	Locations    []*ReadLocationInfo `json:"locations"`
}

func (h *MappingDataHandler) GetReadDetailsTable(qname string) *ReadDetailsDto {

	details := &ReadDetailsDto{
		IsAcceptedR1: false,
		IsAcceptedR2: false,
		Locations:    make([]*ReadLocationInfo, 0),
	}

	qnameCluster := h.QnameCluster[qname]

	// numClusterR1 := len(qnameCluster.ClusterR1.SimilarRecords)
	// numClusterR2 := len(qnameCluster.ClusterR2.SimilarRecords)

	mappedBy := make(map[string]bool)

	for _, cR1 := range qnameCluster.ClusterR1.SimilarRecords {

		mappedByCluster := make(map[string]bool)
		for _, r := range cR1 {
			mappedBy[h.MapperInfos[r.MapperIndex].MapperName] = true
			mappedByCluster[h.MapperInfos[r.MapperIndex].MapperName] = true
		}
		mappedByClusterList := make([]string, 0, len(mappedByCluster))
		for mapperName := range mappedByCluster {
			mappedByClusterList = append(mappedByClusterList, mapperName)
		}
		slices.Sort(mappedByClusterList)

		locationInfo := ReadLocationInfo{
			Pair:                 "first",
			Contig:               cR1[0].Rname,
			Strand:               !cR1[0].Flag.IsReverseStrand(),
			Position:             cR1[0].Pos,
			Cigar:                cR1[0].CigarObj.StringUniform(),
			CigarDetailed:        cR1[0].CigarObj.String(),
			NumMismatches:        cR1[0].CigarObj.GetNumMismatches(),
			NumGaps:              cR1[0].CigarObj.GetNumGaps(),
			NumMappedBy:          len(mappedByCluster),
			MappedBy:             mappedByClusterList,
			ReadIndices:          make([]int, 0),
			ReadIndicesInMappers: make([]int, 0),
			IsAccepted:           cR1[0].IsAccepted,
			TargetRegionOverlap:  cR1[0].TargetRegionOverlap,
		}

		if cR1[0].IsAccepted {
			details.IsAcceptedR1 = true
		}

		for _, r := range cR1 {
			locationInfo.ReadIndices = append(locationInfo.ReadIndices, r.Index)
			locationInfo.ReadIndicesInMappers = append(locationInfo.ReadIndicesInMappers, r.IndexInSam)
		}

		details.Locations = append(details.Locations, &locationInfo)
	}

	for _, cR2 := range qnameCluster.ClusterR2.SimilarRecords {

		mappedByCluster := make(map[string]bool)
		for _, r := range cR2 {
			mappedBy[h.MapperInfos[r.MapperIndex].MapperName] = true
			mappedByCluster[h.MapperInfos[r.MapperIndex].MapperName] = true
		}
		mappedByClusterList := make([]string, 0, len(mappedByCluster))
		for mapperName := range mappedByCluster {
			mappedByClusterList = append(mappedByClusterList, mapperName)
		}
		slices.Sort(mappedByClusterList)

		locationInfo := ReadLocationInfo{
			Pair:                 "second",
			Contig:               cR2[0].Rname,
			Strand:               !cR2[0].Flag.IsReverseStrand(),
			Position:             cR2[0].Pos,
			Cigar:                cR2[0].CigarObj.StringUniform(),
			CigarDetailed:        cR2[0].CigarObj.String(),
			NumMismatches:        cR2[0].CigarObj.GetNumMismatches(),
			NumGaps:              cR2[0].CigarObj.GetNumGaps(),
			NumMappedBy:          len(mappedByCluster),
			MappedBy:             mappedByClusterList,
			ReadIndices:          make([]int, 0),
			ReadIndicesInMappers: make([]int, 0),
			IsAccepted:           cR2[0].IsAccepted,
			TargetRegionOverlap:  cR2[0].TargetRegionOverlap,
		}

		if cR2[0].IsAccepted {
			details.IsAcceptedR2 = true
		}

		for _, r := range cR2 {
			locationInfo.ReadIndices = append(locationInfo.ReadIndices, r.Index)
			locationInfo.ReadIndicesInMappers = append(locationInfo.ReadIndicesInMappers, r.IndexInSam)
		}

		details.Locations = append(details.Locations, &locationInfo)
	}

	return details
}

type ReadDetailsViewerData struct {
	// Intervals    []*Interval                     `json:"intervals"`
	// IntervalsMap map[*Interval][]*EnhancedRecord `json:"intervalsMap"`
	ViewerConfigs []map[string]any `json:"viewerConfigs"`
}

func (h *MappingDataHandler) GetReadDetailsViewerData(qname string) *ReadDetailsViewerData {

	paddingBases := 500

	intervals := make([]*Interval, 0)
	intervalsMap := make(map[*Interval][]*EnhancedRecord)

	qclust := h.QnameCluster[qname]

	for _, r := range append(qclust.ClusterR1.SimilarRecords, qclust.ClusterR2.SimilarRecords...) {

		interval := GetGenomicIntervalFromRecord(r[0])
		interval.Start -= paddingBases
		interval.End += paddingBases

		overlapping := make([]int, 0)
		for i, e := range intervals {
			fmt.Println(i, e.Start, e.End)
			if interval.Start <= e.End && interval.End >= e.Start {
				overlapping = append(overlapping, i)
			}
		}

		if len(overlapping) > 0 {

			newIntervals := make([]*Interval, 0)
			newIntervalsMap := make(map[*Interval][]*EnhancedRecord)

			minStart := interval.Start
			maxEnd := interval.End
			for _, idx := range overlapping {
				if intervals[idx].Start < minStart {
					minStart = intervals[idx].Start
				}
				if intervals[idx].End > maxEnd {
					maxEnd = intervals[idx].End
				}
			}

			newInterval := &Interval{
				Contig: interval.Contig,
				Start:  minStart,
				End:    maxEnd,
			}

			newIntervals = append(newIntervals, newInterval)
			newIntervalsMap[newInterval] = make([]*EnhancedRecord, 0)
			newIntervalsMap[newInterval] = append(newIntervalsMap[newInterval], r...)

			for i, e := range intervals {
				if slices.Contains(overlapping, i) {
					newIntervalsMap[newInterval] = append(newIntervalsMap[newInterval], intervalsMap[e]...)
				} else {
					newIntervals = append(newIntervals, e)
					newIntervalsMap[e] = intervalsMap[e]
				}
			}

			intervals = newIntervals
			intervalsMap = newIntervalsMap

		} else {

			intervals = append(intervals, interval)
			intervalsMap[interval] = make([]*EnhancedRecord, 0)
			intervalsMap[interval] = append(intervalsMap[interval], r...)
		}
	}

	viewerConfigs := make([]map[string]any, 0)

	for xi, i := range intervals {

		c := intervalsMap[i]

		// identifier for this region within the read
		id := fmt.Sprintf("%s_%d", qname, xi)

		header := fmt.Sprintf(">%s\n", i.Contig)
		seq := dataloader.ExtractSequenceAsStringFromFasta(h.FastaFile, h.FastaIndex, i.Contig, uint32(i.Start), uint32(i.End))

		fastaIndex := fmt.Sprintf("%s\t%d\t%d\t%d\t%d", i.Contig, len(seq), len(header), len(seq), len(seq)+1)

		details := &DetailsViewerData{
			Interval:         i,
			FastaSequence:    header + seq + "\n",
			FastaIndexString: fastaIndex,
			Records:          c,
		}

		h.DetailsViewerData[id] = details

		config := h.GetViewerConfig(id, i.Contig, i.Start, i.End)

		viewerConfigs = append(viewerConfigs, config)
	}

	return &ReadDetailsViewerData{
		ViewerConfigs: viewerConfigs,
	}
}

func (h *MappingDataHandler) GetViewerConfig(id string, contig string, start int, end int) map[string]any {

	genomeConfig := IgvGenomeConfig{
		Id:       "chr " + contig + " (relative)",
		Label:    "BASE_GENOME_LABEL",
		FastaUrl: "http://localhost:8000/api/details/viewer/fasta?id=" + id,
		IndexUrl: "http://localhost:8000/api/details/viewer/fastaIndex?id=" + id,
	}

	track := IgvTrackConfig{
		Name:        "Read Mappings",
		Format:      "sam",
		DisplayMode: "EXPANDED",
		Url:         "http://localhost:8000/api/details/viewer/bam?id=" + id,
		IndexUrl:    "http://localhost:8000/api/details/viewer/bamIndex?id=" + id,
		Type:        "alignment",
		Height:      100 + 20*(len(h.DetailsViewerData[id].Records)),
		MaxHeight:   1000,
		MaxRows:     200,
		ColorBy:     "none",
	}

	tracks := make([]IgvTrackConfig, 0)
	tracks = append(tracks, track)

	// convert global location to local location
	targetRegionStrLocal := fmt.Sprintf("%s:%d-%d", contig, 1, end-start+1)

	response := map[string]any{
		"genomeConfig": genomeConfig,
		"tracks":       tracks,
		"location":     targetRegionStrLocal,
	}

	return response
}

func (h *MappingDataHandler) GetReadGroupSam(groupId string) string {

	details := h.DetailsViewerData[groupId]

	fmt.Println(details.Interval.Contig, details.Interval.Start, details.Interval.End)

	var builder strings.Builder

	builder.WriteString("@HD\tVN:1.6\tSO:unsorted\n")
	builder.WriteString(fmt.Sprintf("@SQ\tSN:%s\tLN:%d\n", details.Interval.Contig, details.Interval.End-details.Interval.Start+1))

	for _, r := range details.Records {

		fmt.Println(r.Pos)

		builder.WriteString(r.StringWithMapperInfo(details.Interval.Start, h.MapperInfos[r.MapperIndex].MapperName, r.Index))
		builder.WriteString("\n")
	}

	return builder.String()
}
