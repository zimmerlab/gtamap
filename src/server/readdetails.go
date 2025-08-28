package server

import (
	"fmt"
	"slices"
	"strconv"
	"strings"

	"github.com/KleinSamuel/gtamap/src/dataloader"
)

type ReadDetailsDto struct {
	IsDiscarded  bool                `json:"isDiscarded"`
	IsAcceptedR1 bool                `json:"isAcceptedR1"`
	IsAcceptedR2 bool                `json:"isAcceptedR2"`
	Locations    []*ReadLocationInfo `json:"locations"`
}

func (h *MappingDataHandler) GetReadDetailsTable(qname string) *ReadDetailsDto {

	details := &ReadDetailsDto{
		IsDiscarded:  false,
		IsAcceptedR1: false,
		IsAcceptedR2: false,
		Locations:    make([]*ReadLocationInfo, 0),
	}

	qnameCluster := h.QnameCluster[qname]

	details.IsDiscarded = qnameCluster.IsDiscarded

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

type ReadDetailsData struct {
	Confidence   int            `json:"confidence"`
	ReadGroups   []*ReadGroup   `json:"readGroups"`
	R1Mapped     []*MapperState `json:"r1Mapped"`
	R2Mapped     []*MapperState `json:"r2Mapped"`
	IsDiscarded  bool           `json:"isDiscarded"`
	IsAcceptedR1 bool           `json:"isAcceptedR1"`
	IsAcceptedR2 bool           `json:"isAcceptedR2"`
}

type MapperState struct {
	MapperName string `json:"mapperName"`
	Mapped     bool   `json:"mapped"`
}

type ReadGroup struct {
	Interval     *Interval           `json:"interval"`
	Name         string              `json:"name"`
	ViewerConfig map[string]any      `json:"viewerConfig"`
	Locations    []*ReadLocationInfo `json:"locations"`
}

func (h *MappingDataHandler) GetReadDetailsData(qname string) *ReadDetailsData {

	paddingBases := 500

	intervals := make([]*Interval, 0)
	intervalsMap := make(map[*Interval][]*EnhancedRecord)

	qclust := h.QnameCluster[qname]

	r1Mappers := make(map[string]bool, 0)
	for _, r := range qclust.ClusterR1.Records {
		mapperName := h.MapperInfos[r.MapperIndex].MapperName
		r1Mappers[mapperName] = true
	}

	r2Mappers := make(map[string]bool, 0)
	for _, r := range qclust.ClusterR2.Records {
		mapperName := h.MapperInfos[r.MapperIndex].MapperName
		r2Mappers[mapperName] = true
	}

	r1Mapped := make([]*MapperState, len(h.MapperInfos))
	r2Mapped := make([]*MapperState, len(h.MapperInfos))
	for i, mapperInfo := range h.MapperInfos {
		r1Mapped[i] = &MapperState{
			MapperName: mapperInfo.MapperName,
			Mapped:     r1Mappers[mapperInfo.MapperName],
		}
		r2Mapped[i] = &MapperState{
			MapperName: mapperInfo.MapperName,
			Mapped:     r2Mappers[mapperInfo.MapperName],
		}
	}

	for _, r := range append(qclust.ClusterR1.SimilarRecords, qclust.ClusterR2.SimilarRecords...) {

		interval := GetGenomicIntervalFromRecord(r[0])
		interval.Start -= paddingBases
		interval.End += paddingBases

		overlapping := make([]int, 0)
		for i, e := range intervals {
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

	readGroups := make([]*ReadGroup, 0)

	for xi, i := range intervals {

		c := intervalsMap[i]

		// identifier for this region within the read
		id := fmt.Sprintf("%s_%d", qname, xi)

		header := fmt.Sprintf(">%s\n", i.Contig)
		seq := dataloader.ExtractSequenceAsStringFromFasta(h.FastaFile, h.FastaIndex, i.Contig, uint32(i.Start), uint32(i.End))

		fastaIndex := fmt.Sprintf("%s\t%d\t%d\t%d\t%d", i.Contig, len(seq), len(header), len(seq), len(seq)+1)

		mapperNames := make([]string, 0)
		recordsByMapper := make([][]*EnhancedRecord, 0)

		// keeps the order of the mapper names
		for _, mapperInfo := range h.MapperInfos {

			mapperIndex := len(mapperNames)
			found := false

			for _, r := range c {
				mapperName := h.MapperInfos[r.MapperIndex].MapperName

				if mapperName != mapperInfo.MapperName {
					continue
				}

				if !found {
					mapperNames = append(mapperNames, mapperName)
					recordsByMapper = append(recordsByMapper, make([]*EnhancedRecord, 0))
				}
				found = true

				recordsByMapper[mapperIndex] = append(recordsByMapper[mapperIndex], r)
			}
		}

		details := &DetailsViewerData{
			Interval:         i,
			FastaSequence:    header + seq + "\n",
			FastaIndexString: fastaIndex,
			Records:          c,
			MapperNames:      mapperNames,
			RecordsByMapper:  recordsByMapper,
		}

		h.DetailsViewerData[id] = details

		config := h.GetViewerConfig(id, i.Contig, i.Start, i.End)

		// read groups
		locations := make([]*ReadLocationInfo, 0)

		for _, r := range c {

			pair := "first"
			if !r.Flag.IsFirstInPair() {
				pair = "second"
			}

			locationInfo := &ReadLocationInfo{
				Pair:                 pair,
				Contig:               r.Rname,
				Strand:               !r.Flag.IsReverseStrand(),
				Position:             r.Pos,
				Cigar:                r.CigarObj.StringUniform(),
				CigarDetailed:        r.CigarObj.String(),
				NumMismatches:        r.CigarObj.GetNumMismatches(),
				NumGaps:              r.CigarObj.GetNumGaps(),
				NumMappedBy:          1,
				MappedBy:             []string{h.MapperInfos[r.MapperIndex].MapperName},
				ReadIndices:          []int{r.Index},
				ReadIndicesInMappers: []int{r.MapperIndex},
				IsAccepted:           r.IsAccepted,
				TargetRegionOverlap:  r.TargetRegionOverlap,
			}

			locations = append(locations, locationInfo)
		}

		group := &ReadGroup{
			Name:         id,
			Interval:     i,
			ViewerConfig: config,
			Locations:    locations,
		}

		readGroups = append(readGroups, group)
	}

	// acceptance state

	return &ReadDetailsData{
		IsDiscarded:  qclust.IsDiscarded,
		IsAcceptedR1: qclust.ClusterR1.AcceptedRecord != nil,
		IsAcceptedR2: qclust.ClusterR2.AcceptedRecord != nil,
		Confidence:   qclust.ConfidenceLevel,
		ReadGroups:   readGroups,
		R1Mapped:     r1Mapped,
		R2Mapped:     r2Mapped,
	}
}

func (h *MappingDataHandler) GetViewerConfig(id string, contig string, start int, end int) map[string]any {

	genomeConfig := IgvGenomeConfig{
		Id:       "chr " + contig + " (relative)",
		Label:    "BASE_GENOME_LABEL",
		FastaUrl: "http://localhost:8000/api/details/viewer/fasta?id=" + id,
		IndexUrl: "http://localhost:8000/api/details/viewer/fastaIndex?id=" + id,
	}

	tracks := make([]IgvTrackConfig, 0)

	for i, mapperName := range h.DetailsViewerData[id].MapperNames {

		track := IgvTrackConfig{
			Name:              mapperName,
			Format:            "sam",
			DisplayMode:       "EXPANDED",
			Url:               "http://localhost:8000/api/details/viewer/bam?id=" + id + "&i=" + strconv.Itoa(i),
			IndexUrl:          "http://localhost:8000/api/details/viewer/bamIndex?id=" + id + "&i=" + strconv.Itoa(i),
			Type:              "alignment",
			Height:            20 + 20*(len(h.DetailsViewerData[id].RecordsByMapper[i])),
			MaxHeight:         1000,
			MaxRows:           200,
			ColorBy:           "none",
			ShowCoverage:      false,
			ShowSoftClips:     true,
			ShowInsertionText: true,
		}

		tracks = append(tracks, track)
	}

	// convert global location to local location
	targetRegionStrLocal := fmt.Sprintf("%s:%d-%d", contig, 1, end-start+1)

	response := map[string]any{
		"genomeConfig": genomeConfig,
		"tracks":       tracks,
		"location":     targetRegionStrLocal,
	}

	return response
}

func (h *MappingDataHandler) GetReadGroupSam(groupId string, mapperIndex int) string {

	details := h.DetailsViewerData[groupId]

	var builder strings.Builder

	builder.WriteString("@HD\tVN:1.6\tSO:unsorted\n")
	builder.WriteString(fmt.Sprintf("@SQ\tSN:%s\tLN:%d\n", details.Interval.Contig, details.Interval.End-details.Interval.Start+1))

	for _, r := range details.RecordsByMapper[mapperIndex] {
		builder.WriteString(r.StringWithMapperInfo(details.Interval.Start, h.MapperInfos[r.MapperIndex].MapperName, r.Index))
		builder.WriteString("\n")
	}

	return builder.String()
}
