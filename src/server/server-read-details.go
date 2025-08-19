package server

import (
	"encoding/json"
	"net/http"
	"slices"
)

func getQnameParam(r *http.Request) string {
	qname := r.URL.Query().Get("qname")
	if qname == "" {
		return ""
	}
	return qname
}

func (s *Server) readDetailsTable(w http.ResponseWriter, r *http.Request) {

	qname := getQnameParam(r)
	if qname == "" {
		w.WriteHeader(http.StatusBadRequest)
		w.Write([]byte("qname parameter is required"))
		return
	}

	info := make([]*ReadOverviewInfo, 1)
	info[0] = s.getReadSummaryTableSingle(qname)

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(info)
}

func (s *Server) getReadSummaryTableSingle(qname string) *ReadOverviewInfo {

	qnameCluster := s.Handler.QnameCluster[qname]

	numClusterR1 := len(qnameCluster.ClusterR1.SimilarRecords)
	numClusterR2 := len(qnameCluster.ClusterR2.SimilarRecords)

	locations := make([]*ReadLocationInfo, 0, numClusterR1+numClusterR2)

	mappedBy := make(map[string]bool)

	for _, cR1 := range qnameCluster.ClusterR1.SimilarRecords {

		mappedByCluster := make(map[string]bool)
		for _, r := range cR1 {
			mappedBy[s.Handler.MapperInfos[r.MapperIndex].MapperName] = true
			mappedByCluster[s.Handler.MapperInfos[r.MapperIndex].MapperName] = true
		}
		mappedByClusterList := make([]string, 0, len(mappedByCluster))
		for mapperName := range mappedByCluster {
			mappedByClusterList = append(mappedByClusterList, mapperName)
		}
		slices.Sort(mappedByClusterList)

		locationInfo := ReadLocationInfo{
			Pair:          "first",
			Contig:        cR1[0].Rname,
			Strand:        !cR1[0].Flag.IsReverseStrand(),
			Position:      cR1[0].Pos,
			Cigar:         cR1[0].CigarObj.StringUniform(),
			CigarDetailed: cR1[0].CigarObj.String(),
			NumMismatches: cR1[0].CigarObj.GetNumMismatches(),
			NumGaps:       cR1[0].CigarObj.GetNumGaps(),
			NumMappedBy:   len(mappedByCluster),
			MappedBy:      mappedByClusterList,
			ReadIndices:   make([]int, 0),
		}

		for _, r := range cR1 {
			locationInfo.ReadIndices = append(locationInfo.ReadIndices, r.IndexInSam)
		}

		locations = append(locations, &locationInfo)
	}

	for _, cR2 := range qnameCluster.ClusterR2.SimilarRecords {

		mappedByCluster := make(map[string]bool)
		for _, r := range cR2 {
			mappedBy[s.Handler.MapperInfos[r.MapperIndex].MapperName] = true
			mappedByCluster[s.Handler.MapperInfos[r.MapperIndex].MapperName] = true
		}
		mappedByClusterList := make([]string, 0, len(mappedByCluster))
		for mapperName := range mappedByCluster {
			mappedByClusterList = append(mappedByClusterList, mapperName)
		}
		slices.Sort(mappedByClusterList)

		locationInfo := ReadLocationInfo{
			Pair:          "second",
			Contig:        cR2[0].Rname,
			Strand:        !cR2[0].Flag.IsReverseStrand(),
			Position:      cR2[0].Pos,
			Cigar:         cR2[0].CigarObj.StringUniform(),
			CigarDetailed: cR2[0].CigarObj.String(),
			NumMismatches: cR2[0].CigarObj.GetNumMismatches(),
			NumGaps:       cR2[0].CigarObj.GetNumGaps(),
			NumMappedBy:   len(mappedByCluster),
			MappedBy:      mappedByClusterList,
			ReadIndices:   make([]int, 0),
		}

		for _, r := range cR2 {
			locationInfo.ReadIndices = append(locationInfo.ReadIndices, r.IndexInSam)
		}

		locations = append(locations, &locationInfo)
	}

	mappedByList := make([]string, 0, len(mappedBy))
	for mapperName := range mappedBy {
		mappedByList = append(mappedByList, mapperName)
	}
	slices.Sort(mappedByList)

	distanceMean := (qnameCluster.ClusterR1.RecordDistanceMean + qnameCluster.ClusterR2.RecordDistanceMean) / 2.0

	return &ReadOverviewInfo{
		Qname:           qname,
		ReadLengthR1:    len(s.Handler.ReadNamesMap[qname].SequenceFirstOfPair),
		ReadLengthR2:    len(s.Handler.ReadNamesMap[qname].SequenceSecondOfPair),
		NumMappedBy:     len(mappedBy),
		MappedBy:        mappedByList,
		NumLocations:    len(locations),
		Locations:       locations,
		DistanceScore:   distanceMean,
		ConfidenceLevel: qnameCluster.ConfidenceLevel,
	}
}
