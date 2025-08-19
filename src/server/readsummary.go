package server

import "slices"

func (h *MappingDataHandler) GetReadSummaryTableData() []ReadOverviewInfo {

	reads := make([]ReadOverviewInfo, 0, len(h.ReadInfos))

	for _, readInfo := range h.ReadInfos {

		qnameCluster := h.QnameCluster[readInfo.Qname]

		numClusterR1 := len(qnameCluster.ClusterR1.SimilarRecords)
		numClusterR2 := len(qnameCluster.ClusterR2.SimilarRecords)

		locations := make([]*ReadLocationInfo, 0, numClusterR1+numClusterR2)

		mappedBy := make(map[string]bool)

		isAcceptedR1 := false
		isAcceptedR2 := false

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
			}

			if cR1[0].IsAccepted {
				isAcceptedR1 = true
			}

			for _, r := range cR1 {
				locationInfo.ReadIndices = append(locationInfo.ReadIndices, r.Index)
				locationInfo.ReadIndicesInMappers = append(locationInfo.ReadIndicesInMappers, r.IndexInSam)
			}

			locations = append(locations, &locationInfo)
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
			}

			if cR2[0].IsAccepted {
				isAcceptedR2 = true
			}

			for _, r := range cR2 {
				locationInfo.ReadIndices = append(locationInfo.ReadIndices, r.Index)
				locationInfo.ReadIndicesInMappers = append(locationInfo.ReadIndicesInMappers, r.IndexInSam)
			}

			locations = append(locations, &locationInfo)
		}

		mappedByList := make([]string, 0, len(mappedBy))
		for mapperName := range mappedBy {
			mappedByList = append(mappedByList, mapperName)
		}
		slices.Sort(mappedByList)

		distanceMean := (qnameCluster.ClusterR1.RecordDistanceMean + qnameCluster.ClusterR2.RecordDistanceMean) / 2.0

		readOverviewInfo := ReadOverviewInfo{
			Qname:           readInfo.Qname,
			ReadLengthR1:    len(readInfo.SequenceFirstOfPair),
			ReadLengthR2:    len(readInfo.SequenceSecondOfPair),
			NumMappedBy:     len(mappedBy),
			MappedBy:        mappedByList,
			NumLocations:    len(locations),
			Locations:       locations,
			DistanceScore:   distanceMean,
			ConfidenceLevel: qnameCluster.ConfidenceLevel,
			IsAccepted:      isAcceptedR1 && isAcceptedR2,
		}

		reads = append(reads, readOverviewInfo)
	}

	return reads
}
