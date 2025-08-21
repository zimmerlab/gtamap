package server

import (
	"fmt"
	"strings"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/utils"
)

func (h *MappingDataHandler) GetReadSelectionTableData() {

}

func (h *MappingDataHandler) GetAcceptedRecordsSam(relativeCoords bool) string {

	var builder strings.Builder

	// read contig length from fasta index
	contigLength := int(h.FastaIndex.Entries[config.GetTargetContig()].Length)
	// target offset within contig to calculate relative coordinates
	targetOffset := 0
	if relativeCoords {
		contigLength = config.GetTargetEnd() - config.GetTargetStart() + 1
		targetOffset = config.GetTargetStart()
	}

	builder.WriteString("@HD\tVN:1.6\tSO:unsorted\n")
	builder.WriteString(fmt.Sprintf("@SQ\tSN:%s\tLN:%d\n", config.GetTargetContig(), contigLength))

	for _, readInfo := range h.ReadInfos {
		qclust := h.QnameCluster[readInfo.Qname]
		acceptedR1 := qclust.ClusterR1.AcceptedRecord
		acceptedR2 := qclust.ClusterR2.AcceptedRecord

		if acceptedR1 != nil && IsInTargetRegion(acceptedR1) == 1 {

			similarRecords := h.GetSimilarRecordsInCluster(acceptedR1)
			sNames := make([]string, 0, len(similarRecords))
			rIndices := make([]int, 0, len(similarRecords))
			for _, r := range similarRecords {
				sNames = append(sNames, h.MapperInfos[r.MapperIndex].MapperName)
				rIndices = append(rIndices, r.Index)
			}
			sNamesString := utils.ArrayStringToString(sNames, ",")
			rIndicesString := utils.ArrayIntToString(rIndices, ",")

			builder.WriteString(acceptedR1.StringFinal(targetOffset, sNamesString, rIndicesString))
			builder.WriteString("\n")
		}
		if acceptedR2 != nil && IsInTargetRegion(acceptedR2) == 1 {

			similarRecords := h.GetSimilarRecordsInCluster(acceptedR2)
			sNames := make([]string, 0, len(similarRecords))
			rIndices := make([]int, 0, len(similarRecords))
			for _, r := range similarRecords {
				sNames = append(sNames, h.MapperInfos[r.MapperIndex].MapperName)
				rIndices = append(rIndices, r.Index)
			}
			sNamesString := utils.ArrayStringToString(sNames, ",")
			rIndicesString := utils.ArrayIntToString(rIndices, ",")

			builder.WriteString(acceptedR2.StringFinal(targetOffset, sNamesString, rIndicesString))
			builder.WriteString("\n")
		}
	}

	return builder.String()
}
