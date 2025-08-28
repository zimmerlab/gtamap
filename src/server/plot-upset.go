package server

import (
	"sort"
	"strconv"

	"github.com/KleinSamuel/gtamap/src/config"
)

type UpsetElement struct {
	Name string   `json:"name"`
	Sets []string `json:"sets"`
}

func (h *MappingDataHandler) GetUpsetData(onlyTargetRegion bool, usePosition bool, useCigar bool) []UpsetElement {

	data := make([]UpsetElement, 0)

	for _, qnameClust := range h.QnameCluster {

		set := make(map[string]map[string]bool)

		for _, c1 := range qnameClust.ClusterR1.Records {

			// skip reads that do not map to the target region if flag is set
			if onlyTargetRegion && (c1.Rname != config.GetTargetContig() ||
				!c1.MappedGenome.Overlaps(config.GetTargetStart(), config.GetTargetEnd())) {
				continue
			}

			readKey := getReadKey(c1, "R1", usePosition, useCigar)

			if _, exists := set[readKey]; !exists {
				set[readKey] = make(map[string]bool)
			}
			set[readKey][h.MapperInfos[c1.MapperIndex].MapperName] = true
		}
		for _, c2 := range qnameClust.ClusterR2.Records {
			// skip reads that do not map to the target region if flag is set
			if onlyTargetRegion && (c2.Rname != config.GetTargetContig() ||
				!c2.MappedGenome.Overlaps(config.GetTargetStart(), config.GetTargetEnd())) {
				continue
			}

			readKey := getReadKey(c2, "R2", usePosition, useCigar)

			if _, exists := set[readKey]; !exists {
				set[readKey] = make(map[string]bool)
			}
			set[readKey][h.MapperInfos[c2.MapperIndex].MapperName] = true
		}

		for qname, mapperNames := range set {

			list := make([]string, 0, len(mapperNames))
			for name := range mapperNames {
				list = append(list, name)
			}
			sort.Strings(list)

			data = append(data, UpsetElement{
				Name: qname,
				Sets: list,
			})
		}
	}

	return data
}

func getReadKey(record *EnhancedRecord, prefix string, usePosition bool, useCigar bool) string {
	readKey := prefix + "::" + record.Qname
	if usePosition {
		readKey += "::" + strconv.Itoa(record.Pos)
	}
	if useCigar {
		readKey += "::" + record.UniformCigar()
	}
	return readKey
}
