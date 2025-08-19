package server

import (
	"fmt"
	"log"
	"math"
	"os"
	"slices"
	"strings"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/dataloader"
	"github.com/KleinSamuel/gtamap/src/formats/fasta"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
)

type MapperInfo struct {
	Index      int // index in MapperNames
	MapperName string
	Target     string // target | genome
}

type ReadInfo struct {
	Qname                string // name of the read
	Index                int    // index in ReadNames
	IsPaired             bool   // true if read is paired-end
	SequenceFirstOfPair  string // sequence of the first read in the pair
	QualityFirstOfPair   string // quality of the first read in the pair
	SequenceSecondOfPair string // sequence of the second read in the pair
	QualitySecondOfPair  string // quality of the second read in the pair
}

func NewReadInfo() *ReadInfo {
	return &ReadInfo{
		Qname:                "",
		Index:                -1,
		IsPaired:             false,
		SequenceFirstOfPair:  "",
		QualityFirstOfPair:   "",
		SequenceSecondOfPair: "",
		QualitySecondOfPair:  "",
	}
}

type TargetRegion struct {
	Contig             string // name of the contig
	ContigLength       int    // the length of the original contig
	Start              int    // start position (0-based)
	End                int    // end position (0-based, exclusive)
	IsForwardStrand    bool   // true if the target region is on the forward strand, false if on the reverse strand
	SequenceDnaForward string // sequence of the target region
	SequenceDnaReverse string // reverse complement sequence of the target region
}

type MappingDataHandler struct {
	FastaFile       *os.File
	FastaIndex      *fasta.Index
	TargetRegion    *TargetRegion
	MapperInfos     []*MapperInfo                  // list of mapper infos in consistent order
	ReadInfos       []*ReadInfo                    // list of read infos
	ReadNamesMap    map[string]*ReadInfo           // map has key = read name, value = read info
	Records         []*EnhancedRecord              // list of all records
	RecordsByMapper [][]*EnhancedRecord            // outer array = mappers by order of MapperNames, inner array = records for each mapper
	QnamesByMapper  []map[string][]*EnhancedRecord // outer array = mappers by order of MapperNames, map has key = qname, value = enhanced record
	QnameCluster    map[string]*QnameCluster       // map has key = qname, value = cluster of records with this qname
}

type QnameCluster struct {
	Qname           string
	MapperPresent   []bool       // boolean array indicating which mappers have this qname
	ClusterR1       *ReadCluster // cluster of records for the first read in the pair (R1)
	ClusterR2       *ReadCluster // cluster of records for the second read in the pair (R2)
	ConfidenceLevel int
}

type ReadCluster struct {
	Records            []*EnhancedRecord // records of the first read in the pair
	Distances          [][]float64       // distances between each pair of R1 records
	DistancesMapper    [][]float64       // distances between each pair of mapper
	MapperDistanceMean float64
	RecordDistanceMean float64
	ConfidenceLevel    int
	SimilarRecords     [][]*EnhancedRecord // similar records for the first read in the pair
	AcceptedRecord     *EnhancedRecord     // the record that is accepted as the representative of the cluster
}

func NewMappingDataHandler(fastaFilePath string, fastaIndexFilePath string) (*MappingDataHandler, error) {

	fastaIndex := fasta.ReadFastaIndexUsingPath(fastaIndexFilePath)
	if fastaIndex == nil {
		return nil, fmt.Errorf("could not read fasta index: %s", fastaIndexFilePath)
	}

	fastaFile, fastaFileErr := os.Open(fastaFilePath)
	if fastaFileErr != nil {
		return nil, fmt.Errorf("error opening fasta file: %w", fastaFileErr)
	}
	// defer fastaFile.Close()

	h := &MappingDataHandler{
		FastaFile:       fastaFile,
		FastaIndex:      fastaIndex,
		MapperInfos:     make([]*MapperInfo, 0),
		ReadInfos:       make([]*ReadInfo, 0),
		ReadNamesMap:    make(map[string]*ReadInfo),
		RecordsByMapper: make([][]*EnhancedRecord, 0),
		QnamesByMapper:  make([]map[string][]*EnhancedRecord, 0),
		QnameCluster:    make(map[string]*QnameCluster),
	}

	seqFw, err := h.ExtractSequence(config.GetTargetContig(), config.GetTargetStart(), config.GetTargetEnd())
	if err != nil {
		return nil, fmt.Errorf("error extracting sequence from fasta: %w", err)
	}
	seqRv := utils.ReverseComplementDNA(seqFw)

	targetRegion := &TargetRegion{
		Contig:             config.GetTargetContig(),
		Start:              config.GetTargetStart(),
		End:                config.GetTargetEnd(),
		IsForwardStrand:    config.GetTargetStrand() == "+",
		SequenceDnaForward: seqFw,
		SequenceDnaReverse: seqRv,
	}

	h.TargetRegion = targetRegion

	return h, nil
}

func (h *MappingDataHandler) ExtractSequence(contig string, start int, end int) (string, error) {

	if h.FastaIndex == nil {
		panic("Fasta index is not initialized")
	}
	if contig == "" || start < 0 || end < 0 || start >= end {
		return "", fmt.Errorf("invalid parameters: contig='%s', start=%d, end=%d", contig, start, end)
	}
	// check if contig exists in fasta index
	if _, exists := h.FastaIndex.Entries[contig]; !exists {
		return "", fmt.Errorf("contig '%s' not found in fasta index", contig)
	}
	// extract sequence from fasta file
	seq := dataloader.ExtractSequenceAsStringFromFasta(h.FastaFile, h.FastaIndex, contig, uint32(start), uint32(end))
	return seq, nil
}

func (h *MappingDataHandler) GetMapperInfoByIndex(index int) *MapperInfo {
	if index < 0 || index >= len(h.MapperInfos) {
		return nil
	}
	return h.MapperInfos[index]
}

func (h *MappingDataHandler) GetMapperInfoByNameAndTarget(mapperName string, target string) *MapperInfo {
	for _, info := range h.MapperInfos {
		if info.MapperName == mapperName && info.Target == target {
			return info
		}
	}
	return nil
}

func (h *MappingDataHandler) AddReadInfo(record sam.Record, mapperIndex int) *ReadInfo {

	readInfo := h.ReadNamesMap[record.Qname]

	// check and return if the info about this read was already added
	if readInfo != nil {
		if readInfo.SequenceFirstOfPair != "" && record.Flag.IsFirstInPair() {
			return readInfo
		}
		if readInfo.SequenceSecondOfPair != "" && record.Flag.IsSecondInPair() {
			return readInfo
		}
	} else {
		readInfo = NewReadInfo()
		h.ReadInfos = append(h.ReadInfos, readInfo)
		h.ReadNamesMap[record.Qname] = readInfo
		readInfo.Qname = record.Qname
		readInfo.IsPaired = record.Flag.IsPaired()
	}

	if record.Flag.IsFirstInPair() {
		readInfo.SequenceFirstOfPair = record.Seq
		readInfo.QualityFirstOfPair = record.Qual
	} else {
		readInfo.SequenceSecondOfPair = record.Seq
		readInfo.QualitySecondOfPair = record.Qual
	}

	return readInfo
}

func (h *MappingDataHandler) AddMapperInfo(mapperName string, mapperSamFilePath string, target string) error {

	// check if mapper with this target was already added
	mapperInfo := h.GetMapperInfoByNameAndTarget(mapperName, target)

	if mapperInfo != nil {
		return fmt.Errorf("mapper with name '%s' and target '%s' already exists", mapperName, target)
	}

	// add new mapper info
	mapperInfo = &MapperInfo{
		Index:      len(h.MapperInfos),
		MapperName: mapperName,
		Target:     target,
	}
	h.MapperInfos = append(h.MapperInfos, mapperInfo)

	// add records of this mapper
	parsedFile, err := sam.ParseSAMFile(mapperSamFilePath)
	if err != nil {
		return fmt.Errorf("error reading sam file: %w", err)
	}

	recordsByMapper := make([]*EnhancedRecord, 0, len(parsedFile.Records))
	qnameMap := make(map[string][]*EnhancedRecord, len(parsedFile.Records))

	for _, record := range parsedFile.Records {

		h.AddReadInfo(record, mapperInfo.Index)

		if record.Rname == "*" {
			continue // skip records that are unmapped
		}

		enhancedRecord, errRecord := h.NewEnhancedRecord(record, mapperInfo.Index)
		if errRecord != nil {
			return fmt.Errorf("error creating enhanced record: %w", errRecord)
		}

		// check records for duplicates
		if len(recordsByMapper) > 0 {
			duplicateFound := false
			for _, existingRecord := range recordsByMapper {
				if existingRecord.IsEqual(enhancedRecord) {
					duplicateFound = true
					break
				}
			}
			if duplicateFound {
				logrus.Warn("duplicate record found: ", record.Qname, " in mapper: ", mapperName)
				continue // skip duplicate records
			}
		}

		// add the records
		enhancedRecord.Index = len(h.Records)
		h.Records = append(h.Records, enhancedRecord)

		recordsByMapper = append(recordsByMapper, enhancedRecord)

		if _, exists := qnameMap[record.Qname]; !exists {
			qnameMap[record.Qname] = make([]*EnhancedRecord, 0)
		}
		qnameMap[record.Qname] = append(qnameMap[record.Qname], enhancedRecord)

		if _, exists := h.QnameCluster[record.Qname]; !exists {
			h.QnameCluster[record.Qname] = &QnameCluster{
				Qname: record.Qname,
				ClusterR1: &ReadCluster{
					Records:         make([]*EnhancedRecord, 0),
					SimilarRecords:  make([][]*EnhancedRecord, 0),
					Distances:       make([][]float64, 0),
					DistancesMapper: make([][]float64, 0),
				},
				ClusterR2: &ReadCluster{
					Records:         make([]*EnhancedRecord, 0),
					SimilarRecords:  make([][]*EnhancedRecord, 0),
					Distances:       make([][]float64, 0),
					DistancesMapper: make([][]float64, 0),
				},
			}
		}
		qnameCluster := h.QnameCluster[record.Qname]

		if enhancedRecord.Flag.IsFirstInPair() {
			qnameCluster.ClusterR1.Records = append(qnameCluster.ClusterR1.Records, enhancedRecord)
		} else {
			qnameCluster.ClusterR2.Records = append(qnameCluster.ClusterR2.Records, enhancedRecord)
		}

		enhancedRecord.QnameCluster = qnameCluster
	}

	h.RecordsByMapper = append(h.RecordsByMapper, recordsByMapper)
	h.QnamesByMapper = append(h.QnamesByMapper, qnameMap)

	return nil
}

func (h *MappingDataHandler) ComputeQnameClusters() {

	for _, readInfo := range h.ReadInfos {

		// if readInfo.Qname != "11" {
		// 	continue
		// }

		cluster := h.QnameCluster[readInfo.Qname]

		// fmt.Println("cluster: ", cluster.Qname, " = ", len(cluster.ClusterR1.Records), ",", len(cluster.ClusterR2.Records))

		h.ComputeReadClusters(cluster.ClusterR1, true)
		h.ComputeReadClusters(cluster.ClusterR2, false)

		cluster.ConfidenceLevel = int(math.Min(float64(cluster.ClusterR1.ConfidenceLevel), float64(cluster.ClusterR2.ConfidenceLevel)))

		// PrintMatrix2D(cluster.ClusterR1.Distances)
		// PrintMatrix2D(cluster.ClusterR1.DistancesMapper)
		// fmt.Println("score: ", cluster.ClusterR1.RecordDistanceMean)
		// fmt.Println("record distance mean: ", cluster.ClusterR1.MapperDistanceMean)

		// PrintMatrix2D(cluster.ClusterR2.Distances)
		// PrintMatrix2D(cluster.ClusterR2.DistancesMapper)
		// fmt.Println("score: ", cluster.ClusterR2.RecordDistanceMean)
		// fmt.Println("record distance mean: ", cluster.ClusterR2.MapperDistanceMean)

		// fmt.Println(cluster.ClusterR1.RecordDistanceMean, cluster.ClusterR2.RecordDistanceMean, (cluster.ClusterR1.RecordDistanceMean+cluster.ClusterR2.RecordDistanceMean)/2)

	}
}

func (h *MappingDataHandler) ComputeReadClusters(cluster *ReadCluster, isR1 bool) {

	for i := 0; i < len(cluster.Records); i++ {
		cluster.Distances = append(cluster.Distances, make([]float64, len(cluster.Records)))
	}

	// distances between mappers
	for i := 0; i < len(h.MapperInfos); i++ {
		cluster.DistancesMapper = append(cluster.DistancesMapper, make([]float64, len(h.MapperInfos)))
	}
	for i := 0; i < len(h.MapperInfos); i++ {
		for j := i + 1; j < len(h.MapperInfos); j++ {
			cluster.DistancesMapper[i][j] = 1.0
			cluster.DistancesMapper[j][i] = 1.0
		}
	}

	for i := 0; i < len(cluster.Records); i++ {
		for j := i + 1; j < len(cluster.Records); j++ {
			dist, err := h.ComputeDistance(cluster.Records[i], cluster.Records[j], isR1)
			if err != nil {
				logrus.Errorf("error computing distance between records: %v", err)
				continue
			}

			cluster.Distances[i][j] = dist
			cluster.Distances[j][i] = dist

			mIndexA := cluster.Records[i].MapperIndex
			mIndexB := cluster.Records[j].MapperIndex

			if mIndexA != mIndexB {
				if cluster.DistancesMapper[mIndexA][mIndexB] > dist {
					cluster.DistancesMapper[mIndexA][mIndexB] = dist
					cluster.DistancesMapper[mIndexB][mIndexA] = dist
				}
			}
		}
	}

	sumDistRecord := 0.0
	for i := 0; i < len(cluster.Records); i++ {
		for j := i + 1; j < len(cluster.Records); j++ {
			sumDistRecord += cluster.Distances[i][j]
		}
	}
	cluster.RecordDistanceMean = sumDistRecord / float64(len(cluster.Records)*(len(cluster.Records)-1)/2)

	sumDistMapper := 0.0
	for i := 0; i < len(h.MapperInfos); i++ {
		for j := i + 1; j < len(h.MapperInfos); j++ {
			sumDistMapper += cluster.DistancesMapper[i][j]
		}
	}
	cluster.MapperDistanceMean = sumDistMapper / float64(len(h.MapperInfos)*(len(h.MapperInfos)-1)/2)

	// fmt.Println("---")
	// for i, r := range cluster.Records {
	// 	fmt.Println(i, r.MappedRead, r.MappedGenome)
	// }
	// PrintMatrix2D(cluster.Distances)

	// find same records
	visited := make(map[int]bool)

	for i := 0; i < len(cluster.Records); i++ {

		if visited[i] {
			continue
		}

		// do dfs to find connected component with record at index i
		toVisit := make([]int, 0)
		toVisit = append(toVisit, i)
		toVisitIndex := 0

		for toVisitIndex < len(toVisit) {
			j := toVisit[toVisitIndex]
			visited[j] = true

			for k := 0; k < len(cluster.Records); k++ {
				if k == j || visited[k] || cluster.Distances[j][k] > 0 {
					continue
				}
				toVisit = append(toVisit, k)
				break
			}
			toVisitIndex++
		}

		// fmt.Println(toVisit)

		simRecords := make([]*EnhancedRecord, 0)
		for _, idx := range toVisit {
			simRecords = append(simRecords, cluster.Records[idx])
		}
		cluster.SimilarRecords = append(cluster.SimilarRecords, simRecords)
	}

	// compute confidence level
	// 5 = all mappers produce exactly one record which is shared
	// 4 = all mappers produce multiple records which are all shared
	// 3 = all mappers produce one shared record but also have others
	// 2 = there is a shared record by multiple mappers
	// 1 = there are no shared records between any mapper

	minLenCluster := len(h.MapperInfos)
	maxLenCluster := 0
	for _, simRecords := range cluster.SimilarRecords {
		if len(simRecords) > maxLenCluster {
			maxLenCluster = len(simRecords)
		}
		if len(simRecords) < minLenCluster {
			minLenCluster = len(simRecords)
		}
	}

	if maxLenCluster == len(h.MapperInfos) {
		// there is a record that is shared by all mappers
		if len(cluster.SimilarRecords) == 1 {
			// all mappers produce exactly one record which is shared
			cluster.ConfidenceLevel = 5
		} else if minLenCluster == len(h.MapperInfos) {
			// all mappers produce multiple records which are all shared
			cluster.ConfidenceLevel = 4
		} else {
			// all mappers produce one shared record but also have others
			cluster.ConfidenceLevel = 3
		}
	} else if maxLenCluster > 1 {
		// there is a shared record by multiple mappers, but not all
		cluster.ConfidenceLevel = 2
	} else {
		// there is no shared record
		cluster.ConfidenceLevel = 1
	}
}

func PrintMatrix2D(matrix [][]float64) {
	for _, row := range matrix {
		for _, val := range row {
			// Format to 2 decimal places and add a tab
			fmt.Printf("%.2f\t", val)
		}
		fmt.Println()
	}
}

func (h *MappingDataHandler) ComputeDistance(recordA *EnhancedRecord, recordB *EnhancedRecord, isR1 bool) (float64, error) {

	weightMappedDifferent := 0.5
	weightMappedSame := 1.0
	weightGapBounds := 0.5

	// fmt.Println(recordA.MapperIndex, recordA.Qname, recordA.Rname, recordA.Pos, recordA.Cigar)
	// fmt.Println(recordB.MapperIndex, recordB.Qname, recordB.Rname, recordB.Pos, recordB.Cigar)

	if recordA == nil || recordB == nil {
		return 0, fmt.Errorf("one of the records is nil")
	}

	// records on different contigs have maximum distance
	if recordA.Rname != recordB.Rname {
		return 1, nil
	}

	// fmt.Println(recordA.MappedGenome.Regions)
	// fmt.Println(recordA.MappedRead.Regions)
	// fmt.Println(recordB.MappedGenome.Regions)
	// fmt.Println(recordB.MappedRead.Regions)

	// records on the same contig have a distance based on the overlap of their
	// mapped regions
	var readLength int
	if isR1 {
		readLength = len(h.ReadNamesMap[recordA.Qname].SequenceFirstOfPair)
	} else {
		readLength = len(h.ReadNamesMap[recordA.Qname].SequenceSecondOfPair)
	}

	regionIndexA := 0
	subIndexA := 0
	regionIndexB := 0
	subIndexB := 0

	numMappedOnlyA := 0
	numMappedOnlyB := 0
	// number of positions in the read that were mapped in both records
	numMappedBoth := 0
	// number of positions that were mapped in both records and to the same
	// genomic position
	numMappedSame := 0

	doneA := false
	doneB := false

	for i := 0; i < readLength; i++ {

		for !doneA && recordA.MappedRead.Regions[regionIndexA].Start+subIndexA < i {
			subIndexA++
			if subIndexA >= recordA.MappedRead.Regions[regionIndexA].Length() {
				subIndexA = 0
				regionIndexA++
			}
			if regionIndexA >= len(recordA.MappedRead.Regions) {
				doneA = true
				break
			}
		}

		isMappedA := !doneA && recordA.MappedRead.Regions[regionIndexA].Start+subIndexA == i
		// if !isMappedA {
		// 	continue
		// }

		for !doneB && recordB.MappedRead.Regions[regionIndexB].Start+subIndexB < i {
			subIndexB++
			if subIndexB >= recordB.MappedRead.Regions[regionIndexB].Length() {
				subIndexB = 0
				regionIndexB++
			}
			if regionIndexB >= len(recordB.MappedRead.Regions) {
				doneB = true
				break
			}
		}

		isMappedB := !doneB && recordB.MappedRead.Regions[regionIndexB].Start+subIndexB == i
		// if !isMappedB {
		// 	continue
		// }

		// if done {
		// 	break
		// }

		if !isMappedA && !isMappedB {
			continue
		}

		if isMappedA && !isMappedB {
			numMappedOnlyA++
		} else if !isMappedA && isMappedB {
			numMappedOnlyB++
		} else {

			genomePosA := recordA.MappedGenome.Regions[regionIndexA].Start + subIndexA
			genomePosB := recordB.MappedGenome.Regions[regionIndexB].Start + subIndexB

			numMappedBoth++
			if genomePosA == genomePosB {
				numMappedSame++
			}
		}
	}

	percentMappedOnlyOne := float64(numMappedOnlyA+numMappedOnlyB) / float64(readLength)

	percentMappedBoth := 0.0
	if numMappedBoth > 0 {
		percentMappedBoth = float64(numMappedSame) / float64(numMappedBoth)
	}

	impliedGapsA := recordA.MappedGenome.GetGaps()
	impliedGapsA.SortInPlace()
	impliedGapsB := recordB.MappedGenome.GetGaps()
	impliedGapsB.SortInPlace()

	// fmt.Println("implied gaps A:", impliedGapsA)
	// fmt.Println("implied gaps B:", impliedGapsB)

	numGapsA := len(impliedGapsA.Regions)
	numGapsB := len(impliedGapsB.Regions)

	numGapStartsSame := 0
	numGapEndsSame := 0

	for _, gapA := range impliedGapsA.Regions {
		startOk := false
		endOk := false
		for _, gapB := range impliedGapsB.Regions {
			if !startOk && gapA.Start == gapB.Start {
				numGapStartsSame++
				startOk = true
			}
			if !endOk && gapA.End == gapB.End {
				numGapEndsSame++
				endOk = true
			}
			if startOk && endOk {
				break
			}
		}
	}

	percentGapBoundsSame := 1.0
	if (numGapsA + numGapsB) > 0 {
		percentGapBoundsSame = float64(numGapStartsSame+numGapEndsSame) / float64(numGapsA+numGapsB)
	}

	// fmt.Println("percent mapped only one:", percentMappedOnlyOne)
	// fmt.Println("both same:", percentMappedBoth)
	// fmt.Println("percent gaps", percentGapBoundsSame)

	dist := ((percentMappedOnlyOne * weightMappedDifferent) +
		((1 - percentMappedBoth) * weightMappedSame) +
		((1 - percentGapBoundsSame) * weightGapBounds)) /
		3.0

	return dist, nil
}

type MapperResult struct {
	Name   string `json:"name"`
	Path   string `json:"path"`
	Target string `json:"target"`
}

// GetMapperResults retrieves the results of mappers from the output directory.
// It assumes that the output directory contains files named in the format "mapperName.target.sam".
func (h *MappingDataHandler) GetMapperResults() []MapperResult {

	results := make([]MapperResult, 0)

	files, err := os.ReadDir(config.GetRunDir())
	if err != nil {
		logrus.Error("Error reading output directory: ", err)
		return results
	}

	for _, file := range files {
		if strings.HasSuffix(file.Name(), ".sam") {
			parts := strings.Split(file.Name(), ".")

			if len(parts) < 3 {
				logrus.Warn("Skipping mapper file: ", file.Name())
				continue
			}

			result := MapperResult{
				Name:   parts[1],
				Path:   config.GetRunDir() + "/" + file.Name(),
				Target: parts[2],
			}

			results = append(results, result)
		}
	}

	return results
}

func (h *MappingDataHandler) LoadMapperResults() {

	results := h.GetMapperResults()

	for _, result := range results {

		err := h.AddMapperInfo(result.Name, result.Path, result.Target)

		if err != nil {
			log.Fatalf("Error initializing analysis service: %v", err)
		}

		logrus.WithFields(logrus.Fields{
			"mapperName":        result.Name,
			"mapperSamFilePath": result.Path,
			"target":            result.Target,
		}).Info("added mapper result")
	}

	h.ComputeQnameClusters()
}

func (h *MappingDataHandler) GetAcceptedRecords() []*EnhancedRecord {

	accepted := make([]*EnhancedRecord, 0)

	for _, cluster := range h.QnameCluster {

		if cluster.ClusterR1.AcceptedRecord != nil {
			accepted = append(accepted, cluster.ClusterR1.AcceptedRecord)
		}
		if cluster.ClusterR2.AcceptedRecord != nil {
			accepted = append(accepted, cluster.ClusterR2.AcceptedRecord)
		}
	}

	return accepted
}

func (h *MappingDataHandler) AcceptRecord(r *EnhancedRecord) {

	logrus.WithFields(logrus.Fields{
		"qname": r.Qname,
		"index": r.Index,
	}).Info("accepting record")

	r.IsAccepted = true

	if r.Flag.IsFirstInPair() {
		r.QnameCluster.ClusterR1.AcceptedRecord = r
	} else {
		r.QnameCluster.ClusterR2.AcceptedRecord = r
	}
}

func (h *MappingDataHandler) UnacceptRecord(r *EnhancedRecord) {

	logrus.WithFields(logrus.Fields{
		"qname": r.Qname,
		"index": r.Index,
	}).Info("unaccepting record")

	r.IsAccepted = false

	if r.Flag.IsFirstInPair() {
		r.QnameCluster.ClusterR1.AcceptedRecord = nil
	} else {
		r.QnameCluster.ClusterR2.AcceptedRecord = nil
	}
}

func (h *MappingDataHandler) GetSimilarRecordsInCluster(r *EnhancedRecord) []*EnhancedRecord {

	var clust *ReadCluster
	if r.Flag.IsFirstInPair() {
		clust = r.QnameCluster.ClusterR1
	} else {
		clust = r.QnameCluster.ClusterR2
	}

	for _, s1 := range clust.SimilarRecords {
		contained := slices.Contains(s1, r)
		if !contained {
			continue
		}
		return s1
	}

	return nil
}

func (h *MappingDataHandler) GetMapperNamesThatMappedRecord(r *EnhancedRecord) []string {

	var clust *ReadCluster
	if r.Flag.IsFirstInPair() {
		clust = r.QnameCluster.ClusterR1
	} else {
		clust = r.QnameCluster.ClusterR2
	}

	names := make([]string, 0)

	for _, s1 := range clust.SimilarRecords {
		contained := slices.Contains(s1, r)
		if !contained {
			continue
		}
		for _, r := range s1 {
			mapperName := h.MapperInfos[r.MapperIndex].MapperName
			names = append(names, mapperName)
		}
	}

	return names
}
