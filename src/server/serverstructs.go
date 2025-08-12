package server

import (
	"fmt"
	"log"
	"os"
	"strings"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/dataloader"
	"github.com/KleinSamuel/gtamap/src/formats/cigar"
	"github.com/KleinSamuel/gtamap/src/formats/fasta"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
)

type EnhancedRecord struct {
	sam.Record
	MapperIndex  int                        // index of the mapper in MapperInfos
	MappedGenome *regionvector.RegionVector // the matched intervals in the genome (mis- + matches (M, =, X in cigar))
	MappedRead   *regionvector.RegionVector // the matched intervals in the read (matches (M, =, X) in cigar)
	Mismatches   []int                      // the positions of the mismatches in the read sequence (min = 0, max = read length)
	CigarObj     *cigar.Object
}

func (h *MappingDataHandler) NewEnhancedRecord(record sam.Record, mapperIndex int) (*EnhancedRecord, error) {

	cigarObj, errCigarObj := cigar.NewCigarFromString(record.Cigar)
	if errCigarObj != nil {
		return nil, fmt.Errorf("error creating cigar from record: %w", errCigarObj)
	}

	genomicRegionsMapped := regionvector.NewRegionVector()
	readRegionsMapped := regionvector.NewRegionVector()
	start := record.Pos - 1 // SAM format is 1-based, convert to 0-based
	startRead := 0
	for _, elem := range cigarObj.Elements {
		if elem.IsMatch() {
			genomicRegionsMapped.AddRegion(start, start+elem.Length)
			readRegionsMapped.AddRegion(startRead, startRead+elem.Length)
		}
		if elem.ConsumesRead() {
			startRead += elem.Length
		}
		if elem.ConsumesReference() {
			start += elem.Length
		}
	}

	// number of genomic and read regions should be the same
	if genomicRegionsMapped.NumRegions() != readRegionsMapped.NumRegions() {
		logrus.Errorf("number of genomic regions (%d) does not match number of read regions (%d) in record %s", genomicRegionsMapped.NumRegions(), readRegionsMapped.NumRegions(), record.Qname)
	}

	genomicCombined, readCombined, errCombine := regionvector.CombineRegionVectorsConsecutiveInBoth(genomicRegionsMapped, readRegionsMapped)
	if errCombine != nil {
		fmt.Println(cigarObj.String())
		logrus.Errorf("error combining genomic and read regions: %v", errCombine)
		logrus.Errorf("record: %s for mapper index %d", record.Qname, mapperIndex)
		logrus.Errorf("genomic regions: %s", genomicRegionsMapped)
		logrus.Errorf("read regions: %s", readRegionsMapped)
	}

	// calculate mismatches based on inferred genome sequence, read sequence
	// and cigar string based region vectors
	posInRead := 0
	posInGenome := record.Pos - 1
	detailedCigarElements := make([]cigar.Element, 0)
	for _, elem := range cigarObj.Elements {

		if elem.Type == "S" || elem.Type == "H" || elem.Type == "I" {
			// skip soft-clipped, hard-clipped and inserted bases
			posInRead += elem.Length
			detailedCigarElements = append(detailedCigarElements, cigar.Element{
				Length: elem.Length,
				Type:   elem.Type,
			})
			continue
		}

		if elem.Type == "D" || elem.Type == "N" {
			// skip deleted bases in the genome
			posInGenome += elem.Length
			detailedCigarElements = append(detailedCigarElements, cigar.Element{
				Length: elem.Length,
				Type:   elem.Type,
			})
			continue
		}

		if elem.Type == "X" || elem.Type == "=" {
			posInGenome += elem.Length
			posInRead += elem.Length
			detailedCigarElements = append(detailedCigarElements, cigar.Element{
				Length: elem.Length,
				Type:   elem.Type,
			})
			continue
		}

		if elem.Type == "M" {

			// compare the read sequence with the genome sequence to find
			// mismatches, genome sequence is consecutive

			endInGenome := posInGenome + elem.Length

			isInTargetSeq := h.TargetRegion.Contig == record.Rname && h.TargetRegion.Start <= posInGenome && endInGenome <= h.TargetRegion.End

			var seq string

			if isInTargetSeq {
				seq = h.TargetRegion.SequenceDnaForward[posInGenome-h.TargetRegion.Start : endInGenome-h.TargetRegion.Start]
			} else {
				var errSeq error
				seq, errSeq = h.ExtractSequence(record.Rname, posInGenome, endInGenome)
				if errSeq != nil {
					return nil, fmt.Errorf("error extracting sequence from fasta: %w", errSeq)
				}
			}

			currentLen := 1
			currentIsMatch := record.Seq[posInRead] == seq[0]
			for i := 1; i < elem.Length; i++ {
				isMatch := record.Seq[posInRead+i] == seq[i]

				if currentIsMatch == isMatch {
					currentLen += 1
					continue
				}

				typeChar := "X"
				if currentIsMatch {
					typeChar = "="
				}

				detailedCigarElements = append(detailedCigarElements, cigar.Element{
					Length: currentLen,
					Type:   typeChar,
				})

				currentLen = 1
				currentIsMatch = isMatch
			}

			typeChar := "X"
			if currentIsMatch {
				typeChar = "="
			}

			detailedCigarElements = append(detailedCigarElements, cigar.Element{
				Length: currentLen,
				Type:   typeChar,
			})

			posInRead += elem.Length
			posInGenome += elem.Length

			continue
		}

		return nil, fmt.Errorf("unknown cigar element type: %s", elem.Type)
	}

	return &EnhancedRecord{
		Record:       record,
		MapperIndex:  mapperIndex,
		MappedGenome: genomicCombined,
		MappedRead:   readCombined,
		Mismatches:   make([]int, 0),
		CigarObj:     &cigar.Object{Elements: detailedCigarElements},
	}, nil
}

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
	RecordsByMapper [][]*EnhancedRecord            // outer array = mappers by order of MapperNames, inner array = records for each mapper
	QnamesByMapper  []map[string][]*EnhancedRecord // outer array = mappers by order of MapperNames, map has key = qname, value = enhanced record
	QnameCluster    map[string]*QnameCluster       // map has key = qname, value = cluster of records with this qname
}

type QnameCluster struct {
	Qname         string
	MapperPresent []bool     // boolean array indicating which mappers have this qname
	ClusterR1     *ClusterR1 // cluster of records for the first read in the pair (R1)
	ClusterR2     *ClusterR2 // cluster of records for the second read in the pair (R2)
}

type ClusterR1 struct {
	Records   []*EnhancedRecord // records of the first read in the pair (R1)
	Distances [][]float64       // distances between each pair of R1 records
}

type ClusterR2 struct {
	Records   []*EnhancedRecord // records of the second read in the pair (R2)
	Distances [][]float64       // distances between each pair of R2 records
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

		recordsByMapper = append(recordsByMapper, enhancedRecord)

		if _, exists := qnameMap[record.Qname]; !exists {
			qnameMap[record.Qname] = make([]*EnhancedRecord, 0)
		}
		qnameMap[record.Qname] = append(qnameMap[record.Qname], enhancedRecord)

		if _, exists := h.QnameCluster[record.Qname]; !exists {
			h.QnameCluster[record.Qname] = &QnameCluster{
				Qname: record.Qname,
				ClusterR1: &ClusterR1{
					Records: make([]*EnhancedRecord, 0),
				},
				ClusterR2: &ClusterR2{
					Records: make([]*EnhancedRecord, 0),
				},
			}
		}
		qnameCluster := h.QnameCluster[record.Qname]

		if enhancedRecord.Flag.IsFirstInPair() {
			qnameCluster.ClusterR1.Records = append(qnameCluster.ClusterR1.Records, enhancedRecord)
		} else {
			qnameCluster.ClusterR2.Records = append(qnameCluster.ClusterR2.Records, enhancedRecord)
		}
	}

	h.RecordsByMapper = append(h.RecordsByMapper, recordsByMapper)
	h.QnamesByMapper = append(h.QnamesByMapper, qnameMap)

	return nil
}

func (h *MappingDataHandler) ComputeQnameClusters() {

	for _, readInfo := range h.ReadInfos {

		if readInfo.Qname != "11" {
			continue
		}

		qnameCluster := h.QnameCluster[readInfo.Qname]

		fmt.Println("cluster: ", qnameCluster.Qname, " = ", len(qnameCluster.ClusterR1.Records), ",", len(qnameCluster.ClusterR2.Records))

		// compute distances R1

		for i := 0; i < len(qnameCluster.ClusterR1.Records); i++ {
			for j := i + 1; j < len(qnameCluster.ClusterR1.Records); j++ {
				dist, err := ComputeDistance(qnameCluster.ClusterR1.Records[i], qnameCluster.ClusterR1.Records[j])
				if err != nil {
					logrus.Errorf("error computing distance between records: %v", err)
					continue
				}
				fmt.Println(dist)
			}
		}
	}
}

func ComputeDistance(recordA *EnhancedRecord, recordB *EnhancedRecord) (float64, error) {

	fmt.Println(recordA.MapperIndex, recordA.Qname, recordA.Rname, recordA.Pos, recordA.Cigar)
	fmt.Println(recordB.MapperIndex, recordB.Qname, recordB.Rname, recordB.Pos, recordB.Cigar)

	if recordA == nil || recordB == nil {
		return 0, fmt.Errorf("one of the records is nil")
	}

	// records on different contigs have maximum distance
	if recordA.Rname != recordB.Rname {
		return 1, nil
	}

	// records on the same contig have a distance based on the overlap of their
	// mapped regions
	fmt.Println(recordA.MappedGenome)
	fmt.Println(recordA.MappedRead)
	fmt.Println(recordB.MappedGenome)
	fmt.Println(recordB.MappedRead)

	os.Exit(2)

	return 0.0, nil
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
