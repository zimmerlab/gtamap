package secondpass

import (
	"bufio"
	"errors"
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/interval"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"io"
	"os"
	"strconv"
)

func InferIntrons(samPath string) {
	gapsPerSeq := ReadSAMIntrons(samPath)
	countedGapsPerSeq := countGaps(gapsPerSeq)
	intronClusters := clusterGaps(countedGapsPerSeq)
	fmt.Println(intronClusters)

	file, err := os.Open(samPath)
	if err != nil {
		panic("Error during sam reading")
	}
	samIndex := CreateSAMIndex(file)

	for _, entry := range samIndex {
		fw, _ := ParseSparseSAMEntry(file, entry.first)
		rv, _ := ParseSparseSAMEntry(file, entry.second)
		fmt.Println(fw)
		fmt.Println(rv)
	}

}

func countGaps(gapsPerSeq map[int][]*interval.Interval) map[int]map[interval.Interval]int {
	gapMap := make(map[int]map[interval.Interval]int)
	// count unique gapsPerSeq per SI
	for sid, gaps := range gapsPerSeq {
		for _, g := range gaps {
			// check if we already have seen the sid
			_, exists := gapMap[sid]
			if !exists {
				// if not, init new map
				gapMap[sid] = make(map[interval.Interval]int)
			}
			// check if we have already seen the curr gap in sid
			gap := interval.Interval{Start: g.Start, End: g.End}
			_, exists = gapMap[sid][gap]
			if !exists {
				// init new count with 1
				gapMap[sid][gap] = 1
			} else {
				// if the curr gap is already in gapMap[sid], increment counter
				gapMap[sid][gap] = gapMap[sid][gap] + 1
			}
		}
	}
	return gapMap
}

func ReadSAMIntrons(samfile string) map[int][]*interval.Interval {

	file, err := os.Open(samfile)
	if err != nil {
		panic("Error during sam reading")
	}
	scanner := bufio.NewScanner(file)

	gaps := make(map[int][]*interval.Interval, 0)

	for scanner.Scan() {
		line := scanner.Bytes()

		if len(line) == 0 {
			continue
		}

		if line[0] == '@' {
			continue
		}

		numTabs := 0
		var readStart []byte
		var cigar []byte
		var sidStart int
		var sidEnd int

		for i := 0; i < len(line); i++ {
			if line[i] == '\t' {
				numTabs++
				continue
			}

			// parse start coords of read
			if numTabs == 3 {
				start := i
				for line[i] != '\t' {
					i++
				}
				stop := i
				readStart = make([]byte, stop-start)
				copy(readStart, line[start:stop])
				i++
				numTabs++
			}

			// parse cigar string
			if numTabs == 5 {
				start := i
				for line[i] != '\t' {
					i++
				}
				stop := i
				cigar = make([]byte, stop-start)
				copy(cigar, line[start:stop])

				i++
				numTabs++
			}

			// parse SI attribute
			if numTabs == 11 {
				i = i + 5 // skip chars 'SI:i:'
				sidStart = i
				// parse all chars after SI:i:...
				for i < len(line) {
					i++
				}
				sidEnd = i
				break
			}
		}

		sidInt, err := strconv.Atoi(string(line[sidStart:sidEnd]))
		if err != nil {
			panic("Error parsing read seq SI")
		}

		geneIndex := (sidInt / 2)
		readStartInt, err := strconv.Atoi(string(readStart))
		if err != nil {
			panic("Error parsing read start pos from sam file")
		}

		gapIntervals, err := GetGapPositions(readStartInt, cigar)
		if err != nil {
			fmt.Println(err)
		}

		for _, gap := range gapIntervals {
			gaps[geneIndex] = append(gaps[geneIndex], gap)
		}
	}

	return gaps
}

type CigarOp struct {
	Length int
	Op     byte
}

func ParseCIGAR(cigar []byte) ([]CigarOp, error) {

	var ops []CigarOp

	start := 0
	for i := 0; i < len(cigar); i++ {
		// as soon as we hit a non-numeric byte, we can create a new cigarOp
		if cigar[i] < '0' || cigar[i] > '9' {
			cigraPrefix, err := strconv.Atoi(string(cigar[start:i]))
			if err != nil {
				return nil, err
			}
			ops = append(ops, CigarOp{Op: cigar[i], Length: cigraPrefix})
			start = i + 1
		}
	}

	return ops, nil
}

type SAMIndexEntry struct {
	Name           []byte
	Offset         int64
	FirstOfPair    bool
	DelimPositions [11]uint16 // this is hardcoded since we know how many fields each sam entry has. Using this array, we can quickly jump to a desired field
}

func CreateSAMIndex(samfile *os.File) map[string]*SAMReadPair {

	scanner := bufio.NewScanner(samfile)
	const maxCapacity = 1024 * 1024 // 1MB
	buf := make([]byte, maxCapacity)
	scanner.Buffer(buf, maxCapacity)

	samIndex := make(map[string]*SAMReadPair)

	offset := int64(0)
	for scanner.Scan() {
		line := scanner.Bytes()
		if line[0] == '@' {
			offset += int64(len(line) + 1)
			continue
		}

		name, delimPositions, firstOfPair := parseDelims(line)
		indexEntry := &SAMIndexEntry{
			Name:           name,
			Offset:         offset,
			FirstOfPair:    firstOfPair,
			DelimPositions: delimPositions,
		}
		strName := string(name)
		_, exists := samIndex[strName]
		if !exists {
			samIndex[strName] = &SAMReadPair{}
		}

		if firstOfPair {
			samIndex[strName].first = indexEntry
		} else {
			samIndex[strName].second = indexEntry
		}

		offset += int64(len(line)) + 1 // +1 for newline
	}

	return samIndex
}

type SAMReadPair struct {
	first  *SAMIndexEntry
	second *SAMIndexEntry
}

func parseDelims(line []byte) ([]byte, [11]uint16, bool) {
	tabCount := 0
	delimPositions := [11]uint16{}
	var name []byte
	var flagStart int
	var isFist bool
	for i := 0; i < len(line); i++ {
		if line[i] == '\t' {

			if tabCount == 0 {
				name = make([]byte, i)
				copy(name, line[0:i])
				flagStart = i + 1
			}
			if tabCount == 1 {
				flagBytes := line[flagStart:i]
				flagInt, err := strconv.Atoi(string(flagBytes))
				if err != nil {
					panic("Error parsing flag of read")
				}
				isFist = flagInt&0x40 != 0
			}
			delimPositions[tabCount] = uint16(i)
			tabCount++
		}

	}
	return name, delimPositions, isFist
}

func GetGapPositions(start int, cigar []byte) ([]*interval.Interval, error) {
	ops, err := ParseCIGAR(cigar)
	if err != nil {
		return nil, err
	}

	refPos := start
	var gaps []*interval.Interval

	for _, op := range ops {
		switch op.Op {
		case 'M', 'D', 'N', '=', 'X':
			if op.Op == 'N' {
				gaps = append(gaps, &interval.Interval{Start: refPos, End: refPos + op.Length - 1})
			}
			refPos += op.Length
		}
	}

	return gaps, nil
}

type IntronCluster struct {
	lStart         int // left most coord
	rStop          int // right most coord
	maxEvidence    int // how many gaps had eStart and eStop
	eStart         int // start of gap with max evidence
	eStop          int // stop of gap with max evidence
	inclusionReads int // a count which stores how many reads are contained within the gap. This helps to identify I. potential deletions but also II. wrong gaps
}

func clusterGaps(gapMapPerSeq map[int]map[interval.Interval]int) map[int][]*IntronCluster {
	// cluster overlapping gaps
	clusters := make(map[int][]*IntronCluster, 0)
	for sid, gapMap := range gapMapPerSeq {
		if _, ok := clusters[sid]; !ok {
			clusters[sid] = []*IntronCluster{}
		}
		for currGap, currGapEvidence := range gapMap {
			if len(clusters[sid]) == 0 {
				// init slice with first cluster
				cluster := &IntronCluster{
					lStart:         currGap.Start,
					rStop:          currGap.End,
					maxEvidence:    currGapEvidence,
					eStart:         currGap.Start,
					eStop:          currGap.End,
					inclusionReads: 0, // init with 0
				}
				clusters[sid] = append(clusters[sid], cluster)
				continue
			}
			// keep track if curr gap was added
			addedToExistingCluster := false

			for _, cluster := range clusters[sid] {
				// is the currGap overlapping with the cluster
				// I. -----##########------- cluster
				//    ------######---------- currGap → contained by cluster, add to cluster
				if currGap.Start >= cluster.lStart && currGap.End <= cluster.rStop {
					if cluster.maxEvidence < currGapEvidence {
						cluster.maxEvidence = currGapEvidence
						cluster.eStart = currGap.Start
						cluster.eStop = currGap.End
					}
					addedToExistingCluster = true
					break
				}
				// II -----##########------- cluster
				//    ---#############------ currGap → contains cluster completely
				if currGap.Start <= cluster.lStart && currGap.End >= cluster.rStop {
					if cluster.maxEvidence < currGapEvidence {
						cluster.maxEvidence = currGapEvidence
						cluster.eStart = currGap.Start
						cluster.eStop = currGap.End
						cluster.lStart = currGap.Start
						cluster.rStop = currGap.End
					} else {
						cluster.lStart = currGap.Start
						cluster.rStop = currGap.End
					}
					addedToExistingCluster = true
					break
				}
				// III -----##########------- cluster
				//     ---############------  currGap → partially contains cluster completely
				if currGap.Start <= cluster.lStart && currGap.End >= cluster.lStart && currGap.End <= cluster.rStop {
					if cluster.maxEvidence < currGapEvidence {
						cluster.maxEvidence = currGapEvidence
						cluster.eStart = currGap.Start
						cluster.eStop = currGap.End
					}
					cluster.lStart = currGap.Start
					addedToExistingCluster = true
					break
				}
				// IV  -----##########------- cluster
				//     -----############----  currGap → partially contains cluster completely
				if currGap.Start >= cluster.lStart && currGap.Start <= cluster.rStop && currGap.End >= cluster.rStop {
					if cluster.maxEvidence < currGapEvidence {
						cluster.maxEvidence = currGapEvidence
						cluster.eStart = currGap.Start
						cluster.eStop = currGap.End
					}
					cluster.rStop = currGap.End
					addedToExistingCluster = true
					break
				}
			}

			// if the gap wasn't added to any existing cluster, create a new one
			if !addedToExistingCluster {
				newCluster := &IntronCluster{
					lStart:         currGap.Start,
					rStop:          currGap.End,
					maxEvidence:    currGapEvidence,
					eStart:         currGap.Start,
					eStop:          currGap.End,
					inclusionReads: 0, // init with 0
				}
				clusters[sid] = append(clusters[sid], newCluster)
			}
		}
	}
	return clusters
}

// ParseSparseSAMEntry takes a file and a sparseSAMEntry and parses the corresponding sam entry
func ParseSparseSAMEntry(file *os.File, sparseSAMEntry *SAMIndexEntry) (*mapperutils.ReadMatchResult, error) {
	if sparseSAMEntry == nil || (sparseSAMEntry.Offset < 0) {
		return nil, errors.New("invalid read: no offset available")
	}
	offset := sparseSAMEntry.Offset

	// jump to read pos in sam
	_, err := file.Seek(offset, io.SeekStart)
	if err != nil {
		return nil, err
	}

	reader := bufio.NewReader(file)
	line, err := reader.ReadBytes('\n')
	if err != nil && err != io.EOF {
		return nil, err
	}

	// remove trailing newline if present
	if len(line) > 0 && line[len(line)-1] == '\n' {
		line = line[:len(line)-1]
	}
	sidBytes := line[sparseSAMEntry.DelimPositions[10]+uint16(6) : len(line)]
	startBytes := line[sparseSAMEntry.DelimPositions[2]+1 : sparseSAMEntry.DelimPositions[3]]
	//qualityRange := line[sparseSAMEntry.DelimPositions[3]+1 : sparseSAMEntry.DelimPositions[4]]
	cigarBytes := line[sparseSAMEntry.DelimPositions[4]+1 : sparseSAMEntry.DelimPositions[5]]
	sequenceBytes := line[sparseSAMEntry.DelimPositions[8]+1 : sparseSAMEntry.DelimPositions[9]]
	//qualityStrBytes := line[sparseSAMEntry.DelimPositions[9]+1 : sparseSAMEntry.DelimPositions[10]]

	sid, err := strconv.Atoi(string(sidBytes))
	if err != nil {
		return nil, err
	}
	startPos, err := strconv.Atoi(string(startBytes))
	if err != nil {
		return nil, err
	}

	result := &mapperutils.ReadMatchResult{
		SequenceIndex: sid,
		SecondPass:    false, // this is always false since the read is listed in the output sam
	}

	result.MatchedRead = regionvector.NewRegionVector()
	result.MatchedGenome = regionvector.NewRegionVector()

	// Process CIGAR string to build region vectors and identify mismatches
	err = processCigarForRegions(cigarBytes, startPos, sequenceBytes, result)
	if err != nil {
		return nil, err
	}

	return result, nil
}

// processCigarForRegions analyzes the CIGAR string to build read and genome region vectors
func processCigarForRegions(cigar []byte, refPos int, readSeq []byte, result *mapperutils.ReadMatchResult) error {
	readPos := 0

	// get cigar ops
	cigarOps, err := ParseCIGAR(cigar)
	if err != nil {
		return err
	}

	for _, op := range cigarOps {
		switch op.Op {
		case 'M', '=', 'X':
			result.MatchedRead.AddRegionNonOverlappingPanic(readPos, readPos+op.Length)
			result.MatchedGenome.AddRegionNonOverlappingPanic(refPos, refPos+op.Length)

			if op.Op == 'X' {
				for i := 0; i < op.Length; i++ {
					result.MismatchesRead = append(result.MismatchesRead, readPos+i)
				}
			}

			readPos += op.Length
			refPos += op.Length

		case 'I':
			readPos += op.Length

		case 'D', 'N':
			refPos += op.Length

		case 'S':
			readPos += op.Length

		case 'H':
			// nothing is consumed

		case 'P':
			// nothing is consumed
		}
	}

	return nil
}
