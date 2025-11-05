package main

import (
	"bufio"
	"os"
	"strconv"
	"strings"

	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
)

func main() {
	logrus.SetFormatter(&logrus.TextFormatter{
		FullTimestamp: true,
	})
	logrus.SetLevel(logrus.InfoLevel)

	mappingInfoPath := "/mnt/raidproj/proj/projekte/mapping/gtamap/mapper-diff/dmd-err0/read.mappinginfo"
	fastqR1Path := "/mnt/raidproj/proj/projekte/mapping/gtamap/mapper-diff/dmd-err0/fw.fastq"
	fastqR2Path := "/mnt/raidproj/proj/projekte/mapping/gtamap/mapper-diff/dmd-err0/rw.fastq"

	samPath := "/mnt/raidproj/proj/projekte/mapping/gtamap/mapper-diff/dmd-err0/read.sam"

	fastqReader, err := fastq.InitFromPaths(&fastqR1Path, &fastqR2Path)
	if err != nil {
		logrus.Fatalf("Error initializing FASTQ reader: %v", err)
	}

	// Create SAM output file
	samFile, err := os.Create(samPath)
	if err != nil {
		logrus.Fatalf("Error creating SAM file: %v", err)
	}
	defer samFile.Close()

	// Write SAM header
	_, err = samFile.WriteString("@HD\tVN:1.6\n")
	if err != nil {
		logrus.Fatalf("Error writing SAM header: %v", err)
	}
	_, err = samFile.WriteString("@SQ\tSN:X\tLN:33339609\n")
	if err != nil {
		logrus.Fatalf("Error writing SAM header: %v", err)
	}

	file, err := os.Open(mappingInfoPath)
	if err != nil {
		logrus.Fatalf("Error opening mappinginfo file: %v", err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	lineNumber := 0

	// Skip first line of mappinginfo
	if scanner.Scan() {
		lineNumber++
	}

	// Read FASTQ entries and mappinginfo simultaneously
	for {
		// Read FASTQ entry (4 lines: header, sequence, +, quality)
		readPair, _ := fastqReader.NextRead()
		if readPair == nil {
			break
		}

		// Read mappinginfo entry (1 line)
		if !scanner.Scan() {
			break
		}
		lineNumber++

		line := scanner.Text()
		fields := strings.Split(line, "\t")

		if len(fields) < 4 {
			logrus.Warnf("Incomplete mappinginfo line %d: %s", lineNumber, line)
			continue
		}

		qname := fields[0]
		contig := fields[1]
		fw_range := fields[4]
		rv_range := fields[5]

		fw_rv := rawRangeToRegionVector(fw_range)
		rv_rv := rawRangeToRegionVector(rv_range)

		samR1, samR2 := toSamRecordStrings(qname, contig, fw_rv, rv_rv, readPair.ReadR1, readPair.ReadR2)

		// Write SAM records to file
		_, err = samFile.WriteString(samR1 + "\n")
		if err != nil {
			logrus.Fatalf("Error writing SAM R1 record: %v", err)
		}
		_, err = samFile.WriteString(samR2 + "\n")
		if err != nil {
			logrus.Fatalf("Error writing SAM R2 record: %v", err)
		}

		// fmt.Println(samR1)
		// fmt.Println(samR2)

		// break
	}

	if err := scanner.Err(); err != nil {
		logrus.Fatalf("Error reading mappinginfo file: %v", err)
	}
}

func toSamRecordStrings(
	qname string,
	contig string,
	fwRv *regionvector.RegionVector,
	rvRv *regionvector.RegionVector,
	readR1 *fastq.Read,
	readR2 *fastq.Read,
) (string, string) {
	sbR1 := strings.Builder{}
	sbR2 := strings.Builder{}

	flagR1 := sam.Flag{}
	flagR1.SetPaired()
	flagR1.SetProperlyPaired()
	flagR1.SetFirstInPair()
	flagR1.SetMateReverseStrand()

	fwStart := strconv.Itoa(fwRv.Regions[0].Start)

	flagR2 := sam.Flag{}
	flagR2.SetPaired()
	flagR2.SetProperlyPaired()
	flagR2.SetSecondInPair()
	flagR2.SetReverseStrand()

	rvStart := strconv.Itoa(rvRv.Regions[0].Start)

	sbR1.WriteString(qname)
	sbR1.WriteString("\t")
	sbR1.WriteString(flagR1.String())
	sbR1.WriteString("\t")
	sbR1.WriteString(contig)
	sbR1.WriteString("\t")
	sbR1.WriteString(fwStart)
	sbR1.WriteString("\t254\t")
	sbR1.WriteString("100M")
	sbR1.WriteString("\t")
	sbR1.WriteString("*")
	sbR1.WriteString("\t")
	sbR1.WriteString(rvStart)
	sbR1.WriteString("\t")
	sbR1.WriteString("0")
	sbR1.WriteString("\t")
	sbR1.WriteString(string(*readR1.Sequence))
	sbR1.WriteString("\t")
	sbR1.WriteString(string(*readR1.Quality))

	sbR2.WriteString(qname)
	sbR2.WriteString("\t")
	sbR2.WriteString(flagR2.String())
	sbR2.WriteString("\t")
	sbR2.WriteString(contig)
	sbR2.WriteString("\t")
	sbR2.WriteString(rvStart)
	sbR2.WriteString("\t254\t")
	sbR2.WriteString("100M")
	sbR2.WriteString("\t")
	sbR2.WriteString("*")
	sbR2.WriteString("\t")
	sbR2.WriteString(fwStart)
	sbR2.WriteString("\t")
	sbR2.WriteString("0")
	sbR2.WriteString("\t")
	sbR2.WriteString(utils.ReverseComplementDNA(string(*readR2.Sequence)))
	sbR2.WriteString("\t")
	sbR2.WriteString(utils.ReverseString(string(*readR2.Quality)))

	return sbR1.String(), sbR2.String()
}

func rawRangeToRegionVector(rawRange string) *regionvector.RegionVector {
	rv := regionvector.NewRegionVector()

	fields := strings.Split(rawRange, "|")

	for _, field := range fields {

		startStop := strings.Split(field, "-")

		start, errStart := strconv.Atoi(startStop[0])
		if errStart != nil {
			logrus.Fatalf("Error converting start to int: %v", errStart)
		}

		end, errEnd := strconv.Atoi(startStop[1])
		if errEnd != nil {
			logrus.Fatalf("Error converting end to int: %v", errEnd)
		}

		rv.AddRegion(start, end)
	}

	return rv
}
