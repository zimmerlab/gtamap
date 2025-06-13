package debugout

import (
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
)

const contextWindow = 10

func GenerateAlignmentView(genomeIndex *index.GenomeIndex, mappedRead mapperutils.ReadMatchResult, readMeta *fastq.Read) {
	view := getAlignment(genomeIndex, mappedRead, readMeta)
	// append to file
	file, err := os.OpenFile("output.txt", os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0o644)
	if err != nil {
		fmt.Println("Error creating file:", err)
	}
	defer file.Close()

	file.WriteString("Read ID: ")
	file.WriteString(readMeta.Header)
	file.WriteString("\n")
	file.WriteString("Sequence Index: ")
	file.WriteString(strconv.Itoa(mappedRead.SequenceIndex))
	file.WriteString("\n")
	file.WriteString("Splice Sites (before second pass): ")
	for _, s := range mappedRead.SpliceSitesInfo {
		if s {
			file.WriteString("True ")
		} else {
			file.WriteString("False ")
		}
	}
	file.WriteString("\n")
	file.WriteString(view)
	file.WriteString(string(*genomeIndex.Sequences[mappedRead.SequenceIndex]))
	file.WriteString("\n")
	file.WriteString("\n")
}

func getAlignment(genomeIndex *index.GenomeIndex, mappedRead mapperutils.ReadMatchResult, readMeta *fastq.Read) string {
	var readView strings.Builder
	var mapView strings.Builder
	var geneView strings.Builder
	var coordView strings.Builder
	var builder strings.Builder
	geneSeq := genomeIndex.Sequences[mappedRead.SequenceIndex]
	readStartALignment := mappedRead.MatchedRead.Regions[0].Start

	for _, genomeRegion := range mappedRead.MatchedGenome.GetAlignmentBlocks() {
		geneStart := genomeRegion.Start
		geneStop := genomeRegion.End
		readStartCoord := regionvector.GenomicCoordToReadCoord(readStartALignment, geneStart, mappedRead.MatchedGenome.Regions)
		readStopCoord := regionvector.GenomicCoordToReadCoord(readStartALignment, geneStop, mappedRead.MatchedGenome.Regions)

		// append prefix separator
		for j := geneStart - contextWindow; j < geneStart; j++ {
			readView.WriteByte(' ')
			coordView.WriteByte(' ')
			mapView.WriteByte(' ')
			geneView.WriteByte('-')
		}

		// append prefix only if readMapping starts past contextWindow
		if geneStart > contextWindow {
			for j := geneStart - contextWindow; j < geneStart; j++ {
				readView.WriteByte(' ')
				mapView.WriteByte(' ')
				coordView.WriteByte(' ')
				geneView.WriteByte((*geneSeq)[j])
			}
		}

		// append alignment
		count := 0
		strStart := strconv.Itoa(readStartCoord)
		strStartInGene := strconv.Itoa(geneStart)
		strStop := strconv.Itoa(readStopCoord)
		strStopInGene := strconv.Itoa(geneStop)
		coordLength := len(strStart) + len(strStartInGene) + len(strStop) + len(strStopInGene) + 6 // 6 chars for sep

		for j := genomeRegion.Start; j < genomeRegion.End; j++ {
			readPos := regionvector.GenomicCoordToReadCoord(readStartALignment, j, mappedRead.MatchedGenome.Regions)

			if j == geneStart {
				coordView.WriteByte('|')
				coordView.WriteString(strStart)
				coordView.WriteByte('-')
				coordView.WriteString(strStop)
				coordView.WriteByte(';')
				coordView.WriteString(strStartInGene)
				coordView.WriteByte('-')
				coordView.WriteString(strStopInGene)
			} else if count > coordLength {
				coordView.WriteByte(' ')
			}

			readView.WriteByte((*readMeta.Sequence)[readPos])
			mapView.WriteByte(getMappingByte((*readMeta.Sequence)[readPos], (*geneSeq)[j]))
			geneView.WriteByte((*geneSeq)[j])
			count++
		}

		// append suffix
		if geneStop+contextWindow < len(*geneSeq) {
			for j := geneStop; j < geneStop+contextWindow; j++ {
				readView.WriteByte(' ')
				mapView.WriteByte(' ')
				coordView.WriteByte(' ')
				geneView.WriteByte((*geneSeq)[j])
			}
		}
		for i := 0; i < 3; i++ {
			coordView.WriteByte(' ')
		}
	}
	builder.WriteString(coordView.String())
	builder.WriteByte('\n')
	builder.WriteString(readView.String())
	builder.WriteByte('\n')
	builder.WriteString(mapView.String())
	builder.WriteByte('\n')
	builder.WriteString(geneView.String())
	builder.WriteByte('\n')

	return builder.String()
}

func getMappingByte(a byte, b byte) byte {
	if a == b {
		return '.'
	}
	return '|'
}

func mergeRegions(regions []*regionvector.Region) []regionvector.Region {
	if len(regions) == 0 {
		return nil
	}

	merged := []regionvector.Region{*regions[0]}

	for _, curr := range regions[1:] {
		last := merged[len(merged)-1]

		// If current region overlaps or is adjacent to the last one, merge them
		if curr.Start <= last.End+1 {
			if curr.End > last.End {
				last.End = curr.End
			}
		} else {
			merged = append(merged, *curr)
		}
	}
	return merged
}

func mergeAlignedRegions(readRegions []*regionvector.Region, genomeRegions []*regionvector.Region) ([]regionvector.Region, []regionvector.Region) {
	mergedReads := []regionvector.Region{*readRegions[0]}
	mergedGenomes := []regionvector.Region{*genomeRegions[0]}

	for i := 1; i < len(readRegions); i++ {
		lastGenome := &mergedGenomes[len(mergedGenomes)-1]
		currGenome := genomeRegions[i]

		if currGenome.Start <= lastGenome.End+1 {
			// Merge genome region
			if currGenome.End > lastGenome.End {
				lastGenome.End = currGenome.End
			}

			// Merge corresponding read region
			lastRead := &mergedReads[len(mergedReads)-1]
			if readRegions[i].End > lastRead.End {
				lastRead.End = readRegions[i].End
			}
		} else {
			mergedGenomes = append(mergedGenomes, *genomeRegions[i])
			mergedReads = append(mergedReads, *readRegions[i])
		}
	}

	return mergedReads, mergedGenomes
}
