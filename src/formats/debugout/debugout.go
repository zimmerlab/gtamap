package debugout

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"os"
	"strconv"
	"strings"
)

const contextWindow = 5

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
	file.WriteString(view)
	file.WriteString("\n")
	file.WriteString("\n")
}

func getAlignment(genomeIndex *index.GenomeIndex, mappedRead mapperutils.ReadMatchResult, readMeta *fastq.Read) string {
	var readView strings.Builder
	var mapView strings.Builder
	var geneView strings.Builder
	var coordView strings.Builder
	var builder strings.Builder
	var geneSeq = genomeIndex.Sequences[mappedRead.SequenceIndex]
	// merge regions based on genomic region vector:
	// read := {0, 10}, {10, 120}, {120, 150}}
	// gene := {100, 110}, {110, 230}, {1020, 1050}}
	// mergedRead: [{0 120}, {120 150}]
	// mergedGene: [{100 230}, {1020 1050}]
	mergedReadRegions, mergedGeneRegions := mergeAlignedRegions(mappedRead.MatchedRead.Regions, mappedRead.MatchedGenome.Regions)

	for i := 0; i < len(mergedReadRegions); i++ {
		readRegion := mergedReadRegions[i]
		geneRegion := mergedGeneRegions[i]
		readStart := readRegion.Start
		readStop := readRegion.End
		geneStart := geneRegion.Start
		geneStop := geneRegion.End

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
		strStart := strconv.Itoa(readStart)
		strStartInGene := strconv.Itoa(geneStart)
		strStop := strconv.Itoa(readStop)
		strStopInGene := strconv.Itoa(geneStop)
		coordLength := len(strStart) + len(strStop) + len(strStartInGene) + len(strStopInGene) + 4 // 4 chars for sep
		for j := geneStart; j < geneStop; j++ {
			if j == geneStart {
				coordView.WriteByte('|')
				coordView.WriteString(strStart)
				coordView.WriteByte('-')
				coordView.WriteString(strStop)
				coordView.WriteByte(',')
				coordView.WriteString(strStartInGene)
				coordView.WriteByte('-')
				coordView.WriteString(strStopInGene)
			} else if count > coordLength {
				coordView.WriteByte(' ')
			}

			posInRead := readStart + j - geneStart
			readView.WriteByte((*readMeta.Sequence)[posInRead])
			mapView.WriteByte(getMappingByte((*readMeta.Sequence)[readStart+j-geneStart], (*geneSeq)[j]))
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

		// TODO: this would currently also append ----- to the actual end of a gene
		if i == len(mappedRead.MatchedRead.Regions)-1 {
			// append suffix separator
			for j := geneStop - contextWindow; j < geneStop; j++ {
				// readView.WriteByte(' ')
				// coordView.WriteByte(' ')
				// mapView.WriteByte(' ')
				geneView.WriteByte('-')
			}
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
	if len(readRegions) != len(genomeRegions) {
		panic("read and genome regions must be aligned and equal in length")
	}
	if len(readRegions) == 0 {
		return nil, nil
	}

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
