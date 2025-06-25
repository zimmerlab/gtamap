package analysis

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
)

func CompareResults() {

	gtamapTarget := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.gtamap.target.sam"
	minimap2Genome := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.minimap2.genome.sam"

	gtamapTargetFile, errGtamapTarget := sam.ParseSAMFile(gtamapTarget)
	if errGtamapTarget != nil {
		panic("Error reading GTAMap target file: " + errGtamapTarget.Error())
	}

	gtamapRecordsByQname := make(map[string][]*sam.Record)
	for _, record := range gtamapTargetFile.Records {
		if _, exists := gtamapRecordsByQname[record.Qname]; !exists {
			gtamapRecordsByQname[record.Qname] = []*sam.Record{}
		}
		gtamapRecordsByQname[record.Qname] = append(gtamapRecordsByQname[record.Qname], &record)
	}

	minimap2GenomeFile, errMinimap2Genome := sam.ParseSAMFile(minimap2Genome)
	if errMinimap2Genome != nil {
		panic("Error reading Minimap2 genome file: " + errMinimap2Genome.Error())
	}

	minimap2RecordsByQname := make(map[string][]*sam.Record)
	for _, record := range minimap2GenomeFile.Records {
		if _, exists := minimap2RecordsByQname[record.Qname]; !exists {
			minimap2RecordsByQname[record.Qname] = []*sam.Record{}
		}
		minimap2RecordsByQname[record.Qname] = append(minimap2RecordsByQname[record.Qname], &record)
	}

	fmt.Printf("gtamap qnames:\t\t%d\n", len(gtamapRecordsByQname))
	fmt.Printf("gtamap records:\t\t%d\n", len(gtamapTargetFile.Records))
	fmt.Printf("minimap2 qnames:\t%d\n", len(minimap2RecordsByQname))
	fmt.Printf("minimap2 records:\t%d\n", len(minimap2GenomeFile.Records))

	// +1 if qname in gtamap and in minimap2
	countSameQname := 0
	// +1 if any position of the same qname is found in minimap2
	countSameQnameSamePos := 0
	// +1 if any position and cigar of the same qname is found in minimap2
	countSameQnameSamePosSameCigar := 0
	// +1 if record in gtamap is found in minimap2 with same pos
	countSameRecordSamePos := 0
	// +1 if record in gtamap is found in minimap2 with same pos and same qname
	countSameRecordSamePosSameCigar := 0

	for qname, gtamapRecords := range gtamapRecordsByQname {

		if _, exists := minimap2RecordsByQname[qname]; !exists {
			//fmt.Printf("GTAMap record %s not found in Minimap2 genome\n", record.Qname)
			continue
		}

		minimap2Records := minimap2RecordsByQname[qname]
		if len(minimap2Records) == 0 {
			//fmt.Printf("No records found for Qname %s in Minimap2 genome\n", record.Qname)
			continue
		}

		countSameQname++

		foundAnySamePos := false
		foundAnySamePosSameCigar := false

		for _, gtamapRecord := range gtamapRecords {

			for _, minimap2Record := range minimap2Records {

				if minimap2Record.Rname == gtamapRecord.Rname && minimap2Record.Pos == gtamapRecord.Pos {
					countSameRecordSamePos++
					foundAnySamePos = true

					if minimap2Record.Cigar == gtamapRecord.Cigar {
						countSameRecordSamePosSameCigar++
						foundAnySamePosSameCigar = true
					}

					break
				}
			}
		}

		if foundAnySamePos {
			countSameQnameSamePos++
		}
		if foundAnySamePosSameCigar {
			countSameQnameSamePosSameCigar++
		}
	}

	fmt.Printf("minimap2 found %d of %d qnames\n", countSameQname, len(gtamapRecordsByQname))
	fmt.Printf("minimap2 found %d of %d same qname and same pos\n", countSameQnameSamePos, len(gtamapRecordsByQname))
	fmt.Printf("minimap2 found %d of %d same qname, same pos and same cigar\n", countSameQnameSamePosSameCigar, len(gtamapRecordsByQname))
	fmt.Printf("minimap2 found %d of %d same records with same pos\n", countSameRecordSamePos, len(gtamapTargetFile.Records))
	fmt.Printf("minimap2 found %d of %d same records with same pos and same cigar\n", countSameRecordSamePosSameCigar, len(gtamapTargetFile.Records))
}
