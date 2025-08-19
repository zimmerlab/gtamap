package server

import (
	"fmt"

	"github.com/KleinSamuel/gtamap/src/core/datastructure/regionvector"
	"github.com/KleinSamuel/gtamap/src/formats/cigar"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
	"github.com/sirupsen/logrus"
)

type EnhancedRecord struct {
	sam.Record
	Index        int                        // id to identify the record = count of global records
	MapperIndex  int                        // index of the mapper in MapperInfos
	MappedGenome *regionvector.RegionVector // the matched intervals in the genome (mis- + matches (M, =, X in cigar))
	MappedRead   *regionvector.RegionVector // the matched intervals in the read (matches (M, =, X) in cigar)
	Mismatches   []int                      // the positions of the mismatches in the read sequence (min = 0, max = read length)
	CigarObj     *cigar.Object
	IsAccepted   bool
	QnameCluster *QnameCluster // pointer to the cluster of records with the same qname
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
		QnameCluster: nil,
	}, nil
}

func (r *EnhancedRecord) IsEqual(other *EnhancedRecord) bool {
	if r == nil || other == nil {
		return false
	}
	if r.Qname != other.Qname || r.Rname != other.Rname || r.Pos != other.Pos {
		return false
	}
	if r.CigarObj.String() != other.CigarObj.String() {
		return false
	}
	return true
}

func (r *EnhancedRecord) StringFinal(posOffset int, mapperNames string, recordIndices string) string {

	alignmentLineString := fmt.Sprintf("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
		r.Qname, r.Flag.String(), r.Rname, r.Pos-posOffset, r.Mapq, r.Cigar,
		r.Rnext, r.Pnext, r.Tlen, r.Seq, r.Qual)

	alignmentLineString += fmt.Sprintf("\tXG:Z:%s", recordIndices)

	alignmentLineString += fmt.Sprintf("\tXM:Z:%s", mapperNames)

	return alignmentLineString
}
