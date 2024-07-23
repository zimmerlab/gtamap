package gtf

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

type Annotation struct {
	Genes []*Gene
}

type Gene struct {
	GeneIdEnsembl       string
	Chromosome          string
	IsForwardStrand     bool
	StartGenomic        uint32
	EndGenomic          uint32
	Transcripts         []*Transcript
	ConsensusTranscript *Transcript
}

type Transcript struct {
	TranscriptIdEnsembl string
	// start position of the transcript relative to the genomic start of its parent gene
	StartRelative uint32
	// end position of the transcript relative to the genomic start of its parent gene
	EndRelative uint32
	SequenceDna string
	Exons       []*Exon
}

type Exon struct {
	StartRelative uint32
	EndRelative   uint32
}

func ReadGtfFromFile(filePath string) *Annotation {

	file, err := os.Open(filePath)
	if err != nil {
		return nil
	}
	defer file.Close()

	return ReadGtf(file)
}

func ReadGtf(gtfFile *os.File) *Annotation {

	annot := Annotation{
		Genes: make([]*Gene, 0),
	}

	// buffered reader
	scanner := bufio.NewScanner(gtfFile)

	var currentGene *Gene
	var currentTranscript *Transcript

	for scanner.Scan() {
		line := strings.Split(scanner.Text(), "\t")

		seqName := &line[0]
		feature := &line[2]
		start64, errStart := strconv.ParseUint(line[3], 10, 32)
		end64, errEnd := strconv.ParseUint(line[4], 10, 32)
		if errStart != nil && errEnd != nil {
			fmt.Println("Could not parse start or end position: ", line[3], line[4])
		}
		// the gtf file is 1-based, but the annotation is 0-based
		start := uint32(start64) - 1
		end := uint32(end64)
		isForwardStrand := line[6] == "+"

		switch *feature {
		case "gene":

			gene := Gene{
				GeneIdEnsembl:   extractAttributeValue(line[8], "gene_id"),
				Chromosome:      *seqName,
				IsForwardStrand: isForwardStrand,
				StartGenomic:    start,
				EndGenomic:      end,
				Transcripts:     make([]*Transcript, 0),
			}

			currentGene = &gene
			annot.Genes = append(annot.Genes, &gene)

		case "transcript":

			transcript := Transcript{
				TranscriptIdEnsembl: extractAttributeValue(line[8], "transcript_id"),
				StartRelative:       start - currentGene.StartGenomic,
				EndRelative:         end - currentGene.StartGenomic,
				// the sequence is built later based on the exon structure and the supplied fasta file
				SequenceDna: "",
			}

			currentTranscript = &transcript
			currentGene.Transcripts = append(currentGene.Transcripts, &transcript)
		case "cds":
			exon := Exon{
				StartRelative: start - currentGene.StartGenomic,
				EndRelative:   end - currentGene.StartGenomic,
			}
			currentTranscript.Exons = append(currentTranscript.Exons, &exon)
		}
	}

	if err := scanner.Err(); err != nil {
		return nil
	}

	return &annot
}

func extractAttributeValue(attributeString string, attributeKey string) string {
	for _, item := range strings.Split(attributeString, ";") {
		split := strings.Split(strings.TrimSpace(item), " ")
		if strings.TrimSpace(split[0]) == attributeKey {
			return strings.Replace(strings.TrimSpace(split[1]), "\"", "", -1)
		}
	}
	return ""
}
