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
	GeneIdEnsembl   string
	Chromosome      string
	IsForwardStrand bool
	StartGenomic    uint32
	EndGenomic      uint32
	Transcripts     []*Transcript
}

type Transcript struct {
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
		start := uint32(start64)
		end := uint32(end64)
		isForwardStrand := line[6] == "+"

		switch *feature {
		case "gene":
			gene := Gene{
				GeneIdEnsembl:   "",
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
				StartRelative: start - currentGene.StartGenomic,
				EndRelative:   end - currentGene.StartGenomic,
				SequenceDna:   "",
			}
			currentTranscript = &transcript
			currentGene.Transcripts = append(currentGene.Transcripts, &transcript)
		case "exon":
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
