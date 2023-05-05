package dataloader

import (
	"github.com/KleinSamuel/gtamap/src/fasta"
	"github.com/KleinSamuel/gtamap/src/gtf"
	"os"
	"sort"
	"strings"
)

func panicIfFilesAreMissing(pathGtf string, pathFasta string, pathFastaIndex string) {

	// checks if the fasta index file exists
	_, errGtf := os.Stat(pathGtf)
	if os.IsNotExist(errGtf) {
		panic("gtf file does not exist: " + pathGtf)
	}
	// check if the fasta index file exists
	_, errFasta := os.Stat(pathFasta)
	if os.IsNotExist(errFasta) {
		panic("fasta index file does not exist: " + pathFasta)
	}
	// check if the fasta index file exists
	_, errIndex := os.Stat(pathFastaIndex)
	if os.IsNotExist(errIndex) {
		panic("fasta index file does not exist: " + pathFastaIndex)
	}
}

// GenerateInputForIndex reads in the genome annotation from given gtf file,
// reads the fasta index (required) and extracts the dna sequence for every transcript
func GenerateInputForIndex(pathGtf string, pathFasta string) *gtf.Annotation {

	pathFastaIndex := pathFasta + ".fai"

	panicIfFilesAreMissing(pathGtf, pathFasta, pathFastaIndex)

	fastaIndex := fasta.ReadFastaIndex(pathFastaIndex)
	annotation := gtf.ReadGtfFromFile(pathGtf)

	fastaFile, err := os.Open(pathFasta)
	if err != nil {
		panic(err)
	}
	defer fastaFile.Close()

	for _, gene := range annotation.Genes {

		index := fastaIndex.Entries[gene.Chromosome]

		for _, trans := range gene.Transcripts {

			// sort exons by start position such that transcript sequence is in order
			sort.Slice(trans.Exons, func(i, j int) bool {
				return trans.Exons[i].StartRelative < trans.Exons[j].StartRelative
			})

			var transcriptSequence []byte

			// build transcript sequence by seeking and concatenating the exon sequences
			for _, exon := range trans.Exons {

				// the start position of the transcript relative to its chromosome
				startRelativeChromosome := gene.StartGenomic - 1 + exon.StartRelative
				// the number of newlines contained between the start of the chromosome and the start of the transcript
				offsetNewlines := int64(startRelativeChromosome) / int64(index.Linebases)
				// the offset to the start of the transcript
				offset := index.Offset + int64(startRelativeChromosome) + offsetNewlines

				_, errSeek := fastaFile.Seek(offset, 0)
				if errSeek != nil {
					panic(errSeek)
				}

				// the length of the exon sequence in nucleotides
				lengthBases := exon.EndRelative - exon.StartRelative + 1
				// the number of newlines contained between the start and end of the exon sequence
				lengthNewlines := int64(lengthBases) / int64(index.Linebases)
				// the number of characters in the current line until the next newline
				numCharsInLine := int64(index.Linebases) - (int64(startRelativeChromosome) % int64(index.Linebases))
				// adds a newline if the exon sequence spans across the current line
				if int64(lengthBases) > numCharsInLine {
					lengthNewlines++
				}
				// the total length of the exon sequence in bytes (including newlines)
				length := int64(lengthBases) + lengthNewlines

				buffer := make([]byte, length)

				_, errRead := fastaFile.Read(buffer)
				if errRead != nil {
					panic(errRead)
				}

				transcriptSequence = append(transcriptSequence, buffer...)
			}

			trans.SequenceDna = strings.ReplaceAll(string(transcriptSequence), "\n", "")
		}
	}

	return annotation
}
