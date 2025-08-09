package dataloader

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/KleinSamuel/gtamap/src/core/datastructure"
	"github.com/KleinSamuel/gtamap/src/formats/fasta"
	"github.com/KleinSamuel/gtamap/src/formats/gtf"
	"github.com/sirupsen/logrus"
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

func ExitIfFastaIndexIsMissing(pathFastaIndex string) {
	// check if the fasta index file exists
	_, errIndex := os.Stat(pathFastaIndex)
	if os.IsNotExist(errIndex) {
		logrus.Fatal("fasta index file (.fai) does not exist: " + pathFastaIndex)
	}
}

func GenerateInputForIndexFromFile(gtfFilePath string, fastaFilePath string, fastaIndexFilePath string) *gtf.Annotation {
	panicIfFilesAreMissing(gtfFilePath, fastaFilePath, fastaIndexFilePath)

	gtfFile, errGtf := os.Open(gtfFilePath)
	if errGtf != nil {
		logrus.Fatal("Error reading gtf file", errGtf)
	}

	fastaFile, errFasta := os.Open(fastaFilePath)
	if errFasta != nil {
		logrus.Fatal("Error reading fasta file", errFasta)
	}

	fastaIndexFile, errIndex := os.Open(fastaIndexFilePath)
	if errIndex != nil {
		logrus.Fatal("Error reading fasta index file (.fai)", errIndex)
	}

	return GenerateInputForIndex(gtfFile, fastaFile, fastaIndexFile)
}

// GenerateInputForIndex reads in the genome annotation from given gtf file,
// reads the fasta index (required) and extracts the dna sequence for every transcript
func GenerateInputForIndex(gtfFile *os.File, fastaFile *os.File, fastaIndexFile *os.File) *gtf.Annotation {
	logrus.WithFields(logrus.Fields{
		"gtfFile":        gtfFile.Name(),
		"fastaFile":      fastaFile.Name(),
		"fastaIndexFile": fastaIndexFile.Name(),
	}).Info("Generate input for index")

	fastaIndex := fasta.ReadFastaIndex(fastaIndexFile)
	annotation := gtf.ReadGtf(gtfFile)

	for _, gene := range annotation.Genes {

		logrus.WithFields(logrus.Fields{
			"geneIdEnsembl": gene.GeneIdEnsembl,
		}).Info("Process gene")

		// index := fastaIndex.Entries[gene.Chromosome]

		for _, trans := range gene.Transcripts {

			logrus.WithFields(logrus.Fields{
				"transcriptIdEnsembl": trans.TranscriptIdEnsembl,
			}).Info("Process transcript")

			// sort exons by start position such that transcript sequence is in order
			sort.Slice(trans.Exons, func(i, j int) bool {
				return trans.Exons[i].StartRelative < trans.Exons[j].StartRelative
			})

			var transcriptSequence []byte

			// build transcript sequence by seeking and concatenating the exon sequences
			for _, exon := range trans.Exons {

				logrus.WithFields(logrus.Fields{
					"startRelative": exon.StartRelative,
					"endRelative":   exon.EndRelative,
				}).Info("Extract exon sequence")

				startGenomic := gene.StartGenomic + exon.StartRelative
				endGenomic := gene.StartGenomic + exon.EndRelative

				transcriptSequence = append(transcriptSequence, ExtractSequenceAsBytesFromFasta(fastaFile, fastaIndex,
					gene.Chromosome, startGenomic, endGenomic)...)
			}

			trans.SequenceDna = strings.ReplaceAll(string(transcriptSequence), "\n", "")
		}

		// TODO: currently only supporting the first gene within the gtf file
		break
	}

	return annotation
}

func OpenBlacklistScanner(blackList string) (*bufio.Scanner, error) {
	ext := strings.ToLower(filepath.Ext(blackList))
	file, err := os.Open(blackList)
	if err != nil {
		return nil, err
	}

	if ext == ".gz" {
		gzReader, err := gzip.NewReader(file)
		if err != nil {
			file.Close()
			return nil, err
		}
		// Scanner over gzip stream
		return bufio.NewScanner(gzReader), nil
	}

	return bufio.NewScanner(file), nil
}

type FastaEntry struct {
	Header   string
	Sequence []byte
}

func ReadFasta(fastaFile *os.File) ([]*FastaEntry, error) {
	scanner := bufio.NewScanner(fastaFile)

	entries := make([]*FastaEntry, 0)

	// default buffer of 64KB is too small for long fasta seqs
	const maxCapacity = 10 * 1024 * 1024 // 10 MB
	buf := make([]byte, maxCapacity)
	scanner.Buffer(buf, maxCapacity)

	var header string
	var sequence []byte

	for scanner.Scan() {
		line := scanner.Bytes()

		if len(line) == 0 {
			continue
		}

		if line[0] == '>' {
			if header != "" {
				entries = append(entries, &FastaEntry{Header: header, Sequence: sequence})
			}
			// store new header and remove leading '>' character
			header = string(line)[1:]
			sequence = nil
		} else {
			sequence = append(sequence, bytes.TrimSpace(line)...)
		}
	}
	if header != "" {
		entries = append(entries, &FastaEntry{Header: header, Sequence: sequence})
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return entries, nil
}

func ExtractSequenceAsStringFromFastaUsingPath(fastaPath string, fastaIndex *fasta.Index, chromosome string,
	startGenomic uint32, endGenomic uint32,
) (string, error) {
	fastaFile, err := os.Open(fastaPath)
	if err != nil {
		return "", err
	}
	defer fastaFile.Close()

	return ExtractSequenceAsStringFromFasta(fastaFile, fastaIndex, chromosome, startGenomic, endGenomic), nil
}

func ExtractSequenceAsStringFromFasta(fastaFile *os.File, fastaIndex *fasta.Index, chromosome string,
	startGenomic uint32, endGenomic uint32,
) string {
	buffer := ExtractSequenceAsBytesFromFasta(fastaFile, fastaIndex, chromosome, startGenomic, endGenomic)
	return strings.ReplaceAll(string(buffer), "\n", "")
}

func ExtractBlacklistRegions(start, end int, contig string, blackList string) ([]*datastructure.Bounds, error) {
	bounds := make([]*datastructure.Bounds, 0)

	blackListFile, err := os.Open(blackList)
	if err != nil {
		return nil, err
	}
	defer blackListFile.Close()
	scanner := bufio.NewScanner(blackListFile)

	buf := make([]byte, 1000)
	scanner.Buffer(buf, 1000)

	for scanner.Scan() {
		line := scanner.Bytes()

		if len(line) == 0 {
			continue
		}
		fmt.Println(line)

	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return bounds, nil
}

func ExtractSequenceAsBytesFromFastaUsingPath(fastaPath string, fastaIndex *fasta.Index, chromosome string,
	startGenomic uint32, endGenomic uint32,
) ([]byte, error) {
	fastaFile, err := os.Open(fastaPath)
	if err != nil {
		return nil, err
	}
	defer fastaFile.Close()

	return ExtractSequenceAsBytesFromFasta(fastaFile, fastaIndex, chromosome, startGenomic, endGenomic), nil
}

func ExtractSequenceAsBytesFromFasta(fastaFile *os.File, fastaIndex *fasta.Index, chromosome string,
	startGenomic uint32, endGenomic uint32,
) []byte {
	index := fastaIndex.Entries[chromosome]

	// the start position of the transcript relative to its chromosome
	startRelativeChromosome := startGenomic
	// the number of newlines contained between the start of the chromosome and the start of the transcript
	lineNumber := int64(startRelativeChromosome) / int64(index.Linebases)
	// the number of bases contained in the last line that the exon is on
	offsetWithinLine := int64(startRelativeChromosome) % int64(index.Linebases)

	// the offset to the start of the transcript
	filePosition := index.Offset + lineNumber*int64(index.Linewidth) + offsetWithinLine

	_, errSeek := fastaFile.Seek(filePosition, 0)
	if errSeek != nil {
		panic(errSeek)
	}

	// the length of the exon sequence in nucleotides
	lengthBases := int64(endGenomic - startGenomic)
	// the number of characters in the current line until the next newline
	numCharsInLineLeft := int64(index.Linebases) - offsetWithinLine
	// the number of newlines contained between the start and end of the exon sequence
	numNewlines := int64(0)
	// adds a newline if the exon sequence spans across the current line
	if numCharsInLineLeft < lengthBases {
		numNewlines += 1
	}
	// adds the number of lines that the exon sequence spans across as newlines
	numNewlines += (lengthBases - numCharsInLineLeft) / int64(index.Linebases)

	// the total length of the exon sequence in bytes (including newlines)
	lengthBuffer := lengthBases + numNewlines

	buffer := make([]byte, lengthBuffer)

	_, errRead := fastaFile.Read(buffer)
	if errRead != nil {
		panic(errRead)
	}

	return buffer
}
