package runner

import (
	"os"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/sirupsen/logrus"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

func GetCommandKmerOccurrences() *cobra.Command {
	var genomeFastaFilePath string
	var outputTsvFilePath string

	kmerOccCmd := &cobra.Command{
		Use:   "kmer-occurrences",
		Short: "Compute kmer occurrences",
		Run: func(cmd *cobra.Command, args []string) {
			SetConfigValue("kmer_occurrences.genome_fasta_file_path", genomeFastaFilePath)
			SetConfigValue("kmer_occurrences.output_tsv_file_path", outputTsvFilePath)

			if err := viper.Unmarshal(config.Mapper); err != nil {
				logrus.Fatalf("Unable to decode config: %v", err)
			}

			ExecKmerOccurrences()
		},
	}

	flags := kmerOccCmd.Flags()

	flags.StringVarP(
		&genomeFastaFilePath,
		"fasta",
		"f",
		"",
		"Genome fasta file to count kmers in (required)",
	)
	kmerOccCmd.MarkFlagRequired("genome-fasta")

	flags.StringVarP(
		&outputTsvFilePath,
		"output",
		"o",
		"",
		"Output TSV file path for kmer counts (required)",
	)
	kmerOccCmd.MarkFlagRequired("output-tsv")

	return kmerOccCmd
}

func ExecKmerOccurrences() {
	genomeFastaFile, errGenome := os.Open(config.Mapper.KmerOccurrences.GenomeFastaFilePath)
	if errGenome != nil {
		logrus.Fatalf("Could not open genome fasta file: %v", errGenome)
	}
	defer genomeFastaFile.Close()

	outputFile, errOutput := os.Create(config.Mapper.KmerOccurrences.OutputTsvFilePath)
	if errOutput != nil {
		logrus.Fatalf("Could not create output TSV file: %v", errOutput)
	}

	index.ComputeGenomeKmerOccurrences(genomeFastaFile, outputFile)

	logrus.Info("Done")
}
