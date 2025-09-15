package runner

import (
	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/sirupsen/logrus"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

func GetCommandIndexPre() *cobra.Command {

	var fastaFilePath string
	var fastaIndexFilePath string
	var gtfFilePath string
	var singleFile bool
	var outputDirPath string
	var fastaFileName string
	var geneIds []string
	var upstreamBases int
	var downstreamBases int

	indexPreCmd := &cobra.Command{
		Use:   "index-pre",
		Short: "Extract gene sequences from a genome",
		Run: func(cmd *cobra.Command, args []string) {

			SetConfigValue("index.fasta_file_path", fastaFilePath)
			SetConfigValue("index.fasta_index_file_path", fastaIndexFilePath)
			SetConfigValue("index.gtf_file_path", gtfFilePath)
			SetConfigValue("index.output.single_file", singleFile)
			SetConfigValue("general.output_dir", outputDirPath)
			SetConfigValue("index.output.fasta_file_name", fastaFileName)
			SetConfigValue("index.gene_ids", geneIds)
			SetConfigValue("index.upstream_bases", upstreamBases)
			SetConfigValue("index.downstream_bases", downstreamBases)

			if err := viper.Unmarshal(config.Mapper); err != nil {
				logrus.Fatalf("Unable to decode config: %v", err)
			}

			config.Mapper.SetIndexFastaIndex(config.Mapper.Index.FastaIndexFilePath)

			ExecIndexPre()
		},
	}

	flags := indexPreCmd.Flags()

	flags.StringVarP(
		&fastaFilePath,
		"fasta",
		"f",
		"",
		"Fasta file (required)",
	)
	indexPreCmd.MarkFlagRequired("fasta")

	flags.StringVarP(
		&fastaIndexFilePath,
		"fasta-index",
		"i",
		"",
		"Fasta index file (default: [--fasta].fai)",
	)

	flags.StringVarP(
		&gtfFilePath,
		"gtf",
		"g",
		"",
		"Genome annotation file (.gtf) (required)",
	)
	indexPreCmd.MarkFlagRequired("gtf")

	flags.StringVarP(
		&outputDirPath,
		"output",
		"o",
		"",
		"Output directory",
	)

	flags.BoolVarP(
		&singleFile,
		"single-file",
		"s",
		false,
		"Write all gene sequences to a single fasta file",
	)

	flags.StringVar(
		&fastaFileName,
		"fasta-file-name",
		"",
		"Output FASTA file name (within output directory) (only for "+
			"--single-file) (default: genes.fa)",
	)

	flags.StringSliceVarP(
		&geneIds,
		"gene-ids",
		"l",
		nil,
		"Gene IDs to extract (comma-separated). If not provided, all "+
			"genes are extracted.",
	)

	flags.IntVarP(
		&upstreamBases,
		"upstream",
		"u",
		0,
		"Number of bases to add upstream of the gene start position.",
	)

	flags.IntVarP(
		&downstreamBases,
		"downstream",
		"d",
		0,
		"Number of bases to add downstream of the gene end position.",
	)

	return indexPreCmd
}

func ExecIndexPre() {

	logrus.Info("Extracting gene sequences from genome")

	// convert gene id list to set
	geneIds := make(map[string]struct{})
	for _, gene := range config.Mapper.Index.GeneIds {
		geneIds[gene] = struct{}{}
	}

	index.ExtractGeneSequenceFromGtfAndFastaForIndex(
		config.Mapper.Index.GtfFilePath,
		config.Mapper.Index.FastaFilePath,
		config.Mapper.Index.FastaIndexFilePath,
		config.Mapper.General.OutputDir,
		geneIds,
		config.Mapper.Index.UpstreamBases,
		config.Mapper.Index.DownstreamBases,
		config.Mapper.Index.Output.SingleFile,
	)
}
