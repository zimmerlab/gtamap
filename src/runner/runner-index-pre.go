package runner

import (
	"strings"

	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/akamensky/argparse"
	"github.com/sirupsen/logrus"
)

type ArgsIndexPre struct {
	FastaFile            *string
	GtfFile              *string
	IsSeparateExtraction *bool
	OutputDir            *string
	GeneIds              *string
	UpstreamBases        *int
	DownstreamBases      *int
}

func AddCommandIndexPre(
	parser *argparse.Parser,
) (
	*argparse.Command,
	*ArgsIndexPre,
) {

	var command *argparse.Command = parser.NewCommand(
		"index-pre",
		"Extract gene sequences from genome.",
	)

	argsObj := &ArgsIndexPre{}

	argsObj.FastaFile = command.String(
		"",
		"fasta",
		&argparse.Options{
			Required: true,
			Help: "Nucleotide sequences (FASTA) file (currently only " +
				"non-compressed).",
		},
	)

	argsObj.GtfFile = command.String(
		"",
		"gtf",
		&argparse.Options{
			Required: true,
			Help: "Genome annotation (GTF) file (currently only " +
				"non-compressed).",
		},
	)

	argsObj.IsSeparateExtraction = command.Flag(
		"",
		"splitgenes",
		&argparse.Options{
			Help: "Extract gene sequences into separate fasta files, if " +
				"more than one gene id is specified in --geneids",
		},
	)

	argsObj.OutputDir = command.String(
		"",
		"output",
		&argparse.Options{
			Required: true,
			Help:     "Output directory for extracted gene sequences.",
		},
	)

	// TODO: maybe use StringList instead of manual parsing
	argsObj.GeneIds = command.String(
		"",
		"geneids",
		&argparse.Options{
			Required: false,
			Help:     "Gene IDs to extract (comma-separated).",
		},
	)

	argsObj.UpstreamBases = command.Int(
		"",
		"upstream", &argparse.Options{
			Required: false,
			Help: "Number of bases to add upstream of the gene start " +
				"position.",
			Default: 0,
		},
	)

	argsObj.DownstreamBases = command.Int(
		"",
		"downstream",
		&argparse.Options{
			Required: false,
			Help: "Number of bases to add downstream of the gene end " +
				"position.",
			Default: 0,
		},
	)

	return command, argsObj
}

func ExecIndexPre(argsObj *ArgsIndexPre) {

	logrus.Info("Extracting gene sequences from genome")

	// parse gene ids from comma separated string to map
	geneIds := make(map[string]struct{})
	if *argsObj.GeneIds != "" {
		genes := strings.Split(*argsObj.GeneIds, ",")
		for _, gene := range genes {
			geneIds[gene] = struct{}{}
		}
	}

	index.ExtractGeneSequenceFromGtfAndFastaForIndex(
		*argsObj.GtfFile,
		*argsObj.FastaFile,
		*argsObj.OutputDir,
		geneIds,
		*argsObj.UpstreamBases,
		*argsObj.DownstreamBases,
		*argsObj.IsSeparateExtraction,
	)

	logrus.Info("Done extracting gene sequences from genome")
}
