package sam

import "fmt"

const Version string = "1.6"

type Header struct {
	Version                 string // version of the SAM format (e.g. 1.6)
	ReferenceSequenceName   string // the chromosome of the gene
	ReferenceSequenceLength int    // the length of the chromosome
	GenomeAnnotationVersion string // ensembl version of the gtf file (e.g. 108)
	GenomeAssemblyVersion   string // version of the genome (e.g. hg38)
	OrganismTaxId           string // taxonomy id of the organism (e.g. 9606)
	ToolVersion             string
	Transcripts             []*TranscriptInfo
}

type TranscriptInfo struct {
	Id                  string
	TranscriptEnsemblId string
	TranscriptLength    int
}

func (header *Header) String() string {
	headerString := fmt.Sprintf("@HD\tVN:%s\n", header.Version)
	headerString += fmt.Sprintf("@SQ\tSN:%s\tLN:%d\tSP:%s\tAS:%s\tDS:%s\n",
		header.ReferenceSequenceName, header.ReferenceSequenceLength, header.OrganismTaxId, header.GenomeAssemblyVersion, "annotation ensembl version "+header.GenomeAnnotationVersion)

	for _, transcriptInfo := range header.Transcripts {
		headerString += fmt.Sprintf("@SQ\tSN:%s\tLN:%d\tDS:%s\n", transcriptInfo.Id, transcriptInfo.TranscriptLength, transcriptInfo.TranscriptEnsemblId)
	}

	headerString += fmt.Sprintf("@PG\tID:%s\tPN:%s\tVN:%s\n", "GTAMap", "GTAMap", header.ToolVersion)

	return headerString
}

type Entry struct {
	Qname string // id of the read
	Flag  int    // bitwise flag
	Rname string // id of the reference
	Pos   int    // 1-based leftmost mapping position
	Mapq  int    // mapping quality
	Cigar string // CIGAR string (compact representation of alignment)
	Rnext string // id of the mate reference
	Pnext int    // position of the mate reference
	Tlen  int    // observed template length
	Seq   string // read sequence
	Qual  string // read quality
}

func (entry *Entry) String() string {
	return fmt.Sprintf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s", entry.Qname, entry.Flag, entry.Rname, entry.Pos, entry.Mapq, entry.Cigar, entry.Rnext, entry.Pnext, entry.Tlen, entry.Seq, entry.Qual)
}
