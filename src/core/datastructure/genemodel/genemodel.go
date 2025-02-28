package genemodel

type Gene struct {
	GeneIdEnsembl  string // e.g. "ENSG00000173585"
	Chromosome     string // e.g. "1", "2", "X", "Y", "MT"
	IsFowardStrand bool   // true if on forward strand
	StartGenomic   uint32 // 0-based genomic start location
	EndGenomic     uint32 // exclusive genomic end location
}

type Transcript struct {
	TranscriptIdEnsembl       string // e.g. "ENST00000342992"
	SequenceDnaForward53Index int    // the index of the forward sequence in the sequence list of the suffix tree
	SequenceDnaReverse53Index int    // the index in the reverse sequence in the sequence list of the suffix tree
	SequenceLength            int    // the length of the dna sequence
	Exons                     []*Exon
	SequenceDnaForward        string
	SequenceDnaReverse        string
}

type Exon struct {
	StartRelative uint32 // 0-based start location relative to genomic location of parent gene
	EndRelative   uint32 // exclusive end location relative to genomic location of parent gene
}

type EquivalenceClass struct {
	Id          uint32
	FromGenomic uint32
	ToGenomic   uint32
}

type ExonJunction struct {
	Id          uint32
	FromGenomic uint32
	ToGenomic   uint32
}
