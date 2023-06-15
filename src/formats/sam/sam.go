package sam

type Entry struct {
	Id    int
	Qname string
	Flag  int
	Rname string
	Pos   int
	Mapq  int
	Cigar string
	Rnext string
	Pnext int
	Tlen  int
	Seq   string
	Qual  string
}
