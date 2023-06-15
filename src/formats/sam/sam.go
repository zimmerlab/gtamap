package sam

import "fmt"

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
