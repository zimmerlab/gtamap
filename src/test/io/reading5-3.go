package main

import (
	"bufio"
	"fmt"
	"github.com/sirupsen/logrus"
	"io"
	"os"
	"strings"
	"sync"
	"time"

	"github.com/klauspost/pgzip"
)

type FASTQEntryByte struct {
	Header   *[]byte
	Sequence *[]byte
	Quality  *[]byte
}

type PairedFASTQEntry struct {
	Forward *FASTQEntryByte
	Reverse *FASTQEntryByte
}

func main() {

	//filenameR1 := flag.String("r1", "", "Gzipped FASTQ file to process")
	//filenameR2 := flag.String("r2", "", "Gzipped FASTQ file to process")
	//flag.Parse()

	filenameR1 := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_1.fastq.gz"
	filenameR2 := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_2.fastq.gz"
	bufSize := 4 << 20 // 4 MiB read buffer
	chanDepth := 4     // Small buffered channels, we are IO bound

	// Open gzipped readers with parallel decompression
	r1, err := openGzip(filenameR1, bufSize)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Failed to open R1: %v\n", err)
		os.Exit(1)
	}
	r2, err := openGzip(filenameR2, bufSize)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Failed to open R2: %v\n", err)
		os.Exit(1)
	}

	taskQueueR1 := make(chan FASTQEntryByte, chanDepth)
	taskQueueR2 := make(chan FASTQEntryByte, chanDepth)
	out := make(chan PairedFASTQEntry, chanDepth)

	var wg sync.WaitGroup
	wg.Add(3)

	// Read R1
	go func() {
		defer wg.Done()
		err := readFASTQEntries(r1, taskQueueR1)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error reading R1: %v\n", err)
			os.Exit(1)
		}
		close(taskQueueR1)
	}()

	// Read R2
	go func() {
		defer wg.Done()
		err := readFASTQEntries(r2, taskQueueR2)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error reading R2: %v\n", err)
			os.Exit(1)
		}
		close(taskQueueR2)
	}()

	// Merge paired reads
	go func() {
		defer wg.Done()
		mergeFASTQPairs(taskQueueR1, taskQueueR2, out)
		close(out)
	}()

	// Consume pairs
	total := 0
	start := time.Now()
	for pair := range out {

		//fmt.Println("\t\t--r1--")
		//fmt.Println("\t\t", string(*pair.Forward.Header))
		//fmt.Println("\t\t", string(*pair.Forward.Sequence))
		//fmt.Println("\t\t", string(*pair.Forward.Quality))
		//fmt.Println("\t\t--r2--")
		//fmt.Println("\t\t", string(*pair.Reverse.Header))
		//fmt.Println("\t\t", string(*pair.Reverse.Sequence))
		//fmt.Println("\t\t", string(*pair.Reverse.Quality))

		h1 := strings.Split(string(*pair.Forward.Header), " ")[0]
		h2 := strings.Split(string(*pair.Reverse.Header), " ")[0]

		if h1 != h2 {
			logrus.Errorf("Headers do not match")
			logrus.Errorf("R1: %s", string(*pair.Forward.Header))
			logrus.Errorf("R2: %s", string(*pair.Reverse.Header))
			logrus.Fatal("Headers do not match")
		}

		total++
		if total%1_000_000 == 0 {
			elapsed := time.Since(start)
			rate := float64(total/1_000_000) / elapsed.Minutes()
			fmt.Printf("Processed %d mio pairs (%.2f mio/min)\n", total/1_000_000, rate)
		}
	}

	wg.Wait()
	fmt.Printf("Finished: %d paired FASTQ entries processed\n", total)
}

// Opens a parallel gzip reader with a large buffer
func openGzip(filename string, bufSize int) (*bufio.Reader, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	zr, err := pgzip.NewReader(f)
	if err != nil {
		f.Close()
		return nil, err
	}
	// pgzip.Reader will close `f` when zr.Close() is called
	return bufio.NewReaderSize(zr, bufSize), nil
}

// Reads 4 lines per entry, sends entries to the channel
func readFASTQEntries(r *bufio.Reader, out chan<- FASTQEntryByte) error {
	for {
		header, err := r.ReadSlice('\n')
		if err != nil {
			if err == io.EOF {
				break
			}
			return err
		}
		seq, err := r.ReadSlice('\n')
		if err != nil {
			return err
		}
		plus, err := r.ReadSlice('\n')
		if err != nil {
			return err
		}
		qual, err := r.ReadSlice('\n')
		if err != nil {
			return err
		}

		// Minimal copy: just trim trailing '\n'
		entry := FASTQEntryByte{
			Header:   trimNewline(header),
			Sequence: trimNewline(seq),
			Quality:  trimNewline(qual),
		}

		//fmt.Println("--")
		//fmt.Println(string(entry.Header))
		//fmt.Println(string(entry.Sequence))
		//fmt.Println(string(entry.Quality))

		_ = plus // unused, but needed to read 4 lines
		out <- entry
	}
	return nil
}

// Strips trailing newline
func trimNewline(b []byte) *[]byte {
	// copy the slice to avoid modifying the original
	c := make([]byte, len(b))

	copy(c, b)

	if len(c) > 0 && b[len(c)-1] == '\n' {
		c = c[:len(c)-1]
	}
	return &c
}

// Pairs up R1 and R2 entries
func mergeFASTQPairs(ch1, ch2 <-chan FASTQEntryByte, out chan<- PairedFASTQEntry) {
	for {
		entry1, ok1 := <-ch1
		entry2, ok2 := <-ch2
		if !ok1 || !ok2 {
			break
		}
		out <- PairedFASTQEntry{Forward: &entry1, Reverse: &entry2}
	}
}
