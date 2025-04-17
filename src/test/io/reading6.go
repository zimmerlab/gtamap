package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"runtime"
	"sync"
	"time"

	"github.com/klauspost/pgzip"
)

type FASTQEntryByte struct {
	Header, Sequence, Quality []byte
}

type PairedFASTQEntry struct {
	Forward FASTQEntryByte
	Reverse FASTQEntryByte
}

var entryPool = sync.Pool{
	New: func() interface{} {
		return &PairedFASTQEntry{}
	},
}

func main() {

	filenameR1 := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_1.fastq.gz"
	filenameR2 := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_2.fastq.gz"

	bufSize := 8 << 20 // 8 MiB read buffer
	numWorkers := runtime.NumCPU()

	// ensure pgzip uses all cores
	runtime.GOMAXPROCS(runtime.NumCPU() * 2)

	// open both gzipped FASTQ inputs
	r1 := mustOpenFASTQ(filenameR1, bufSize)
	r2 := mustOpenFASTQ(filenameR2, bufSize)

	// channel of work items
	workChan := make(chan *PairedFASTQEntry, numWorkers*8)

	// spawn workers
	var wg sync.WaitGroup
	var total uint64
	start := time.Now()

	wg.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go func() {
			defer wg.Done()
			for p := range workChan {
				// ─── your processing logic here ───
				// e.g. write p to file, analyze, etc.
				_ = p

				//// count & occasional progress report
				//cnt := atomic.AddUint64(&total, 1)
				//if cnt%1_000_000 == 0 {
				//	elapsed := time.Since(start)
				//	rate := float64(cnt/1_000_000) / elapsed.Minutes()
				//	fmt.Printf("Processed %d M pairs (%.2f M/min)\n",
				//		cnt/1_000_000, rate)
				//}

				// recycle the struct for reuse
				entryPool.Put(p)
			}
		}()
	}

	counter := 0

	// read‑and‑dispatch loop (runs in main goroutine)
	for {
		h1, s1, q1, err1 := readOne(r1)
		h2, s2, q2, err2 := readOne(r2)
		if err1 == io.EOF || err2 == io.EOF {
			break
		}
		if err1 != nil || err2 != nil {
			fmt.Fprintf(os.Stderr, "read error: %v / %v\n", err1, err2)
			os.Exit(1)
		}

		// get a pooled entry, copy in data
		p := entryPool.Get().(*PairedFASTQEntry)
		p.Forward.Header = append(p.Forward.Header[:0], h1...)
		p.Forward.Sequence = append(p.Forward.Sequence[:0], s1...)
		p.Forward.Quality = append(p.Forward.Quality[:0], q1...)
		p.Reverse.Header = append(p.Reverse.Header[:0], h2...)
		p.Reverse.Sequence = append(p.Reverse.Sequence[:0], s2...)
		p.Reverse.Quality = append(p.Reverse.Quality[:0], q2...)

		// hand off to a worker
		workChan <- p

		counter++
		if counter%1_000_000 == 0 {
			elapsed := time.Since(start)
			rate := float64(counter/1_000_000) / elapsed.Minutes()
			fmt.Printf("Read %d M pairs (%.2f M/min)\n",
				counter/1_000_000, rate)
		}
	}

	// no more work
	close(workChan)
	// wait for all workers to finish
	wg.Wait()

	fmt.Printf("Done: %d paired FASTQ entries in %v\n", total, time.Since(start))
}

// mustOpenFASTQ opens a .gz file with a parallel pgzip reader + bufio buffer.
func mustOpenFASTQ(path string, bufSize int) *bufio.Reader {
	f, err := os.Open(path)
	if err != nil {
		fmt.Fprintf(os.Stderr, "open %s: %v\n", path, err)
		os.Exit(1)
	}
	zr, err := pgzip.NewReaderN(f, bufSize, 32)
	if err != nil {
		fmt.Fprintf(os.Stderr, "gzip %s: %v\n", path, err)
		os.Exit(1)
	}
	// let pgzip spawn one worker per CPU, 1 MiB block size
	//zr.SetConcurrency(runtime.NumCPU(), 1<<20)
	// closing zr will also close f
	return bufio.NewReaderSize(zr, bufSize)
}

// readOne reads exactly 4 lines (header, seq, plus, qual) and returns hdr/seq/qual without '\n'.
func readOne(r *bufio.Reader) (hdr, seq, qual []byte, err error) {
	var line []byte
	if line, err = r.ReadSlice('\n'); err != nil {
		return
	}
	hdr = trimNL(line)
	if line, err = r.ReadSlice('\n'); err != nil {
		return
	}
	seq = trimNL(line)
	if _, err = r.ReadSlice('\n'); err != nil { // skip '+'
		return
	}
	if line, err = r.ReadSlice('\n'); err != nil {
		return
	}
	qual = trimNL(line)
	return
}

// trimNL removes a trailing newline byte if present.
func trimNL(b []byte) []byte {
	if len(b) > 0 && b[len(b)-1] == '\n' {
		return b[:len(b)-1]
	}
	return b
}
