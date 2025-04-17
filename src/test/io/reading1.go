package main

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"sync"
	"time"
)

// FASTQEntry represents a single FASTQ entry with 4 lines
type FASTQEntry struct {
	Header   string
	Sequence string
	Quality  string
}

// processChunk processes FASTQ entries from a byte chunk
func processChunk(chunk []byte) []FASTQEntry {
	entries := []FASTQEntry{}
	lines := []string{}
	start := 0

	// Split chunk into lines
	for i := 0; i < len(chunk); i++ {
		if chunk[i] == '\n' {
			lines = append(lines, string(chunk[start:i]))
			start = i + 1
		}
	}
	if start < len(chunk) {
		lines = append(lines, string(chunk[start:]))
	}

	// Group lines into FASTQ entries (4 lines per entry)
	for i := 0; i+3 < len(lines); i += 4 {
		entry := FASTQEntry{
			Header:   lines[i],
			Sequence: lines[i+1],
			Quality:  lines[i+3],
		}
		entries = append(entries, entry)
	}

	return entries
}

func main() {
	//filename := flag.String("file", "", "Gzipped FASTQ file to process")
	//concurrency := flag.Int("concurrency", runtime.NumCPU(), "Number of worker goroutines")
	//chunkSize := flag.Int("chunk-size", 1024*1024, "Size of chunks to process")
	//countOnly := flag.Bool("count", false, "Only count entries without storing them")
	//flag.Parse()

	filename := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_1.fastq.gz"
	concurrency := 10
	chunkSize := 1024 * 1024
	countOnly := true

	startTime := time.Now()

	// Open the gzipped file
	file, err := os.Open(filename)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening file: %v\n", err)
		os.Exit(1)
	}
	defer file.Close()

	// Create a gzip reader
	gzReader, err := gzip.NewReader(file)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating gzip reader: %v\n", err)
		os.Exit(1)
	}
	defer gzReader.Close()

	// Use a buffered reader for better performance
	reader := bufio.NewReaderSize(gzReader, chunkSize*2)

	// Channel for distributing work
	chunks := make(chan []byte, concurrency)
	// Channel for collecting results
	results := make(chan []FASTQEntry, concurrency)
	// Channel for signaling completion
	done := make(chan struct{})

	var wg sync.WaitGroup

	// Start worker goroutines
	for i := 0; i < concurrency; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for chunk := range chunks {
				results <- processChunk(chunk)
			}
		}()
	}

	// Start a goroutine to close the results channel when all workers are done
	go func() {
		wg.Wait()
		close(results)
		done <- struct{}{}
	}()

	// Start a goroutine to read chunks and send them to workers
	go func() {
		for {
			chunk := make([]byte, chunkSize)
			n, err := reader.Read(chunk)
			if err != nil && err != io.EOF {
				fmt.Fprintf(os.Stderr, "Error reading file: %v\n", err)
				close(chunks)
				return
			}
			if n == 0 {
				close(chunks)
				return
			}
			chunks <- chunk[:n]
		}
	}()

	// Collect and process results
	var totalEntries int
	var allEntries []FASTQEntry

	go func() {
		for entries := range results {
			totalEntries += len(entries)
			if !countOnly {
				allEntries = append(allEntries, entries...)
			}
		}
		done <- struct{}{}
	}()

	// Wait for both collection and processing to complete
	<-done
	<-done

	elapsed := time.Since(startTime)
	fmt.Printf("Processed %d FASTQ entries in %v\n", totalEntries, elapsed)
	fmt.Printf("Speed: %.2f entries/second\n", float64(totalEntries)/elapsed.Seconds())
}
