package main

import (
	"bytes"
	"flag"
	"fmt"
	"github.com/klauspost/pgzip" // For parallel gzip decompression
	"io"
	"os"
	"runtime"
	"sync"
	"time"
)

// decompressGzipFile decompresses a gzipped file and returns the raw content
func decompressGzipFile(filename string) ([]byte, error) {
	// Open the file
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("error opening file: %v", err)
	}
	defer file.Close()

	// Get file info for size estimation (compressed size)
	stat, err := file.Stat()
	if err != nil {
		return nil, fmt.Errorf("error getting file info: %v", err)
	}

	// Use pgzip for parallel decompression
	gz, err := pgzip.NewReader(file)
	if err != nil {
		return nil, fmt.Errorf("error creating gzip reader: %v", err)
	}
	defer gz.Close()

	// For large gzipped files, the decompressed size can be much larger
	// We'll estimate it at 5x the compressed size initially
	estimatedSize := stat.Size() * 5
	buffer := bytes.NewBuffer(make([]byte, 0, estimatedSize))

	// Copy the decompressed content to the buffer
	_, err = io.Copy(buffer, gz)
	if err != nil {
		return nil, fmt.Errorf("error decompressing file: %v", err)
	}

	return buffer.Bytes(), nil
}

// processMemoryMappedData processes a byte slice of uncompressed FASTQ data
func processMemoryMappedData(data []byte, numWorkers int) []FASTQEntry {
	// Find line breaks to split the data into chunks for parallel processing
	lineBreaks := make([]int, 0, len(data)/50) // Estimate based on avg line length
	for i := 0; i < len(data); i++ {
		if data[i] == '\n' {
			lineBreaks = append(lineBreaks, i)
		}
	}

	// Make sure we have a complete number of 4-line entries
	numLines := len(lineBreaks)
	completeEntriesLines := (numLines / 4) * 4

	// Calculate chunk boundaries ensuring each chunk contains complete FASTQ entries
	chunkSize := completeEntriesLines / numWorkers
	if chunkSize < 4 {
		chunkSize = 4 // At least one entry per worker
	}
	chunkSize = (chunkSize / 4) * 4 // Make sure it's a multiple of 4

	// Prepare chunks for workers
	var chunks [][]int // Start and end positions for each chunk
	for i := 0; i < completeEntriesLines; i += chunkSize {
		end := i + chunkSize
		if end > completeEntriesLines {
			end = completeEntriesLines
		}

		// Get actual byte positions
		startPos := 0
		if i > 0 {
			startPos = lineBreaks[i-1] + 1 // Start after previous line break
		}

		endPos := lineBreaks[end-1]
		chunks = append(chunks, []int{startPos, endPos})
	}

	// Process chunks in parallel
	var wg sync.WaitGroup
	results := make([][]FASTQEntry, len(chunks))

	for i, chunk := range chunks {
		wg.Add(1)
		go func(idx int, start, end int) {
			defer wg.Done()
			results[idx] = processChunk(data[start : end+1])
		}(i, chunk[0], chunk[1])
	}

	wg.Wait()

	// Combine results
	var allEntries []FASTQEntry
	for _, entries := range results {
		allEntries = append(allEntries, entries...)
	}

	return allEntries
}

func main() {
	// Parse command line arguments
	forwardFile := flag.String("forward", "", "Forward reads (gzipped FASTQ)")
	reverseFile := flag.String("reverse", "", "Reverse reads (gzipped FASTQ)")
	workers := flag.Int("workers", runtime.NumCPU(), "Number of worker goroutines")
	countOnly := flag.Bool("count", false, "Only count entries without storing them")
	validatePairs := flag.Bool("validate", true, "Validate that read pairs match")
	flag.Parse()

	if *forwardFile == "" || *reverseFile == "" {
		fmt.Println("Please provide both forward and reverse FASTQ files")
		flag.Usage()
		os.Exit(1)
	}

	startTime := time.Now()

	// Process files concurrently
	var fwdEntries, revEntries []FASTQEntry
	var fwdErr, revErr error
	var wg sync.WaitGroup

	wg.Add(2)
	go func() {
		defer wg.Done()
		// Decompress the gzipped file first
		fmt.Println("Decompressing forward reads...")
		fwdData, err := decompressGzipFile(*forwardFile)
		if err != nil {
			fwdErr = err
			return
		}

		// Process the decompressed data
		fmt.Println("Processing forward reads...")
		fwdEntries = processMemoryMappedData(fwdData, *workers)
	}()

	go func() {
		defer wg.Done()
		// Decompress the gzipped file first
		fmt.Println("Decompressing reverse reads...")
		revData, err := decompressGzipFile(*reverseFile)
		if err != nil {
			revErr = err
			return
		}

		// Process the decompressed data
		fmt.Println("Processing reverse reads...")
		revEntries = processMemoryMappedData(revData, *workers)
	}()

	wg.Wait()

	if fwdErr != nil {
		fmt.Fprintf(os.Stderr, "Error processing forward reads: %v\n", fwdErr)
		os.Exit(1)
	}

	if revErr != nil {
		fmt.Fprintf(os.Stderr, "Error processing reverse reads: %v\n", revErr)
		os.Exit(1)
	}

	// Now pair the reads
	var totalPairs int
	var allPairs []ReadPair

	// If validation is disabled, we can optimize by assuming reads are in the same order
	if !*validatePairs {
		totalPairs = len(fwdEntries)
		if !*countOnly {
			allPairs = make([]ReadPair, totalPairs)
			for i := 0; i < totalPairs; i++ {
				allPairs[i] = ReadPair{
					Forward: fwdEntries[i],
					Reverse: revEntries[i],
				}
			}
		}
	} else {
		// Pair the reads with validation
		pairs, err := pairReads(fwdEntries, revEntries)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error pairing reads: %v\n", err)
			os.Exit(1)
		}

		totalPairs = len(pairs)
		if !*countOnly {
			allPairs = pairs
		}
	}

	elapsed := time.Since(startTime)
	fmt.Printf("Processed %d read pairs in %v\n", totalPairs, elapsed)
	fmt.Printf("Speed: %.2f pairs/second\n", float64(totalPairs)/elapsed.Seconds())

	// Example of how to access the first few pairs
	if !*countOnly && len(allPairs) > 0 {
		fmt.Println("\nFirst read pair example:")
		fmt.Println("Forward:")
		fmt.Println("  Header:", allPairs[0].Forward.Header)
		fmt.Println("  Sequence:", allPairs[0].Forward.Sequence)
		fmt.Println("  Quality:", allPairs[0].Forward.Quality)
		fmt.Println("Reverse:")
		fmt.Println("  Header:", allPairs[0].Reverse.Header)
		fmt.Println("  Sequence:", allPairs[0].Reverse.Sequence)
		fmt.Println("  Quality:", allPairs[0].Reverse.Quality)
	}
}
