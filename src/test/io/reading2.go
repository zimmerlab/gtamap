package main

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strings"
	"sync"
	"time"
)

// ReadPair represents a pair of forward and reverse reads
type ReadPair struct {
	Forward FASTQEntry
	Reverse FASTQEntry
}

// openGzippedFASTQ opens a gzipped FASTQ file and returns a buffered reader
func openGzippedFASTQ(filename string, bufferSize int) (*bufio.Reader, *os.File, *gzip.Reader, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, nil, nil, fmt.Errorf("error opening file: %v", err)
	}

	gzReader, err := gzip.NewReader(file)
	if err != nil {
		file.Close()
		return nil, nil, nil, fmt.Errorf("error creating gzip reader: %v", err)
	}

	reader := bufio.NewReaderSize(gzReader, bufferSize)
	return reader, file, gzReader, nil
}

// pairReads matches forward and reverse reads based on their headers
func pairReads(forwardEntries, reverseEntries []FASTQEntry) ([]ReadPair, error) {
	if len(forwardEntries) != len(reverseEntries) {
		return nil, fmt.Errorf("mismatched number of entries: forward=%d, reverse=%d",
			len(forwardEntries), len(reverseEntries))
	}

	pairs := make([]ReadPair, len(forwardEntries))

	// For best performance, we assume the reads are in the same order in both files
	for i := 0; i < len(forwardEntries); i++ {
		fwdID := getReadID(forwardEntries[i].Header)
		revID := getReadID(reverseEntries[i].Header)

		// Optional validation that the pairs match (can be disabled for speed)
		if fwdID != revID {
			return nil, fmt.Errorf("read pair mismatch at position %d: %s vs %s",
				i, fwdID, revID)
		}

		pairs[i] = ReadPair{
			Forward: forwardEntries[i],
			Reverse: reverseEntries[i],
		}
	}

	return pairs, nil
}

// getReadID extracts the read identifier from a FASTQ header
func getReadID(header string) string {
	// Common format: "@READID/1" or "@READID/2" for paired reads
	// Strip the trailing /1 or /2 to get the base ID
	parts := strings.Split(strings.TrimSpace(header), " ")
	baseID := parts[0]

	// Remove the trailing /1 or /2 if present
	if idx := strings.LastIndex(baseID, "/"); idx != -1 {
		return baseID[:idx]
	}

	return baseID
}

func main() {
	// Parse command line arguments
	//forwardFile := flag.String("forward", "", "Forward reads (gzipped FASTQ)")
	//reverseFile := flag.String("reverse", "", "Reverse reads (gzipped FASTQ)")
	//concurrency := flag.Int("concurrency", runtime.NumCPU(), "Number of worker goroutines")
	//chunkSize := flag.Int("chunk-size", 1024*1024, "Size of chunks to process")
	//countOnly := flag.Bool("count", false, "Only count entries without storing them")
	//validatePairs := flag.Bool("validate", true, "Validate that read pairs match")
	//flag.Parse()

	forwardFile := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_1.fastq.gz"
	reverseFile := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_2.fastq.gz"
	concurrency := 10
	chunkSize := 1024 * 1024
	//countOnly := true
	//validatePairs := true

	startTime := time.Now()

	// Open both FASTQ files
	fwdReader, fwdFile, fwdGz, err := openGzippedFASTQ(forwardFile, chunkSize*2)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error with forward reads: %v\n", err)
		os.Exit(1)
	}
	defer fwdFile.Close()
	defer fwdGz.Close()

	revReader, revFile, revGz, err := openGzippedFASTQ(reverseFile, chunkSize*2)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error with reverse reads: %v\n", err)
		os.Exit(1)
	}
	defer revFile.Close()
	defer revGz.Close()

	// Channels for distributing work
	fwdChunks := make(chan []byte, concurrency)
	revChunks := make(chan []byte, concurrency)

	// Channels for collecting results
	fwdResults := make(chan []FASTQEntry, concurrency)
	revResults := make(chan []FASTQEntry, concurrency)

	// Channel for final paired results
	//pairedResults := make(chan []ReadPair, concurrency)

	// Channel for signaling completion
	done := make(chan struct{})

	var wg sync.WaitGroup

	// Start worker goroutines for forward reads
	for i := 0; i < concurrency; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for chunk := range fwdChunks {
				fwdResults <- processChunk(chunk)
			}
		}()
	}

	// Start worker goroutines for reverse reads
	for i := 0; i < concurrency; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for chunk := range revChunks {
				revResults <- processChunk(chunk)
			}
		}()
	}

	// Start a goroutine to close the results channels when all workers are done
	go func() {
		wg.Wait()
		close(fwdResults)
		close(revResults)
	}()

	// Start a goroutine to read forward chunks
	go func() {
		for {
			chunk := make([]byte, chunkSize)
			n, err := fwdReader.Read(chunk)
			if err != nil && err != io.EOF {
				fmt.Fprintf(os.Stderr, "Error reading forward file: %v\n", err)
				close(fwdChunks)
				return
			}
			if n == 0 {
				close(fwdChunks)
				return
			}
			fwdChunks <- chunk[:n]
		}
	}()

	// Start a goroutine to read reverse chunks
	go func() {
		for {
			chunk := make([]byte, chunkSize)
			n, err := revReader.Read(chunk)
			if err != nil && err != io.EOF {
				fmt.Fprintf(os.Stderr, "Error reading reverse file: %v\n", err)
				close(revChunks)
				return
			}
			if n == 0 {
				close(revChunks)
				return
			}
			revChunks <- chunk[:n]
		}
	}()

	// Collect and buffer all forward entries
	var allForwardEntries []FASTQEntry
	go func() {
		for entries := range fwdResults {
			allForwardEntries = append(allForwardEntries, entries...)
		}
		done <- struct{}{}
	}()

	// Collect and buffer all reverse entries
	var allReverseEntries []FASTQEntry
	go func() {
		for entries := range revResults {
			allReverseEntries = append(allReverseEntries, entries...)
		}
		done <- struct{}{}
	}()

	// Wait for both collections to complete
	<-done
	<-done

	// Now pair the reads
	var totalPairs int
	//var allPairs []ReadPair

	// If validation is disabled, we can optimize by assuming reads are in the same order
	//if !validatePairs {
	//	totalPairs = len(allForwardEntries)
	//	if !countOnly {
	//		allPairs = make([]ReadPair, totalPairs)
	//		for i := 0; i < totalPairs; i++ {
	//			allPairs[i] = ReadPair{
	//				Forward: allForwardEntries[i],
	//				Reverse: allReverseEntries[i],
	//			}
	//		}
	//	}
	//} else {
	// Pair the reads with validation
	pairs, err := pairReads(allForwardEntries, allReverseEntries)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error pairing reads: %v\n", err)
		os.Exit(1)
	}

	totalPairs = len(pairs)
	//if !countOnly {
	//	allPairs = pairs
	//}
	//}

	elapsed := time.Since(startTime)
	fmt.Printf("Processed %d read pairs in %v\n", totalPairs, elapsed)
	fmt.Printf("Speed: %.2f pairs/second\n", float64(totalPairs)/elapsed.Seconds())

	//// Example of how to access the first few pairs
	//if !countOnly && len(allPairs) > 0 {
	//	fmt.Println("\nFirst read pair example:")
	//	fmt.Println("Forward:")
	//	fmt.Println("  Header:", allPairs[0].Forward.Header)
	//	fmt.Println("  Sequence:", allPairs[0].Forward.Sequence)
	//	fmt.Println("  Quality:", allPairs[0].Forward.Quality)
	//	fmt.Println("Reverse:")
	//	fmt.Println("  Header:", allPairs[0].Reverse.Header)
	//	fmt.Println("  Sequence:", allPairs[0].Reverse.Sequence)
	//	fmt.Println("  Quality:", allPairs[0].Reverse.Quality)
	//}
}
