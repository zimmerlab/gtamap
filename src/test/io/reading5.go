package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"os/exec"
	"sync"
	"syscall"
)

type FASTQEntryByte struct {
	Header   []byte
	Sequence []byte
	Quality  []byte
}

func main() {
	// Parse command-line arguments
	//filename := flag.String("file", "", "Gzipped file to process")
	//chunkSize := flag.Int("chunk-size", 4, "Number of lines per chunk")
	//outputFile := flag.String("output", "", "Output file (empty for stdout)")
	//bufferSizeMB := flag.Int("buffer", 8, "Buffer size in MB")
	//flag.Parse()
	filename := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_1.fastq.gz"
	chunkSize := 4
	bufferSizeMB := 1024

	// Execute zcat command with larger buffer for the pipe
	cmd := exec.Command("zcat", filename)

	// Set a larger buffer size for zcat output (helps with performance)
	// For some systems, this can significantly increase throughput
	if cmd.SysProcAttr == nil {
		cmd.SysProcAttr = &syscall.SysProcAttr{}
	}

	stdout, err := cmd.StdoutPipe()
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating stdout pipe: %v\n", err)
		os.Exit(1)
	}

	// Calculate buffer size in bytes
	bs := bufferSizeMB * 1024 * 1024

	// Create a buffered channel to avoid blocking
	taskQueue := make(chan FASTQEntryByte, 1000)

	// Create a channel to signal completion
	done := make(chan struct{})

	var wg sync.WaitGroup
	wg.Add(1)

	// Start processing in a goroutine
	go func() {
		defer wg.Done()
		defer close(taskQueue) // Close the channel when done processing
		err := processLines(stdout, chunkSize, taskQueue, bs)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error processing lines: %v\n", err)
			os.Exit(1)
		}
	}()

	// Count entries in a separate goroutine
	var totalEntries int
	go func() {
		for range taskQueue {
			totalEntries++
			
			if totalEntries%1_000_000 == 0 {
				fmt.Printf("Processed %d mio entries\n", totalEntries/1_000_000)
			}
		}
		close(done) // Signal that counting is complete
	}()

	// Start the command
	if err := cmd.Start(); err != nil {
		fmt.Fprintf(os.Stderr, "Error starting zcat: %v\n", err)
		os.Exit(1)
	}

	// Wait for processing to complete
	wg.Wait()

	// Wait for the command to finish
	if err := cmd.Wait(); err != nil {
		fmt.Fprintf(os.Stderr, "Error waiting for zcat: %v\n", err)
		os.Exit(1)
	}

	// Wait for counting to complete
	<-done

	fmt.Printf("Processed %d FASTQ entries\n", totalEntries)
}

func processLines(r io.Reader, chunkSize int, taskQueue chan<- FASTQEntryByte, bs int) error {
	// Use the provided buffer size
	scanner := bufio.NewScanner(r)
	buf := make([]byte, bs)
	scanner.Buffer(buf, bs)

	// Pre-allocate slice to reduce allocations
	lines := make([][]byte, 0, chunkSize)

	// Process the file line by line
	for scanner.Scan() {
		// Make a copy of the bytes to avoid issues with scanner reusing the buffer
		lineBytes := make([]byte, len(scanner.Bytes()))
		copy(lineBytes, scanner.Bytes())

		lines = append(lines, lineBytes)

		// When we have enough lines, create a FASTQ entry
		if len(lines) >= chunkSize {
			// Create fastq entry
			fastqEntry := FASTQEntryByte{
				Header:   lines[0],
				Sequence: lines[1],
				Quality:  lines[3],
			}

			// Send the entry to the task queue
			taskQueue <- fastqEntry

			// Clear the lines buffer while preserving capacity
			lines = lines[:0]
		}
	}

	// Don't forget the last chunk if it's not full but has complete entry
	if len(lines) == chunkSize {
		fastqEntry := FASTQEntryByte{
			Header:   lines[0],
			Sequence: lines[1],
			Quality:  lines[3],
		}

		// Send the entry to the task queue
		taskQueue <- fastqEntry
	}

	if err := scanner.Err(); err != nil {
		return fmt.Errorf("error scanning file: %w", err)
	}

	return nil
}
