package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"os/exec"
	"sync"
	"syscall"
	"time"
)

func main() {
	// Parse command-line arguments
	//filename := flag.String("file", "", "Gzipped file to process")
	//chunkSize := flag.Int("chunk-size", 4, "Number of lines per chunk")
	//outputFile := flag.String("output", "", "Output file (empty for stdout)")
	//bufferSizeMB := flag.Int("buffer", 8, "Buffer size in MB")
	//flag.Parse()
	filenameR1 := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_1.fastq.gz"
	filenameR2 := "/home/sam/Data/genomes/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3/NG-25876_HGT1_TAS2R4ko_lib434869_7080_3_2.fastq.gz"
	chunkSize := 4
	bufferSizeMB := 1024

	// Execute zcat command with larger buffer for the pipe
	cmdR1 := exec.Command("zcat", filenameR1)
	cmdR2 := exec.Command("zcat", filenameR2)

	// Set a larger buffer size for zcat output (helps with performance)
	// For some systems, this can significantly increase throughput
	if cmdR1.SysProcAttr == nil {
		cmdR1.SysProcAttr = &syscall.SysProcAttr{}
	}
	if cmdR2.SysProcAttr == nil {
		cmdR2.SysProcAttr = &syscall.SysProcAttr{}
	}

	stdoutR1, err := cmdR1.StdoutPipe()
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating stdout pipe: %v\n", err)
		os.Exit(1)
	}

	stdoutR2, err := cmdR2.StdoutPipe()
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating stdout pipe: %v\n", err)
		os.Exit(1)
	}

	// Calculate buffer size in bytes
	bs := bufferSizeMB * 1024 * 1024

	channelSize := 10_000

	// Create a buffered channel to avoid blocking
	taskQueueR1 := make(chan FASTQEntryByte, channelSize)
	taskQueueR2 := make(chan FASTQEntryByte, channelSize)

	out := make(chan PairedFASTQEntry, channelSize)

	// Create a channel to signal completion
	done := make(chan struct{})

	wgCombiner := sync.WaitGroup{}
	wgCombiner.Add(1)

	go mergeFASTQPairs(taskQueueR1, taskQueueR2, out, &wgCombiner)

	var wgR1 sync.WaitGroup
	wgR1.Add(1)

	var wgR2 sync.WaitGroup
	wgR2.Add(1)

	// Start processing in a goroutine
	go func() {
		defer wgR1.Done()
		defer close(taskQueueR1) // Close the channel when done processing
		err := processLines(stdoutR1, chunkSize, taskQueueR1, bs)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error processing lines: %v\n", err)
			os.Exit(1)
		}
	}()
	go func() {
		defer wgR2.Done()
		defer close(taskQueueR2) // Close the channel when done processing
		err := processLines(stdoutR2, chunkSize, taskQueueR2, bs)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error processing lines: %v\n", err)
			os.Exit(1)
		}
	}()

	//// Count entries in a separate goroutine
	//var totalEntries int
	//go func() {
	//	for range taskQueue {
	//		totalEntries++
	//
	//		if totalEntries%1_000_000 == 0 {
	//			fmt.Printf("Processed %d mio entries\n", totalEntries/1_000_000)
	//		}
	//	}
	//	close(done) // Signal that counting is complete
	//}()

	totalEntries := 0

	timeStart := time.Now()

	go func() {
		for range out {
			// Process the paired entry
			totalEntries++

			if totalEntries%1_000_000 == 0 {

				timeElapsed := time.Since(timeStart)

				rate := float64(totalEntries/1_000_000) / timeElapsed.Minutes()

				fmt.Printf("Processed %d mio pairs %.2f/min\n", totalEntries/1_000_000, rate)
			}
		}
		close(done) // Signal that counting is complete
	}()

	// Start the command
	if err := cmdR1.Start(); err != nil {
		fmt.Fprintf(os.Stderr, "Error starting zcat: %v\n", err)
		os.Exit(1)
	}
	if err := cmdR2.Start(); err != nil {
		fmt.Fprintf(os.Stderr, "Error starting zcat: %v\n", err)
		os.Exit(1)
	}

	// Wait for processing to complete
	wgR1.Wait()
	wgR2.Wait()

	// Wait for the command to finish
	if err := cmdR1.Wait(); err != nil {
		fmt.Fprintf(os.Stderr, "Error waiting for zcat: %v\n", err)
		os.Exit(1)
	}
	if err := cmdR2.Wait(); err != nil {
		fmt.Fprintf(os.Stderr, "Error waiting for zcat: %v\n", err)
		os.Exit(1)
	}

	// Wait for counting to complete
	<-done

	fmt.Printf("Processed %d FASTQ entries\n", totalEntries)
}

type FASTQEntryByte struct {
	Header   []byte
	Sequence []byte
	Quality  []byte
}

// PairedFASTQEntry represents a pair of FASTQ entries (e.g., forward and reverse reads)
type PairedFASTQEntry struct {
	Forward FASTQEntryByte
	Reverse FASTQEntryByte
}

// mergeFASTQPairs reads entries from two channels, pairs them, and sends them to an output channel
func mergeFASTQPairs(ch1, ch2 <-chan FASTQEntryByte, out chan<- PairedFASTQEntry, wg *sync.WaitGroup) {
	defer wg.Done()
	defer close(out)

	for {
		// Read from first channel
		entry1, ok1 := <-ch1
		if !ok1 {
			// Channel 1 is closed, we're done
			fmt.Println("Channel 1 closed, stopping merger")
			break
		}

		// Read from second channel
		entry2, ok2 := <-ch2
		if !ok2 {
			// Channel 2 is closed, we're done
			fmt.Println("Channel 2 closed, stopping merger")
			break
		}

		// Create paired entry
		pairedEntry := PairedFASTQEntry{
			Forward: entry1,
			Reverse: entry2,
		}

		// Send paired entry to output channel
		out <- pairedEntry
	}
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
