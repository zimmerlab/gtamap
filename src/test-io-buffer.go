package main

import (
	"bufio"
	"fmt"
	"os"
	"runtime"
	"sync"
	"time"
)

const (
	file1Name     = "/home/sam/Data/sra/SRR29933931.fastq"
	file2Name     = "/home/sam/Data/sra/SRR29933931.copy.fastq"
	recordLines   = 4
	bufferSize    = 64 * 1024 * 1024 // MB buffer
	numProcessors = 4
)

type Record struct {
	file1Lines []string
	file2Lines []string
}

func main() {

	startTime := time.Now()

	runtime.GOMAXPROCS(runtime.NumCPU())

	recordChan := make(chan Record, 100)
	var wg sync.WaitGroup

	// Start worker goroutines
	for i := 0; i < numProcessors; i++ {
		wg.Add(1)
		go processRecords(recordChan, &wg)
	}

	// Read both files concurrently
	go readFiles(file1Name, file2Name, recordChan)

	// Wait for all records to be processed
	wg.Wait()

	fmt.Printf("Time taken: %v\n", time.Since(startTime))
}

func readFiles(file1Name, file2Name string, recordChan chan<- Record) {
	file1, err1 := os.Open(file1Name)
	file2, err2 := os.Open(file2Name)
	if err1 != nil || err2 != nil {
		fmt.Printf("Error opening files: %v, %v\n", err1, err2)
		close(recordChan)
		return
	}
	defer file1.Close()
	defer file2.Close()

	reader1 := bufio.NewReaderSize(file1, bufferSize)
	reader2 := bufio.NewReaderSize(file2, bufferSize)

	for {
		record := Record{
			file1Lines: make([]string, recordLines),
			file2Lines: make([]string, recordLines),
		}

		eof := false
		for i := 0; i < recordLines; i++ {
			line1, err1 := reader1.ReadString('\n')
			line2, err2 := reader2.ReadString('\n')

			if err1 != nil || err2 != nil {
				eof = true
				break
			}

			record.file1Lines[i] = line1
			record.file2Lines[i] = line2
		}

		if eof {
			break
		}

		recordChan <- record
	}

	close(recordChan)
}

func processRecords(records <-chan Record, wg *sync.WaitGroup) {
	defer wg.Done()
	for record := range records {
		if len(record.file1Lines) != recordLines || len(record.file2Lines) != recordLines {
			fmt.Printf("Invalid record: %v\n", record)
		}
	}
}
