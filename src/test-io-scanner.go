package main

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"runtime"
	"sync"
	"time"
)

const (
	numProcessors2 = 4
)

func main() {

	startTime := time.Now()

	file1Name := "/home/sam/Data/sra/SRR29933931.fastq"
	file2Name := "/home/sam/Data/sra/SRR29933931.copy.fastq"

	runtime.GOMAXPROCS(runtime.NumCPU())

	reader, _ := fastq.InitFromPaths(&file1Name, &file2Name)

	recordChan := make(chan fastq.ReadPair, 100)
	var wg sync.WaitGroup

	// Start worker goroutines
	for i := 0; i < numProcessors2; i++ {
		wg.Add(1)
		go processRecords2(recordChan, &wg)
	}

	// add each read pair as a mapping task to the task queue
	for readPair, _ := reader.NextRead(); readPair != nil; readPair, _ = reader.NextRead() {
		recordChan <- *readPair
	}

	fmt.Printf("Time taken: %v\n", time.Since(startTime))
}

func processRecords2(records <-chan fastq.ReadPair, wg *sync.WaitGroup) {
	defer wg.Done()
	for record := range records {
		if record.ReadR1 == nil || record.ReadR2 == nil {
			fmt.Printf("Invalid record: %v\n", record)
		}
	}
}
