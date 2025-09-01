package mapper

import (
	"encoding/csv"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"sync"
	"time"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
)

type ProgressSnapshot struct {
	Timestamp            time.Time
	ReadsProcessed       uint64
	ReadsAfterFiltering  uint64
	ReadsMapped          uint64
	NumMappingLocations  uint64
	NumConfidentMappings uint64
	AllocKB              uint64
	HeapAllocKB          uint64
	StackAllocKB         uint64
	SysKB                uint64
	TotalBytes           uint64
	BytesProcessed       uint64
}

type Event struct {
	Type uint16
	Data uint64
}

// enum defining event types
const (
	EventTypeReadsProcessed       = 1
	EventTypeReadsAfterFiltering  = 2
	EventTypeReadsMapped          = 3
	EventTypeNumMappingLocations  = 4
	EventTypeNumConfidentMappings = 5
	EventFileSize                 = 6
	EventBytesProcessed           = 7
)

type MemorySample struct {
	Timestamp    time.Time
	AllocKB      uint64
	HeapAllocKB  uint64
	StackAllocKB uint64
	SysKB        uint64
	Goroutines   int
}

func GetMemorySample() MemorySample {
	var m runtime.MemStats
	runtime.ReadMemStats(&m)

	return MemorySample{
		Timestamp:    time.Now(),
		AllocKB:      m.Alloc / 1024,
		HeapAllocKB:  m.HeapAlloc / 1024,
		StackAllocKB: m.StackInuse / 1024,
		SysKB:        m.Sys / 1024,
		Goroutines:   runtime.NumGoroutine(),
	}
}

func ProgressWorker(progressChan <-chan Event, wg *sync.WaitGroup) {
	defer wg.Done()

	logrus.Debug("Started progressWorker")

	timerStart := time.Now()
	// Ensure parent dirs exist
	if err := os.MkdirAll(filepath.Dir(config.LogOut), 0o755); err != nil {
		logrus.Errorf("Failed to create directories for %s: %v", config.LogOut, err)
		return
	}

	tsvFile, err := os.Create(config.LogOut)
	if err != nil {
		logrus.Errorf("Failed to create TSV file: %v", err)
		return
	}
	defer tsvFile.Close()

	tsvWriter := csv.NewWriter(tsvFile)
	tsvWriter.Comma = '\t'
	defer tsvWriter.Flush()
	var numWriterRecords int = 0

	header := []string{
		"timestamp", "readsProcessed", "readsAfterFiltering", "readsMapped",
		"confidentMappings", "mappingLocations", "percentConfidentOfMapped",
		"percentFiltered", "percentMappedTotal", "percentMappedFilter",
		"meanLocationsPerRead", "allocKB", "heapAllocKB", "stackAllocKB",
		"totalFileSize", "bytesProcessed", "percentFileProcessed",
	}
	tsvWriter.Write(header)
	hasStarted := false

	var fileSize uint64 = 0
	var bytesProcessed uint64 = 0

	var totalReadsProcessed uint64 = 0
	var totalReadsAfterFiltering uint64 = 0
	var totalReadsMapped uint64 = 0
	var numMappingLocations uint64 = 0
	var numConfidentMappings uint64 = 0

	snapshots := make([]ProgressSnapshot, 0)

	ticker := time.NewTicker(3 * time.Second)
	defer ticker.Stop()

	for {
		select {
		case event, ok := <-progressChan:
			if !ok {
				logrus.Debug("Finished progressWorker")

				// write final progress snapshot
				memorySample := GetMemorySample()
				record := []string{
					time.Now().Format("2006-01-02 15:04:05"),
					strconv.FormatUint(totalReadsProcessed, 10),
					strconv.FormatUint(totalReadsAfterFiltering, 10),
					strconv.FormatUint(totalReadsMapped, 10),
					strconv.FormatUint(numConfidentMappings, 10),
					strconv.FormatUint(numMappingLocations, 10),
					strconv.FormatFloat(math.Round(float64(numConfidentMappings)/float64(totalReadsMapped)*10000)/100, 'f', 2, 64),
					strconv.FormatFloat(math.Round(float64(totalReadsAfterFiltering)/float64(totalReadsProcessed)*10000)/100, 'f', 2, 64),
					strconv.FormatFloat(math.Round(float64(totalReadsMapped)/float64(totalReadsProcessed)*10000)/100, 'f', 2, 64),
					strconv.FormatFloat(math.Round(float64(totalReadsMapped)/float64(totalReadsAfterFiltering)*10000)/100, 'f', 2, 64),
					strconv.FormatFloat(math.Round(float64(numMappingLocations)/float64(totalReadsMapped)*100)/100, 'f', 2, 64),
					strconv.FormatUint(memorySample.AllocKB, 10),
					strconv.FormatUint(memorySample.HeapAllocKB, 10),
					strconv.FormatUint(memorySample.StackAllocKB, 10),
					strconv.FormatUint(fileSize, 10),
					strconv.FormatUint(bytesProcessed, 10),
					strconv.FormatFloat(math.Round(float64(bytesProcessed)/float64(fileSize)*10000)/100, 'f', 2, 64),
				}
				tsvWriter.Write(record)

				return
			}

			switch event.Type {
			case EventFileSize:
				fileSize = event.Data
			case EventBytesProcessed:
				bytesProcessed = event.Data
			case EventTypeReadsProcessed:
				totalReadsProcessed += event.Data
			case EventTypeReadsAfterFiltering:
				totalReadsAfterFiltering += event.Data
			case EventTypeReadsMapped:
				totalReadsMapped += event.Data
			case EventTypeNumMappingLocations:
				numMappingLocations += event.Data
			case EventTypeNumConfidentMappings:
				numConfidentMappings += event.Data
			default:
				logrus.Warnf("Unknown event type: %d with data: %d", event.Type, event.Data)
			}

		case <-ticker.C:

			memorySample := GetMemorySample()

			snapshots = append(snapshots, ProgressSnapshot{
				Timestamp:            time.Now(),
				ReadsProcessed:       totalReadsProcessed,
				ReadsAfterFiltering:  totalReadsAfterFiltering,
				ReadsMapped:          totalReadsMapped,
				NumMappingLocations:  numMappingLocations,
				NumConfidentMappings: numConfidentMappings,
				AllocKB:              memorySample.AllocKB,
				HeapAllocKB:          memorySample.HeapAllocKB,
				StackAllocKB:         memorySample.StackAllocKB,
				SysKB:                memorySample.SysKB,
				TotalBytes:           fileSize,
				BytesProcessed:       bytesProcessed,
			})

			if !hasStarted {
				// initial datapoint
				record := []string{
					time.Now().Format("2006-01-02 15:04:05"),
					"0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0",
					strconv.FormatUint(totalReadsProcessed, 10),
					"0", "0",
				}
				tsvWriter.Write(record)
				hasStarted = true
			}

			percentFileRead := strconv.FormatFloat(math.Round(float64(bytesProcessed)/float64(fileSize)*10000)/100, 'f', 2, 64)

			record := []string{
				time.Now().Format("2006-01-02 15:04:05"),
				strconv.FormatUint(totalReadsProcessed, 10),
				strconv.FormatUint(totalReadsAfterFiltering, 10),
				strconv.FormatUint(totalReadsMapped, 10),
				strconv.FormatUint(numConfidentMappings, 10),
				strconv.FormatUint(numMappingLocations, 10),
				strconv.FormatFloat(math.Round(float64(numConfidentMappings)/float64(totalReadsMapped)*10000)/100, 'f', 2, 64),
				strconv.FormatFloat(math.Round(float64(totalReadsAfterFiltering)/float64(totalReadsProcessed)*10000)/100, 'f', 2, 64),
				strconv.FormatFloat(math.Round(float64(totalReadsMapped)/float64(totalReadsProcessed)*10000)/100, 'f', 2, 64),
				strconv.FormatFloat(math.Round(float64(totalReadsMapped)/float64(totalReadsAfterFiltering)*10000)/100, 'f', 2, 64),
				strconv.FormatFloat(math.Round(float64(numMappingLocations)/float64(totalReadsMapped)*100)/100, 'f', 2, 64),
				strconv.FormatUint(memorySample.AllocKB, 10),
				strconv.FormatUint(memorySample.HeapAllocKB, 10),
				strconv.FormatUint(memorySample.StackAllocKB, 10),
				strconv.FormatUint(fileSize, 10),
				strconv.FormatUint(bytesProcessed, 10),
				percentFileRead,
			}
			tsvWriter.Write(record)
			numWriterRecords += 1
			// if numWriterRecords%1 == 0 {
			tsvWriter.Flush()
			//}

			//logrus.WithFields(logrus.Fields{
			//	"readsProcessed":           totalReadsProcessed,
			//	"readsAfterFiltering":      totalReadsAfterFiltering,
			//	"readsMapped":              totalReadsMapped,
			//	"condfidentMappings":       numConfidentMappings,
			//	"mappingLocations":         numMappingLocations,
			//	"percentConfidentOfMapped": math.Round(float64(numConfidentMappings)/float64(totalReadsMapped)*10000) / 100,
			//	"percentFiltered":          math.Round(float64(totalReadsAfterFiltering)/float64(totalReadsProcessed)*10000) / 100,
			//	"percentMappedTotal":       math.Round(float64(totalReadsMapped)/float64(totalReadsProcessed)*10000) / 100,
			//	"percentMappedFilter":      math.Round(float64(totalReadsMapped)/float64(totalReadsAfterFiltering)*10000) / 100,
			//	"meanLocationsPerRead":     math.Round(float64(numMappingLocations)/float64(totalReadsMapped)*100) / 100,
			//	"allocKB":                  memorySample.AllocKB,
			//	"heapAllocKB":              memorySample.HeapAllocKB,
			//	"stackAllocKB":             memorySample.StackAllocKB,
			//}).Info("Progress update")

			elapsed := time.Since(timerStart)
			rate := float64(totalReadsProcessed) / elapsed.Minutes() / 1_000_000 // in million reads per minute

			logrus.WithFields(logrus.Fields{
				"reads": strconv.FormatFloat(math.Round(float64(totalReadsProcessed)/1_000_000*100)/100, 'f', 2, 64) + "M",
				"time":  utils.FormatDuration(elapsed),
				"rate":  strconv.FormatFloat(rate, 'f', 2, 64) + "M/min",
				"done":  percentFileRead + "%",
			}).Info("Progress update")
		}
	}
}
