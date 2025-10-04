package mapper

import (
	"encoding/csv"
	"math"
	"os"
	"runtime"
	"strconv"
	"sync"
	"time"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/mapper/events"
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

func ProgressWorker(progressChan <-chan events.Event, wg *sync.WaitGroup) {
	defer wg.Done()

	logrus.Debug("Started progressWorker")

	timerStart := time.Now()

	tsvFile, err := os.Create(config.Mapper.GetMapperProgressLogPath())
	if err != nil {
		logrus.Errorf("Failed to create TSV file: %v", err)
		return
	}
	defer tsvFile.Close()

	tsvFileStageFile, err := os.Create(config.Mapper.GetMapperProgressStageTimeLogPath())
	if err != nil {
		logrus.Errorf("Failed to create TSV file: %v", err)
		return
	}
	defer tsvFileStageFile.Close()

	tsvWriter := csv.NewWriter(tsvFile)
	tsvWriter.Comma = '\t'
	defer tsvWriter.Flush()
	var numWriterRecords int = 0

	tsvStageFileWriter := csv.NewWriter(tsvFileStageFile)
	tsvStageFileWriter.Comma = '\t'
	defer tsvStageFileWriter.Flush()

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

	stage_timer_header := []string{
		"mappingTaskProducerTime",
		"mapperWorkerTime",
		"confidentWorkerTime",
		"secondPassWorkerTime",
		"outputWorkerTime",
		"mappingTaskProducerStart",
		"mappingTaskProducerEnd",
		"mapperWorkerStart",
		"mapperWorkerEnd",
		"confidentWorkerStart",
		"confidentWorkerEnd",
		"secondPassWorkerStart",
		"secondPassWorkerEnd",
		"outputWorkerStart",
		"outputWorkerEnd",
	}
	tsvStageFileWriter.Write(stage_timer_header)

	var mappingTaskProducerTime uint64 = 0
	var mapperWorkerTime uint64 = 0
	var confidentWorkerTime uint64 = 0
	var secondPassWorkerTime uint64 = 0
	var outputWorkerTime uint64 = 0

	var mappingTaskProducerStart time.Time
	var mappingTaskProducerEnd time.Time
	var mapperWorkerStart time.Time
	var mapperWorkerEnd time.Time
	var confidentWorkerStart time.Time
	var confidentWorkerEnd time.Time
	var secondPassWorkerStart time.Time
	// var secondPassWorkerEnd time.Time
	var outputWorkerStart time.Time
	var outputWorkerEnd time.Time

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

				stageRecord := []string{
					strconv.FormatUint(mappingTaskProducerTime, 10),
					strconv.FormatUint(mapperWorkerTime, 10),
					strconv.FormatUint(confidentWorkerTime, 10),
					strconv.FormatUint(secondPassWorkerTime, 10),
					strconv.FormatUint(outputWorkerTime, 10),
					mappingTaskProducerStart.Format("15:04:05"),
					mappingTaskProducerEnd.Format("15:04:05"),
					mapperWorkerStart.Format("15:04:05"),
					mapperWorkerEnd.Format("15:04:05"),
					confidentWorkerStart.Format("15:04:05"),
					confidentWorkerEnd.Format("15:04:05"),
					secondPassWorkerStart.Format("15:04:05"),
					outputWorkerStart.Format("15:04:05"), // second pass end == output start (had to solve like this due to go routines being hard to time when called in loop like in second pass. It is easier this way)
					outputWorkerStart.Format("15:04:05"),
					outputWorkerEnd.Format("15:04:05"),
				}

				tsvStageFileWriter.Write(stageRecord)

				return
			}

			switch event.Type {
			case events.EventFileSize:
				fileSize = event.Data
			case events.EventBytesProcessed:
				bytesProcessed = event.Data
			case events.EventTypeReadsProcessed:
				totalReadsProcessed += event.Data
			case events.EventTypeReadsAfterFiltering:
				totalReadsAfterFiltering += event.Data
			case events.EventTypeReadsMapped:
				totalReadsMapped += event.Data
			case events.EventTypeNumMappingLocations:
				numMappingLocations += event.Data
			case events.EventTypeNumConfidentMappings:
				numConfidentMappings += event.Data
			case events.EventTypeMapperProducerTime:
				mappingTaskProducerTime = event.Data
				mappingTaskProducerStart = event.Start
				mappingTaskProducerEnd = event.End
			case events.EventTypeMapperWorkerTime:
				mapperWorkerTime = event.Data
				mapperWorkerStart = event.Start
				mapperWorkerEnd = event.End
			case events.EventTypeConfidentWorkerTime:
				confidentWorkerTime = event.Data
				confidentWorkerStart = event.Start
				confidentWorkerEnd = event.End
			case events.EventTypeSecondPassWorkerTime:
				secondPassWorkerTime = event.Data
				secondPassWorkerStart = event.Start
				// secondPassWorkerEnd = event.End
			case events.EventTypeOutputWorkerTime:
				outputWorkerTime = event.Data
				outputWorkerStart = event.Start
				outputWorkerEnd = event.End
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
				"reads":    strconv.FormatFloat(math.Round(float64(totalReadsProcessed)/1_000_000*100)/100, 'f', 2, 64) + "M",
				"duration": utils.FormatDuration(elapsed),
				"rate":     strconv.FormatFloat(rate, 'f', 2, 64) + "M/min",
				"done":     percentFileRead + "%",
			}).Info("Progress update")
		}
	}
}
