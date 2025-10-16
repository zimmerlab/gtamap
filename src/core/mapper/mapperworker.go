package mapper

import (
	"sync"
	"time"

	"github.com/KleinSamuel/gtamap/src/core/mapper/confidentmappingpass"
	"github.com/KleinSamuel/gtamap/src/core/mapper/events"
	"github.com/KleinSamuel/gtamap/src/core/mapper/secondpass"

	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/timer"
)

type ProgressStats struct {
	ReadsProcessed       uint64
	ReadsAfterFiltering  uint64
	ReadsMapped          uint64
	NumMappingLocations  uint64
	NumConfidentMappings uint64
}

func MapperWorker(
	workerId int,
	genomeIndex *index.GenomeIndex,
	wg *sync.WaitGroup,
	taskChan <-chan MappingTask,
	secondpassChan *secondpass.SecondPassChannel,
	confidentMappingChan *confidentmappingpass.ConfidentPassChan,
	// paralogMappingChan chan<- *mapperutils.ReadPairMatchResults,
	progressChan chan<- events.Event,
	timerChan chan<- *timer.Timer,
	// mu *sync.Mutex,
) {
	defer wg.Done()

	//logrus.WithFields(logrus.Fields{
	//	"workerId": workerId,
	//}).Debug("Started MapperWorker")

	progressStats := &ProgressStats{
		ReadsProcessed:       0,
		ReadsAfterFiltering:  0,
		ReadsMapped:          0,
		NumMappingLocations:  0,
		NumConfidentMappings: 0,
	}

	var startMainMapping time.Time
	firstItem := true

	for task := range taskChan {
		if firstItem {
			firstItem = false
			startMainMapping = time.Now()
		}

		// logrus.WithFields(logrus.Fields{
		// 	"workerId": workerId,
		// 	"task":     task.ID,
		// }).Debug("Processing task")

		// send progress update for every chunk of reads processed
		// mu.Lock()
		progressStats.ReadsProcessed++
		// mu.Unlock()
		if progressStats.ReadsProcessed >= 1000 {
			progressChan <- events.Event{
				Type: events.EventTypeReadsProcessed,
				Data: uint64(progressStats.ReadsProcessed),
			}
			progressStats.ReadsProcessed = 0
		}

		MapReadPair(
			task.ReadPair,
			genomeIndex,
			secondpassChan,
			confidentMappingChan,
			timerChan,
			progressChan,
			progressStats,
			// mu,
		)
	}
	durationMainMapping := time.Since(startMainMapping)
	endMainMapping := time.Now()
	progressChan <- events.Event{
		Type:  events.EventTypeMapperWorkerTime,
		Data:  uint64(durationMainMapping),
		Start: startMainMapping,
		End:   endMainMapping,
	}
	if progressStats.ReadsProcessed > 0 {
		progressChan <- events.Event{
			Type: events.EventTypeReadsProcessed,
			Data: uint64(progressStats.ReadsProcessed),
		}
		progressStats.ReadsProcessed = 0
	}
	if progressStats.ReadsMapped > 0 {
		progressChan <- events.Event{
			Type: events.EventTypeReadsMapped,
			Data: progressStats.ReadsMapped,
		}
		progressStats.ReadsMapped = 0
	}
	if progressStats.ReadsAfterFiltering > 0 {
		progressChan <- events.Event{
			Type: events.EventTypeReadsAfterFiltering,
			Data: progressStats.ReadsAfterFiltering,
		}
		progressStats.ReadsAfterFiltering = 0
	}
	if progressStats.NumMappingLocations > 0 {
		progressChan <- events.Event{
			Type: events.EventTypeNumMappingLocations,
			Data: progressStats.NumMappingLocations,
		}
		progressStats.NumMappingLocations = 0
	}
	if progressStats.NumConfidentMappings > 0 {
		progressChan <- events.Event{
			Type: events.EventTypeNumConfidentMappings,
			Data: progressStats.NumConfidentMappings,
		}
		progressStats.NumConfidentMappings = 0
	}
	//logrus.WithFields(logrus.Fields{
	//	"workerId": workerId,
	//}).Debug("Finished MapperWorker")
}
