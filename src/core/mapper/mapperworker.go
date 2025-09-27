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

	start := time.Now()
	for task := range taskChan {

		// logrus.WithFields(logrus.Fields{
		// 	"workerId": workerId,
		// 	"task":     task.ID,
		// }).Debug("Processing task")

		// send progress update for every chunk of reads processed
		progressStats.ReadsProcessed++
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
		)
	}
	end := time.Now()
	duration := time.Since(start)

	progressChan <- events.Event{
		Type:  events.EventTypeMapperWorkerTime,
		Data:  uint64(duration),
		Start: start,
		End:   end,
	}

	//logrus.WithFields(logrus.Fields{
	//	"workerId": workerId,
	//}).Debug("Finished MapperWorker")
}
