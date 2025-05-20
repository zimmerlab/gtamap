package mapper

import (
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/core/timer"
	"github.com/sirupsen/logrus"
)

func MapperWorker(workerId int, genomeIndex *index.GenomeIndex,
	wg *sync.WaitGroup,
	taskChan <-chan MappingTask,
	unmappedChan *mapperutils.UnmappedChannel,
	outputChan chan<- string,
	progressChan chan<- bool,
	timerChan chan<- *timer.Timer,
) {
	defer wg.Done()

	logrus.WithFields(logrus.Fields{
		"workerId": workerId,
	}).Debug("Started MapperWorker")

	for task := range taskChan {

		logrus.WithFields(logrus.Fields{
			"workerId": workerId,
			"task":     task.ID,
		}).Debug("Processing task")

		result, isMappable := MapReadPair(task.ReadPair, genomeIndex, unmappedChan, timerChan)

		if isMappable {
			outputChan <- result
		}

		progressChan <- true
	}

	logrus.WithFields(logrus.Fields{
		"workerId": workerId,
	}).Debug("Finished MapperWorker")
}
