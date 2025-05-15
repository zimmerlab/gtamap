package mapper

import (
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/mapper/postmappingpass"
	"github.com/KleinSamuel/gtamap/src/core/mapper/unmappedpass"

	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/timer"
	"github.com/sirupsen/logrus"
)

func MapperWorker(workerId int, genomeIndex *index.GenomeIndex,
	wg *sync.WaitGroup,
	taskChan <-chan MappingTask,
	unmappedChan *unmappedpass.UnmappedChannel,
	resultChan chan<- *postmappingpass.ReadPairMatchResults,
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

		readPairMappings, isMappable := MapReadPair(task.ReadPair, genomeIndex, unmappedChan, timerChan)

		if isMappable {
			resultChan <- readPairMappings
		}

		progressChan <- true
	}

	logrus.WithFields(logrus.Fields{
		"workerId": workerId,
	}).Debug("Finished MapperWorker")
}
