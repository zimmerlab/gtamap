package mapper

import (
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mappedreadpair"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/core/timer"
	"github.com/sirupsen/logrus"
	"sync"
)

func MapperWorker(workerId int, genomeIndex *index.GenomeIndex,
	wg *sync.WaitGroup,
	taskChan <-chan MappingTask,
	secondPassChan *mapperutils.FourthPassChannel,
	resultChan chan<- *mappedreadpair.ReadPairMatchResult,
	progressChan chan<- bool,
	timerChan chan<- *timer.Timer) {

	defer wg.Done()

	logrus.WithFields(logrus.Fields{
		"workerId": workerId,
	}).Debug("Started MapperWorker")

	for task := range taskChan {

		logrus.WithFields(logrus.Fields{
			"workerId": workerId,
			"task":     task.ID,
		}).Debug("Processing task")

		mappingReadPairs, isMappable := MapReadPair(task.ReadPair, genomeIndex, secondPassChan, timerChan)

		if isMappable {
			for _, maps := range mappingReadPairs {
				resultChan <- &maps
			}
		}

		progressChan <- true
	}

	logrus.WithFields(logrus.Fields{
		"workerId": workerId,
	}).Debug("Finished MapperWorker")
}
