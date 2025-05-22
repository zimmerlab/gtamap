package mapper

import (
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/mapper/confidentmappingpass"
	"github.com/KleinSamuel/gtamap/src/core/mapper/incompletemappingpass"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"

	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/timer"
	"github.com/sirupsen/logrus"
)

func MapperWorker(workerId int, genomeIndex *index.GenomeIndex,
	wg *sync.WaitGroup,
	taskChan <-chan MappingTask,
	incompleteMappingChan *incompletemappingpass.IncompleteMappingChannel,
	confidentMappingChan *confidentmappingpass.ConfidentMappingChannel,
	paralogMappingChan chan<- *mapperutils.ReadPairMatchResults,
	resultChan chan<- *mapperutils.ReadPairMatchResults,
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

		readPairMappings, isMappable := MapReadPair(task.ReadPair, genomeIndex, incompleteMappingChan, confidentMappingChan, paralogMappingChan, timerChan)

		if isMappable {
			resultChan <- readPairMappings
		}

		progressChan <- true
	}

	logrus.WithFields(logrus.Fields{
		"workerId": workerId,
	}).Debug("Finished MapperWorker")
}
