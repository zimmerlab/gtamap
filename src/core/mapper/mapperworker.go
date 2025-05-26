package mapper

import (
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/mapper/confidentmappingpass"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/core/mapper/secondpass"

	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/timer"
	"github.com/sirupsen/logrus"
)

func MapperWorker(workerId int, genomeIndex *index.GenomeIndex,
	wg *sync.WaitGroup,
	taskChan <-chan MappingTask,
	secondpassChan *secondpass.SecondPassChannel,
	confidentMappingChan chan<- *confidentmappingpass.ConfidentMappingTask,
	paralogMappingChan chan<- *mapperutils.ReadPairMatchResults,
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

		MapReadPair(task.ReadPair, genomeIndex, secondpassChan, confidentMappingChan, timerChan, paralogMappingChan)

		progressChan <- true
	}

	logrus.WithFields(logrus.Fields{
		"workerId": workerId,
	}).Debug("Finished MapperWorker")
}
