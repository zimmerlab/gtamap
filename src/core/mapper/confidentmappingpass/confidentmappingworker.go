package confidentmappingpass

import (
	"sync"

	"github.com/sirupsen/logrus"
)

func ConfidentMappingWorker(confidentChan *ConfidentPassChan, wgConfidentMapping *sync.WaitGroup) {
	defer wgConfidentMapping.Done()

	for {
		task, ok := confidentChan.Receive()
		if !ok {
			break
		}
		logrus.Debugf("Using ReadPairMapping %s as confident readPair to determine introns and strandedness.", task.ReadPair.ReadR1.Header)
	}
	logrus.Info("Done with Annotation")
}
