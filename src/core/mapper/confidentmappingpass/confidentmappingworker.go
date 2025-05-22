package confidentmappingpass

import (
	"sync"

	"github.com/sirupsen/logrus"
)

func ConfidentMappingWorker(confidentMappingChan *ConfidentMappingChannel, wgConfidentMapping *sync.WaitGroup) {
	defer wgConfidentMapping.Done()

	for {
		mappedReadPair, ok := confidentMappingChan.Receive()
		if !ok {
			break
		}
		logrus.Infof("Using ReadPairMapping %s as confident readPair to determine introns and strandedness.", mappedReadPair.ReadPair.ReadR1.Header)
	}
}
