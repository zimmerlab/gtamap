package confidentmappingpass

import (
	"sync"

	"github.com/sirupsen/logrus"
)

func ConfidentMappingWorker(confidentMappingChan <-chan *ConfidentMappingTask, wgConfidentMapping *sync.WaitGroup) {
	defer wgConfidentMapping.Done()

	for mappedReadPair := range confidentMappingChan {
		logrus.Infof("Using ReadPairMapping %s as confident readPair to determine introns and strandedness.", mappedReadPair.ReadPair.ReadR1.Header)
	}
}
