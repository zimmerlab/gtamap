package incompletemappingpass

import (
	"sync"

	"github.com/sirupsen/logrus"
)

func IncompleteMappingWorker(incompleteMappingChan *IncompleteMappingChannel, wgIncompleteMapping *sync.WaitGroup) {
	defer wgIncompleteMapping.Done()

	for {
		task, ok := incompleteMappingChan.Receive()
		if !ok {
			break
		}

		logrus.Debugf("Incomplete mapping pass for readPair: %s", task.ReadPair.ReadR1.Header)
	}
}
