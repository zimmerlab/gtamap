package secondpass

import (
	"sync"

	"github.com/sirupsen/logrus"
)

func SecondpassMappingWorker(incompleteMappingChan *SecondPassChannel, wgIncompleteMapping *sync.WaitGroup) {
	// in here we receive all non confident readpairs. Some have multiple maps for fw and rv and some are
	// not completely mapped yet.
	defer wgIncompleteMapping.Done()

	for {
		task, ok := incompleteMappingChan.Receive()
		if !ok {
			break
		}

		logrus.Debugf("Incomplete mapping pass for readPair: %s", task.ReadPair.ReadR1.Header)
	}
}
