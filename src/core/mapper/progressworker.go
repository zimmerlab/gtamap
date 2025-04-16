package mapper

import (
	"github.com/sirupsen/logrus"
	"sync"
)

func ProgressWorker(progressChan <-chan bool, wg *sync.WaitGroup) {
	logrus.Debug("Started progressWorker")

	numTasksDone := 0

	for _ = range progressChan {
		numTasksDone++
		if numTasksDone%1_000_000 == 0 {
			logrus.WithFields(logrus.Fields{
				"readPairsDone": numTasksDone,
			}).Info("Progress update")
		}
	}

	wg.Done()

	logrus.Debug("Finished progressWorker")
}
