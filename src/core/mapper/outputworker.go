package mapper

import (
	"github.com/KleinSamuel/gtamap/src/datawriter"
	"github.com/sirupsen/logrus"
	"sync"
)

func OutputWorker(taskQueue <-chan string, wg *sync.WaitGroup, writer *datawriter.Writer) {

	logrus.Debug("Started writeOutputWorker")

	for task := range taskQueue {
		writer.Write(task)
	}

	wg.Done()

	logrus.Debug("Finished writeOutputWorker")
}
