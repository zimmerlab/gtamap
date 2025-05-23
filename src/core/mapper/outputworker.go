package mapper

import (
	"sync"

	"github.com/KleinSamuel/gtamap/src/datawriter"
	"github.com/sirupsen/logrus"
)

func OutputWorker(taskQueue <-chan string, wg *sync.WaitGroup, writer *datawriter.Writer) {
	logrus.Debug("Started writeOutputWorker")
	defer wg.Done()

	for task := range taskQueue {
		writer.Write(task)
	}

	logrus.Debug("Finished writeOutputWorker")
}
