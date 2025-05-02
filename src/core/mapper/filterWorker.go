package mapper

import (
	"github.com/KleinSamuel/gtamap/src/datawriter"
	"github.com/sirupsen/logrus"
	"sync"
)

func FilterWorker(filterChannel <-chan string, wg *sync.WaitGroup, writer *datawriter.Writer) {

	logrus.Debug("Started writeOutputWorker")

	for id := range filterChannel {
		writer.Write("\n" + id)
	}

	wg.Done()

	logrus.Debug("Finished filterOutput")
}
