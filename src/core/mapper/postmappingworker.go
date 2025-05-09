package mapper

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/sirupsen/logrus"
	"sync"
)

func PostMappingWorker(taskQueue <-chan *mapperutils.ReadPairMatchResult, wg *sync.WaitGroup) {

	logrus.Debug("Started accumulating mapped readpairs")

	results := make(map[int][]*mapperutils.ReadPairMatchResult)
	for task := range taskQueue {
		results[task.Fw.SequenceIndex/2] = append(results[task.Fw.SequenceIndex/2], task)
	}
	fmt.Println(len(results))

	wg.Done()
}
