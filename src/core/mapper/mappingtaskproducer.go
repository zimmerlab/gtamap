package mapper

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"strings"
	"sync"
)

func MappingTaskProducer(reader *fastq.Reader, taskChan chan<- MappingTask,
	maxReads int, specificQname string) {

	defer close(taskChan)

	// TODO: remove after testing
	taskCounter := 0
	foundSpecificQname := false

	// add each read pair as a mapping task to the task queue
	for readPair := reader.NextRead(); readPair != nil; readPair = reader.NextRead() {

		if specificQname != "" {
			name := strings.Split(readPair.ReadR1.Header, " ")[0]
			if name != specificQname {
				continue
			}
			foundSpecificQname = true
		}

		mappingTask := MappingTask{
			ID:       taskCounter,
			ReadPair: readPair,
		}
		taskChan <- mappingTask
		taskCounter++

		// TODO: remove after testing
		if maxReads > 0 && taskCounter >= maxReads {
			break
		}
		if foundSpecificQname {
			break
		}
	}

	fmt.Println("done producing tasks")
}

func MappingTaskWaiter(wgFirstPass *sync.WaitGroup, allTasksProcessedChan chan struct{}) {
	wgFirstPass.Wait()
	close(allTasksProcessedChan)
	fmt.Println("Closed taskChan as all mapper workers finished")
}
