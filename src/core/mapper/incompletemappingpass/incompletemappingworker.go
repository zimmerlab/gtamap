package incompletemappingpass

import (
	"fmt"
	"sync"
)

func IncompleteMappingWorker(incompleteMappingChan *IncompleteMappingChannel, wgIncompleteMapping *sync.WaitGroup) {
	defer wgIncompleteMapping.Done()

	for {
		task, ok := incompleteMappingChan.Receive()
		if !ok {
			break
		}

		fmt.Println("incomplete mapping pass: ", task.ReadPair.ReadR1.Header)
	}
}
