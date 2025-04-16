package mapper

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"sync"
)

func SecondPassWorker(secondPassChan *mapperutils.SecondPassChannel, wgSecondPass *sync.WaitGroup) {

	defer wgSecondPass.Done()

	for {
		task, ok := secondPassChan.Receive()
		if !ok {
			break
		}

		fmt.Println("second pass: ", task.ReadPair.ReadR1.Header)
	}
}
