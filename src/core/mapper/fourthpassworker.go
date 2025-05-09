package mapper

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"sync"
)

// FourthPassWorker receives read paris that were not properly mapped in the first pass
func FourthPassWorker(fourthPassChan *mapperutils.FourthPassChannel, wgFouthsPass *sync.WaitGroup) {

	defer wgFouthsPass.Done()

	for {
		task, ok := fourthPassChan.Receive()
		if !ok {
			break
		}

		fmt.Println("fourth pass: ", task.ReadPair.ReadR1.Header)
	}
}
