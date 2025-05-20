package mapper

import (
	"fmt"
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
)

func unmappedWorker(unmappedChan *mapperutils.UnmappedChannel, wgUnmapped *sync.WaitGroup) {
	defer wgUnmapped.Done()

	for {
		task, ok := unmappedChan.Receive()
		if !ok {
			break
		}

		fmt.Println("unmapped pass: ", task.ReadPair.ReadR1.Header)
	}
}
