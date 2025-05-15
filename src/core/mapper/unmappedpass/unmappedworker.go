package unmappedpass

import (
	"fmt"
	"sync"
)

func UnmappedWorker(unmappedChan *UnmappedChannel, wgUnmapped *sync.WaitGroup) {
	defer wgUnmapped.Done()

	for {
		task, ok := unmappedChan.Receive()
		if !ok {
			break
		}

		fmt.Println("unmapped pass: ", task.ReadPair.ReadR1.Header)
	}
}
