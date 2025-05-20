package main

import (
	"fmt"
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
)

func main() {
	fmt.Println("hello")

	taskChan := make(chan int)
	unmappedChan := mapperutils.NewUnmappedChannel()

	go Producer(taskChan)

	var wgFirstPass sync.WaitGroup

	for i := 0; i < 4; i++ {
		wgFirstPass.Add(1)
		go ConsumerFirst(i, taskChan, unmappedChan, &wgFirstPass)
	}

	var wgUnmapped sync.WaitGroup
	wgUnmapped.Add(1)

	fmt.Println("waiting for first pass")

	wgFirstPass.Wait()

	unmappedChan.Close()

	go ConsumerSecond(unmappedChan, &wgUnmapped)

	fmt.Println("waiting for unmapped pass")

	wgUnmapped.Wait()

	fmt.Println("done all")
}

func Producer(taskChan chan<- int) {
	fmt.Println("Producer")

	defer close(taskChan)

	for i := 0; i < 1_000_000_000; i++ {
		taskChan <- i
	}

	fmt.Println("done producing")
}

func ConsumerFirst(id int, taskChan <-chan int, unmappedChan *mapperutils.UnmappedChannel,
	wg *sync.WaitGroup,
) {
	defer wg.Done()

	for task := range taskChan {
		fmt.Println("task: ", task)

		if task%2 == 0 {
			unmappedChan.Send(&mapperutils.UnmappedTask{
				ReadPair: nil,
				ResultFw: nil,
				ResultRv: nil,
			})
		}
	}

	fmt.Println("producer done", id)
}

func ConsumerSecond(unmappedChan *mapperutils.UnmappedChannel, wg *sync.WaitGroup) {
	fmt.Println("ConsumerSecond")

	defer wg.Done()

	for {
		task, ok := unmappedChan.Receive()
		if !ok {
			break
		}
		fmt.Println("unmapped pass: ", task)
	}

	fmt.Println("done unmapped pass")
}
