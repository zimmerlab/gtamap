package main

import (
	"fmt"
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/mapper/incompletemappingpass"
)

func main() {
	fmt.Println("hello")

	taskChan := make(chan int)
	incompleteMappingChan := incompletemappingpass.NewIncompleteMappingChannel()

	go Producer(taskChan)

	var wgFirstPass sync.WaitGroup

	for i := 0; i < 4; i++ {
		wgFirstPass.Add(1)
		go ConsumerFirst(i, taskChan, incompleteMappingChan, &wgFirstPass)
	}

	var wgIncompleteMapping sync.WaitGroup
	wgIncompleteMapping.Add(1)

	fmt.Println("waiting for first pass")

	wgFirstPass.Wait()

	incompleteMappingChan.Close()

	go ConsumerSecond(incompleteMappingChan, &wgIncompleteMapping)

	fmt.Println("waiting for incomplete mapping pass")

	wgIncompleteMapping.Wait()

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

func ConsumerFirst(id int, taskChan <-chan int, incompleteMappingChan *incompletemappingpass.IncompleteMappingChannel,
	wg *sync.WaitGroup,
) {
	defer wg.Done()

	for task := range taskChan {
		fmt.Println("task: ", task)

		if task%2 == 0 {
			incompleteMappingChan.Send(&incompletemappingpass.IncompleteMappingTask{
				ReadPair: nil,
				ResultFw: nil,
				ResultRv: nil,
			})
		}
	}

	fmt.Println("producer done", id)
}

func ConsumerSecond(incompleteMappingChan *incompletemappingpass.IncompleteMappingChannel, wg *sync.WaitGroup) {
	fmt.Println("ConsumerSecond")

	defer wg.Done()

	for {
		task, ok := incompleteMappingChan.Receive()
		if !ok {
			break
		}
		fmt.Println("unmapped pass: ", task)
	}

	fmt.Println("done unmapped pass")
}
