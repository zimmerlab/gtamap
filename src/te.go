package main

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"sync"
)

func main() {

	fmt.Println("hello")

	taskChan := make(chan int)
	secondPassChan := mapperutils.NewFourthPassChannel()

	go Producer(taskChan)

	var wgFirstPass sync.WaitGroup

	for i := 0; i < 4; i++ {
		wgFirstPass.Add(1)
		go ConsumerFirst(i, taskChan, secondPassChan, &wgFirstPass)
	}

	var wgSecondPass sync.WaitGroup
	wgSecondPass.Add(1)

	fmt.Println("waiting for first pass")

	wgFirstPass.Wait()

	secondPassChan.Close()

	go ConsumerSecond(secondPassChan, &wgSecondPass)

	fmt.Println("waiting for second pass")

	wgSecondPass.Wait()

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

func ConsumerFirst(id int, taskChan <-chan int, secondPassChan *mapperutils.FourthPassChannel,
	wg *sync.WaitGroup) {

	defer wg.Done()

	for task := range taskChan {
		fmt.Println("task: ", task)

		if task%2 == 0 {
			secondPassChan.Send(&mapperutils.FourthPassTask{
				ReadPair: nil,
				ResultFw: nil,
				ResultRv: nil,
			})
		}
	}

	fmt.Println("producer done", id)
}

func ConsumerSecond(secondPassChan *mapperutils.FourthPassChannel, wg *sync.WaitGroup) {
	fmt.Println("ConsumerSecond")

	defer wg.Done()

	for {
		task, ok := secondPassChan.Receive()
		if !ok {
			break
		}
		fmt.Println("second pass: ", task)
	}

	fmt.Println("done second pass")
}
