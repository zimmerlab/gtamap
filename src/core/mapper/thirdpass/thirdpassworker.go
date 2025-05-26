package thirdpass

import (
	"sync"

	"github.com/sirupsen/logrus"
)

func ThirdPassWorker(thirdPassChan *ThirdPassChannel, wgThirdPass *sync.WaitGroup) {
	defer wgThirdPass.Done()

	for {
		task, ok := thirdPassChan.Receive()
		if !ok {
			break
		}

		logrus.Debugf("Third pass: %s", task)
	}
}
