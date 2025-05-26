package thirdpass

import (
	"sync"

	"github.com/sirupsen/logrus"
)

func ThirdPassWorker(thirdPassChan <-chan *ThirdPassTask, wgThirdPass *sync.WaitGroup) {
	defer wgThirdPass.Done()

	for mappedReadPair := range thirdPassChan {
		logrus.Debugf("Third pass: %s", mappedReadPair.ReadPair.ReadR1.Header)
	}
}
