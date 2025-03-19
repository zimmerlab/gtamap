package mapper

import (
	"github.com/KleinSamuel/gtamap/src/core/timer"
	"github.com/sirupsen/logrus"
	"sync"
	"time"
)

func TimerWorker(timerChannel <-chan *timer.Timer, wg *sync.WaitGroup) {
	defer wg.Done()

	totalMapR1 := time.Duration(0)
	totalMapR2 := time.Duration(0)
	totalDetermineReadLocation := time.Duration(0)
	totalExactMatch := time.Duration(0)

	for t := range timerChannel {
		totalMapR1 += t.MapR1
		totalMapR2 += t.MapR2
		totalDetermineReadLocation += t.DetermineReadLocation
		totalExactMatch += t.ExactMatch
	}

	logrus.Info("Total mapping time: ", totalMapR1+totalMapR2)
	logrus.Info("Total mapping time R1: ", totalMapR1)
	logrus.Info("Total mapping time R2: ", totalMapR2)
	logrus.Info("Total determine read location time: ", totalDetermineReadLocation)
	logrus.Info("Total exact match time: ", totalExactMatch)
}
