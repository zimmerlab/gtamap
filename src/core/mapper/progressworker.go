package mapper

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
	"strconv"
	"sync"
	"time"
)

func ProgressWorker(progressChan <-chan bool, wg *sync.WaitGroup) {
	logrus.Debug("Started progressWorker")

	numTasksDone := 0

	timerStart := time.Now()
	timerLast := time.Now()

	for _ = range progressChan {
		numTasksDone++
		if numTasksDone%1_000_000 == 0 {

			timerEnd := time.Now()
			elapsed := timerEnd.Sub(timerStart)
			elapsedSinceLast := timerEnd.Sub(timerLast)
			timerLast = timerEnd

			millionTasks := numTasksDone / 1_000_000
			rate := float64(millionTasks) / elapsed.Minutes()

			logrus.WithFields(logrus.Fields{
				"done":          strconv.Itoa(millionTasks) + "M",
				"delta total":   utils.FormatDuration(elapsed),
				"delta last 1M": utils.FormatDuration(elapsedSinceLast),
				"rate":          fmt.Sprintf("%.2fM per minute", rate),
			}).Info("Progress update")
		}
	}

	wg.Done()

	logrus.Debug("Finished progressWorker")
}
