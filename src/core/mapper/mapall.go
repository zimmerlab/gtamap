package mapper

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/core/timer"
	"github.com/KleinSamuel/gtamap/src/datawriter"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
	"github.com/sirupsen/logrus"
	"runtime"
	"sync"
	"time"
)

func MapAll(genomeIndex *index.GenomeIndex, reader *fastq.Reader, writer *datawriter.Writer,
	numThreads *int) {

	runtime.GOMAXPROCS(runtime.NumCPU())

	timerStartTotal := time.Now()

	sequenceInfos := genomeIndex.GetSequenceInfos()

	samHeader := sam.Header{
		Version:                 sam.Version,
		SequenceInfos:           sequenceInfos,
		GenomeAnnotationVersion: "",
		GenomeAssemblyVersion:   "",
		OrganismTaxId:           "",
		ToolVersion:             config.ToolVersion(),
	}

	writer.Write(samHeader.String())

	numWorkers := runtime.NumCPU()

	if *numThreads > 0 {
		numWorkers = *numThreads
	}

	// the size of the task queue buffer
	bufferSizeMultiplier := 4

	// TODO: remove after testing
	maxTasks := 0
	specificQname := ""
	//specificQname := "@A00604:202:HLYW3DSXY:3:1257:26775:6496 2:N:0:AACTCGGA+TCTGGACA"
	//specificQname := "@43"

	// reads to be debugged:
	// @43_fw		needs to be extended to the right

	// contains the read pairs that need to be mapped
	taskChan := make(chan MappingTask, numWorkers*bufferSizeMultiplier)
	// contains all read pairs that could not be mapped in the first pass
	secondPassChan := mapperutils.NewSecondPassChannel()
	// contains the string results of the mapping
	outputChan := make(chan string)
	// contains information about the duration of each step
	timerChan := make(chan *timer.Timer)

	var waitgroupWriter sync.WaitGroup
	waitgroupWriter.Add(1)
	go OutputWorker(outputChan, &waitgroupWriter, writer)

	var waitGroupTimer sync.WaitGroup
	waitGroupTimer.Add(1)
	go TimerWorker(timerChan, &waitGroupTimer)

	// wait group that keeps track of the mapping goroutines that are still running
	var wgFirstPass sync.WaitGroup
	var wgSecondPass sync.WaitGroup

	// start the mapping worker goroutine pool
	for i := 0; i < numWorkers; i++ {
		wgFirstPass.Add(1)
		go MapperWorker(i, genomeIndex, &wgFirstPass, taskChan, secondPassChan, outputChan, timerChan)
	}

	wgSecondPass.Add(1)

	go MappingTaskProducer(reader, taskChan, maxTasks, specificQname)

	wgFirstPass.Wait()

	logrus.Info("Done with first pass mapping")

	secondPassChan.Close()

	go SecondPassWorker(secondPassChan, &wgSecondPass)

	wgSecondPass.Wait()

	close(outputChan)
	close(timerChan)

	waitgroupWriter.Wait()
	waitGroupTimer.Wait()

	writer.Close()

	fmt.Println("mapping finished")

	totalDuration := time.Since(timerStartTotal)

	logrus.WithFields(logrus.Fields{
		"duration": totalDuration,
		"io":       reader.Duration,
	}).Info("Finished mapping")

	fmt.Println("total mapping time: ", totalDuration)
}
