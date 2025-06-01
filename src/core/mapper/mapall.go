package mapper

import (
	"runtime"
	"sync"
	"time"

	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/confidentmappingpass"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"github.com/KleinSamuel/gtamap/src/core/mapper/secondpass"
	"github.com/KleinSamuel/gtamap/src/core/mapper/thirdpass"
	"github.com/KleinSamuel/gtamap/src/core/timer"
	"github.com/KleinSamuel/gtamap/src/datawriter"
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
	"github.com/KleinSamuel/gtamap/src/utils"
	"github.com/sirupsen/logrus"
)

func MapAll(genomeIndex *index.GenomeIndex, reader *fastq.Reader, writer *datawriter.Writer,
	numThreads *int,
) {
	runtime.GOMAXPROCS(runtime.NumCPU())

	timerStartTotal := time.Now()

	samHeader := sam.Header{
		Version:       sam.Version,
		SequenceInfos: genomeIndex.GetSequenceInfos(),
		ToolVersion:   config.ToolVersion(),
	}

	writer.Write(samHeader.String())

	numWorkers := runtime.NumCPU()
	// TODO: REMOVE DEBUG
	// os.Create("output.txt")

	if *numThreads > 0 {
		numWorkers = *numThreads
	}

	// the size of the task queue buffer
	bufferSizeMultiplier := 10000

	// TODO: remove after testing
	maxTasks := 0

	// fw with left normalization required
	// specificQname := "A00604:202:HLYW3DSXY:3:2114:17020:15562"
	// rv with left normalization required
	// specificQname := "A00604:202:HLYW3DSXY:3:2169:25527:10316"
	// adding region where start == end
	// specificQname := "A00604:202:HLYW3DSXY:3:1103:32217:18991"
	// cigar string contains 1D
	// specificQname := "A00604:202:HLYW3DSXY:3:1157:15121:15468"

	specificQname := ""

	// reads to be debugged:
	// @43_fw		needs to be extended to the right

	// contains the read pairs that need to be mapped
	taskChan := make(chan MappingTask, numWorkers*bufferSizeMultiplier)
	// contains all read pairs (multimapped and parially mapped) except for confident read pairs
	secondpassChan := secondpass.NewSecondPassChannel()
	// contains all read pairs (multimapped and parially mapped) except for confident read pairs
	thirdpassChan := thirdpass.NewThirdPassChannel()
	// contains all confident read pairs that map unique and almost perfect
	confidentMappingChan := confidentmappingpass.NewConfidentChannel()
	// contains annotation per target region
	annotationChan := make(chan map[int]*mapperutils.TargetAnnotation)
	// contains all read pairs that need to be mapped to the paralog index (all)
	paralogMappingChan := make(chan *mapperutils.ReadPairMatchResults)
	// contains the string results of the mapping
	outputChan := make(chan string)
	// contains information about the duration of each step
	timerChan := make(chan *timer.Timer)
	// contains information about the progress of the mapping
	progressChan := make(chan bool)

	var waitgroupProgress sync.WaitGroup
	waitgroupProgress.Add(1)
	go ProgressWorker(progressChan, &waitgroupProgress)

	var waitgroupWriter sync.WaitGroup
	waitgroupWriter.Add(1)
	go OutputWorker(outputChan, &waitgroupWriter, writer)

	var waitgroupParalog sync.WaitGroup
	waitgroupParalog.Add(1)
	go ParalogMappingWorker(paralogMappingChan, &waitgroupParalog, genomeIndex, thirdpassChan)

	var waitGroupTimer sync.WaitGroup
	waitGroupTimer.Add(1)
	go TimerWorker(timerChan, &waitGroupTimer)

	// wait group that keeps track of the mapping goroutines that are still running

	var wgMainMappingPass sync.WaitGroup
	// start the mapping worker goroutine pool
	for i := 0; i < numWorkers; i++ {
		wgMainMappingPass.Add(1)
		go MapperWorker(i, genomeIndex, &wgMainMappingPass, taskChan, secondpassChan, confidentMappingChan, paralogMappingChan, progressChan, timerChan)
	}

	go MappingTaskProducer(reader, taskChan, maxTasks, specificQname)

	wgMainMappingPass.Wait()
	close(paralogMappingChan)

	confidentMappingChan.Close()
	var waitgroupConfidentMap sync.WaitGroup
	waitgroupConfidentMap.Add(1)
	go confidentmappingpass.ConfidentMappingWorker(confidentMappingChan, &waitgroupConfidentMap, annotationChan, genomeIndex)

	var wgThirdPass sync.WaitGroup
	wgThirdPass.Add(1)
	go thirdpass.ThirdPassWorker(thirdpassChan, &wgThirdPass, outputChan, genomeIndex)

	var wgSecondpass sync.WaitGroup
	wgSecondpass.Add(1)
	go secondpass.SecondpassMappingWorker(secondpassChan, &wgSecondpass, annotationChan, thirdpassChan)

	waitgroupConfidentMap.Wait()
	close(annotationChan)

	waitgroupParalog.Wait()
	secondpassChan.Close()
	wgSecondpass.Wait()

	thirdpassChan.Close()
	wgThirdPass.Wait()

	close(outputChan)
	waitgroupWriter.Wait()

	close(progressChan)
	waitgroupProgress.Wait()

	close(timerChan)
	waitGroupTimer.Wait()

	writer.Close()

	totalDuration := time.Since(timerStartTotal)

	logrus.WithFields(logrus.Fields{
		"duration": utils.FormatDuration(totalDuration),
		"io":       utils.FormatDuration(*reader.Duration),
	}).Info("Finished mapping")
}
