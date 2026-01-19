package mapper

import (
	"strings"

	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/sirupsen/logrus"
)

func MappingTaskProducer(
	reader *fastq.Reader,
	taskChan chan<- MappingTask,
	progressChan chan<- Event,
	maxReads int,
	specificQname string,
) {
	defer close(taskChan)

	var byteProgress int64 = 0
	// chunk size for progress updates (1e6 = 1MB)
	progressChunkSize := int64(1e6)

	// TODO: remove after testing
	taskCounter := 0
	foundSpecificQname := false

	var startProducer time.Time
	firstItem := true

	// add each read pair as a mapping task to the task queue
	for readPair, _ := reader.NextRead(); readPair != nil; readPair, _ = reader.NextRead() {
		if firstItem {
			startProducer = time.Now()
			firstItem = false
		}

		if specificQname != "" {
			name := strings.Split(readPair.ReadR1.Header, " ")[0]
			if name != specificQname {
				continue
			}
			foundSpecificQname = true
		}

		byteProgress += reader.ProgressReaderR1.BytesRead + reader.ProgressReaderR2.BytesRead

		// send progress event if enough bytes have been processed
		if byteProgress >= progressChunkSize {
			progressChan <- events.Event{
				Type: events.EventBytesProcessed,
				Data: uint64(byteProgress),
			}
			byteProgress = 0
		}

		mappingTask := MappingTask{
			ID:       taskCounter,
			ReadPair: readPair,
		}
		taskChan <- mappingTask
		taskCounter++

		// TODO: remove after testing
		if maxReads > 0 && taskCounter >= maxReads {
			break
		}
		if foundSpecificQname {
			break
		}
	}
	durationProducer := time.Since(startProducer)
	endProducer := time.Now()
	progressChan <- events.Event{
		Type:  events.EventTypeMapperProducerTime,
		Data:  uint64(durationProducer),
		Start: startProducer,
		End:   endProducer,
	}
	logrus.Debug("MappingTaskProducer finished")
}
