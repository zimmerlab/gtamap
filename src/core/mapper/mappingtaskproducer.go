package mapper

import (
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"github.com/sirupsen/logrus"
	"strings"
)

func MappingTaskProducer(reader *fastq.Reader, taskChan chan<- MappingTask, progressChan chan<- Event,
	maxReads int, specificQname string) {

	defer close(taskChan)

	// TODO: remove after testing
	taskCounter := 0
	foundSpecificQname := false

	// add each read pair as a mapping task to the task queue
	for readPair, _ := reader.NextRead(); readPair != nil; readPair, _ = reader.NextRead() {

		if specificQname != "" {
			name := strings.Split(readPair.ReadR1.Header, " ")[0]
			if name != specificQname {
				continue
			}
			foundSpecificQname = true
		}

		progressChan <- Event{
			Type: EventBytesProcessed,
			Data: uint64(reader.ProgressReaderR1.BytesRead + reader.ProgressReaderR2.BytesRead),
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

	logrus.Debug("MappingTaskProducer finished")
}
