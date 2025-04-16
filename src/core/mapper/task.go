package mapper

import (
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"strings"
)

func startTaskProducer(reader *fastq.Reader, taskQueueMapping chan<- MappingTask, maxReads int, specificQname string) {
	go func() {

		// TODO: remove after testing
		taskCounter := 0
		foundSpecificQname := false

		// add each read pair as a mapping task to the task queue
		for readPair := reader.NextRead(); readPair != nil; readPair = reader.NextRead() {

			if specificQname != "" {
				name := strings.Split(readPair.ReadR1.Header, " ")[0]
				if name != specificQname {
					continue
				}
				foundSpecificQname = true
			}

			mappingTask := MappingTask{
				ID:       taskCounter,
				ReadPair: readPair,
			}
			taskQueueMapping <- mappingTask
			taskCounter++

			// TODO: remove after testing
			if maxReads > 0 && taskCounter >= maxReads {
				break
			}
			if foundSpecificQname {
				break
			}
		}

		// Close the channel when all reads are processed
		close(taskQueueMapping)
	}()
}
