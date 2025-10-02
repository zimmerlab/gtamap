package events

import "time"

// enum defining event types
const (
	EventTypeReadsProcessed       = 1
	EventTypeReadsAfterFiltering  = 2
	EventTypeReadsMapped          = 3
	EventTypeNumMappingLocations  = 4
	EventTypeNumConfidentMappings = 5
	EventFileSize                 = 6
	EventBytesProcessed           = 7
	EventTypeMapperProducerTime   = 8
	EventTypeMapperWorkerTime     = 9
	EventTypeConfidentWorkerTime  = 10
	EventTypeSecondPassWorkerTime = 11
	EventTypeOutputWorkerTime     = 12
)

type Event struct {
	Type  uint16
	Data  uint64
	Start time.Time // optional
	End   time.Time // optional
}
