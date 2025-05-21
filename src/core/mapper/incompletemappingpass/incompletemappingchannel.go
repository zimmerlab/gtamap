package incompletemappingpass

import (
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"

	"github.com/KleinSamuel/gtamap/src/formats/fastq"
)

type IncompleteMappingTask struct {
	ReadPair *fastq.ReadPair
	ResultFw []*mapperutils.ReadMatchResult
	ResultRv []*mapperutils.ReadMatchResult
}

type IncompleteMappingChannel struct {
	in     chan *IncompleteMappingTask
	out    chan *IncompleteMappingTask
	buffer []*IncompleteMappingTask
	mu     sync.Mutex
	closed bool
}

func NewIncompleteMappingChannel() *IncompleteMappingChannel {
	channel := &IncompleteMappingChannel{
		in:     make(chan *IncompleteMappingTask),
		out:    make(chan *IncompleteMappingTask),
		buffer: make([]*IncompleteMappingTask, 0),
	}

	go channel.process()

	return channel
}

func (s *IncompleteMappingChannel) process() {
	var (
		outCh chan<- *IncompleteMappingTask
		next  *IncompleteMappingTask
	)

	for {
		if len(s.buffer) > 0 {
			outCh = s.out
			next = s.buffer[0]
		} else {
			outCh = nil
		}

		select {
		case item, ok := <-s.in:
			if !ok {
				s.mu.Lock()
				s.closed = true
				s.mu.Unlock()

				if len(s.buffer) == 0 {
					close(s.out)
					return
				}

				s.in = nil
				continue
			}

			s.mu.Lock()
			s.buffer = append(s.buffer, item)
			s.mu.Unlock()

		case outCh <- next:
			s.mu.Lock()
			s.buffer = s.buffer[1:]

			if s.closed && len(s.buffer) == 0 {
				close(s.out)
				s.mu.Unlock()
				return
			}

			s.mu.Unlock()
		}
	}
}

func (s *IncompleteMappingChannel) Send(task *IncompleteMappingTask) {
	s.in <- task
}

func (s *IncompleteMappingChannel) Receive() (*IncompleteMappingTask, bool) {
	item, ok := <-s.out
	return item, ok
}

func (s *IncompleteMappingChannel) Close() {
	close(s.in)
}
