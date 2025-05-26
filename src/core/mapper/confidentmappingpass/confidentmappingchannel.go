package confidentmappingpass

import (
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"

	"github.com/KleinSamuel/gtamap/src/formats/fastq"
)

type ConfidentTask struct {
	ReadPair *fastq.ReadPair
	ResultFw *mapperutils.ReadMatchResult
	ResultRv *mapperutils.ReadMatchResult
}

type ConfidentPassChan struct {
	in     chan *ConfidentTask
	out    chan *ConfidentTask
	buffer []*ConfidentTask
	mu     sync.Mutex
	closed bool
}

func NewConfidentChannel() *ConfidentPassChan {
	channel := &ConfidentPassChan{
		in:     make(chan *ConfidentTask),
		out:    make(chan *ConfidentTask),
		buffer: make([]*ConfidentTask, 0),
	}

	go channel.process()

	return channel
}

func (s *ConfidentPassChan) process() {
	var (
		outCh chan<- *ConfidentTask
		next  *ConfidentTask
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

func (s *ConfidentPassChan) Send(task *ConfidentTask) {
	s.in <- task
}

func (s *ConfidentPassChan) Receive() (*ConfidentTask, bool) {
	item, ok := <-s.out
	return item, ok
}

func (s *ConfidentPassChan) Close() {
	close(s.in)
}
