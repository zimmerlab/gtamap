package unmappedpass

import (
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
	"sync"

	"github.com/KleinSamuel/gtamap/src/formats/fastq"
)

type UnmappedTask struct {
	ReadPair *fastq.ReadPair
	ResultFw []*mapperutils.ReadMatchResult
	ResultRv []*mapperutils.ReadMatchResult
}

type UnmappedChannel struct {
	in     chan *UnmappedTask
	out    chan *UnmappedTask
	buffer []*UnmappedTask
	mu     sync.Mutex
	closed bool
}

func NewUnmappedChannel() *UnmappedChannel {
	channel := &UnmappedChannel{
		in:     make(chan *UnmappedTask),
		out:    make(chan *UnmappedTask),
		buffer: make([]*UnmappedTask, 0),
	}

	go channel.process()

	return channel
}

func (s *UnmappedChannel) process() {
	var (
		outCh chan<- *UnmappedTask
		next  *UnmappedTask
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

func (s *UnmappedChannel) Send(task *UnmappedTask) {
	s.in <- task
}

func (s *UnmappedChannel) Receive() (*UnmappedTask, bool) {
	item, ok := <-s.out
	return item, ok
}

func (s *UnmappedChannel) Close() {
	close(s.in)
}
