package thirdpass

import (
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"

	"github.com/KleinSamuel/gtamap/src/formats/fastq"
)

type ThirdPassTask struct {
	// TODO
	ReadPair *fastq.ReadPair
	ResultFw *mapperutils.ReadMatchResult
	ResultRv *mapperutils.ReadMatchResult
}

type ThirdPassChannel struct {
	in     chan *ThirdPassTask
	out    chan *ThirdPassTask
	buffer []*ThirdPassTask
	mu     sync.Mutex
	closed bool
}

func NewThirdPassChannel() *ThirdPassChannel {
	channel := &ThirdPassChannel{
		in:     make(chan *ThirdPassTask),
		out:    make(chan *ThirdPassTask),
		buffer: make([]*ThirdPassTask, 0),
	}

	go channel.process()

	return channel
}

func (s *ThirdPassChannel) process() {
	var (
		outCh chan<- *ThirdPassTask
		next  *ThirdPassTask
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

func (s *ThirdPassChannel) Send(task *ThirdPassTask) {
	s.in <- task
}

func (s *ThirdPassChannel) Receive() (*ThirdPassTask, bool) {
	item, ok := <-s.out
	return item, ok
}

func (s *ThirdPassChannel) Close() {
	close(s.in)
}
