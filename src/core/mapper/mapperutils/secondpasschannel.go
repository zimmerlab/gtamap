package mapperutils

import (
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"sync"
)

type SecondPassTask struct {
	ReadPair *fastq.ReadPair
	ResultFw *[]ReadMatchResult
	ResultRv *[]ReadMatchResult
}

type SecondPassChannel struct {
	in     chan *SecondPassTask
	out    chan *SecondPassTask
	buffer []*SecondPassTask
	mu     sync.Mutex
	closed bool
}

func NewSecondPassChannel() *SecondPassChannel {
	channel := &SecondPassChannel{
		in:     make(chan *SecondPassTask),
		out:    make(chan *SecondPassTask),
		buffer: make([]*SecondPassTask, 0),
	}

	go channel.process()

	return channel
}

func (s *SecondPassChannel) process() {
	var (
		outCh chan<- *SecondPassTask
		next  *SecondPassTask
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

func (s *SecondPassChannel) Send(task *SecondPassTask) {
	s.in <- task
}

func (s *SecondPassChannel) Receive() (*SecondPassTask, bool) {
	item, ok := <-s.out
	return item, ok
}

func (s *SecondPassChannel) Close() {
	close(s.in)
}
