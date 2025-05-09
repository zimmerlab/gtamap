package mapperutils

import (
	"github.com/KleinSamuel/gtamap/src/formats/fastq"
	"sync"
)

type FourthPassTask struct {
	ReadPair *fastq.ReadPair
	ResultFw *[]ReadMatchResult
	ResultRv *[]ReadMatchResult
}

type FourthPassChannel struct {
	in     chan *FourthPassTask
	out    chan *FourthPassTask
	buffer []*FourthPassTask
	mu     sync.Mutex
	closed bool
}

func NewFourthPassChannel() *FourthPassChannel {
	channel := &FourthPassChannel{
		in:     make(chan *FourthPassTask),
		out:    make(chan *FourthPassTask),
		buffer: make([]*FourthPassTask, 0),
	}

	go channel.process()

	return channel
}

func (s *FourthPassChannel) process() {
	var (
		outCh chan<- *FourthPassTask
		next  *FourthPassTask
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

func (s *FourthPassChannel) Send(task *FourthPassTask) {
	s.in <- task
}

func (s *FourthPassChannel) Receive() (*FourthPassTask, bool) {
	item, ok := <-s.out
	return item, ok
}

func (s *FourthPassChannel) Close() {
	close(s.in)
}
