package secondpass

import (
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
)

type SecondPassChannel struct {
	in     chan *mapperutils.ReadPairMatchResults
	out    chan *mapperutils.ReadPairMatchResults
	buffer []*mapperutils.ReadPairMatchResults
	mu     sync.Mutex
	closed bool
}

func NewSecondPassChannel() *SecondPassChannel {
	channel := &SecondPassChannel{
		in:     make(chan *mapperutils.ReadPairMatchResults),
		out:    make(chan *mapperutils.ReadPairMatchResults),
		buffer: make([]*mapperutils.ReadPairMatchResults, 0),
	}

	go channel.process()

	return channel
}

func (s *SecondPassChannel) process() {
	var (
		outCh chan<- *mapperutils.ReadPairMatchResults
		next  *mapperutils.ReadPairMatchResults
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

func (s *SecondPassChannel) Send(task *mapperutils.ReadPairMatchResults) {
	s.in <- task
}

func (s *SecondPassChannel) Receive() (*mapperutils.ReadPairMatchResults, bool) {
	item, ok := <-s.out
	return item, ok
}

func (s *SecondPassChannel) Close() {
	s.mu.Lock()
	defer s.mu.Unlock()
	if !s.closed {
		close(s.in)
	}
}
