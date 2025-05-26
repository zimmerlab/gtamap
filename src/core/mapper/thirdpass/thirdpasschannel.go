package thirdpass

import (
	"fmt"
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"
)

type ThirdPassTask struct {
	ReadPairId  string
	TargetInfo  *mapperutils.ReadPairMatchResults
	ParalogInfo map[string]map[int]*mapperutils.ValidReadPairCombination
}

type ThirdPassChannel struct {
	in     chan *ThirdPassTask
	out    chan *ThirdPassTask
	buffer map[string]*ThirdPassTask
	mu     sync.Mutex
	closed bool
}

func NewThirdPassChannel() *ThirdPassChannel {
	channel := &ThirdPassChannel{
		in:     make(chan *ThirdPassTask),
		out:    make(chan *ThirdPassTask),
		buffer: make(map[string]*ThirdPassTask),
	}

	go channel.process()

	return channel
}

func (s *ThirdPassChannel) process() {
	for {
		select {
		case task, ok := <-s.in:
			if !ok {
				s.mu.Lock()
				s.closed = true
				if len(s.buffer) == 0 {
					close(s.out)
					s.mu.Unlock()
					return
				}
				s.mu.Unlock()
				s.in = nil
				continue
			}

			s.mu.Lock()
			existing, found := s.buffer[task.ReadPairId]
			if !found {
				s.buffer[task.ReadPairId] = task
			} else {
				if task.TargetInfo != nil {
					existing.TargetInfo = task.TargetInfo
				}
				if task.ParalogInfo != nil {
					existing.ParalogInfo = task.ParalogInfo
				}
			}

			merged := s.buffer[task.ReadPairId]
			if merged.TargetInfo != nil && merged.ParalogInfo != nil {
				s.out <- merged
				delete(s.buffer, task.ReadPairId)
			}
			s.mu.Unlock()
		}
	}
}

func (s *ThirdPassChannel) Receive() (*ThirdPassTask, bool) {
	item, ok := <-s.out
	return item, ok
}

func (s *ThirdPassChannel) Send(task *ThirdPassTask) {
	s.mu.Lock()
	if s.closed {
		s.mu.Unlock()
		fmt.Printf("Cannot send task %s: channel is closed\n", task.ReadPairId)
		return
	}
	s.mu.Unlock()
	s.in <- task
}

func (s *ThirdPassChannel) Close() {
	s.mu.Lock()
	defer s.mu.Unlock()
	if !s.closed {
		close(s.in)
	}
}
