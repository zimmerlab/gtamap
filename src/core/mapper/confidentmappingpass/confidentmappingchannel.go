package confidentmappingpass

import (
	"strconv"
	"strings"
	"sync"

	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mapperutils"

	"github.com/KleinSamuel/gtamap/src/formats/fastq"
)

type ConfidentMappingTask struct {
	ReadPair *fastq.ReadPair
	ResultFw *mapperutils.ReadMatchResult
	ResultRv *mapperutils.ReadMatchResult
	Index    *index.GenomeIndex
}

type ConfidentMappingChannel struct {
	in     chan *ConfidentMappingTask
	out    chan *ConfidentMappingTask
	buffer []*ConfidentMappingTask
	mu     sync.Mutex
	closed bool
}

func NewConfidentMappingChannel() *ConfidentMappingChannel {
	channel := &ConfidentMappingChannel{
		in:     make(chan *ConfidentMappingTask),
		out:    make(chan *ConfidentMappingTask),
		buffer: make([]*ConfidentMappingTask, 0),
	}

	go channel.process()

	return channel
}

func (s *ConfidentMappingChannel) process() {
	var (
		outCh chan<- *ConfidentMappingTask
		next  *ConfidentMappingTask
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

func (s *ConfidentMappingChannel) Send(task *ConfidentMappingTask) {
	s.in <- task
}

func (s *ConfidentMappingChannel) Receive() (*ConfidentMappingTask, bool) {
	item, ok := <-s.out
	return item, ok
}

func (s *ConfidentMappingChannel) Close() {
	close(s.in)
}

func (i ConfidentMappingTask) String() string {
	var builder strings.Builder
	builder.Write([]byte("ReadPairR1 Header: "))
	builder.Write([]byte(i.ReadPair.ReadR1.Header))
	builder.Write([]byte("\n"))
	builder.Write([]byte("  <== FW MAPPING ==>"))
	builder.Write([]byte("\n"))
	builder.Write([]byte("\t SeqIndex: "))
	mapping := i.ResultFw
	seqIndex := strconv.Itoa(mapping.SequenceIndex)
	builder.WriteString(seqIndex)
	builder.WriteString("\n")
	builder.WriteString("\t GENOME -> ")
	builder.WriteString(mapping.MatchedGenome.String())
	builder.WriteString("\n")
	builder.WriteString("\t READ   -> ")
	builder.WriteString(mapping.MatchedRead.String())
	builder.WriteString("\n")
	builder.WriteString("\t MISMAT -> ")
	ints := mapping.MismatchesRead
	strs := make([]string, len(ints))
	for i, v := range ints {
		strs[i] = strconv.Itoa(v)
	}
	builder.WriteString(strings.Join(strs, ","))
	builder.WriteString("\n")
	builder.Write([]byte("  <== RV MAPPING ==>"))
	builder.Write([]byte("\n"))
	builder.Write([]byte("\t SeqIndex: "))
	mapping = i.ResultRv
	seqIndex = strconv.Itoa(mapping.SequenceIndex)
	builder.WriteString(seqIndex)
	builder.WriteString("\n")
	builder.WriteString("\t GENOME -> ")
	builder.WriteString(mapping.MatchedGenome.String())
	builder.WriteString("\n")
	builder.WriteString("\t READ   -> ")
	builder.WriteString(mapping.MatchedRead.String())
	builder.WriteString("\n")
	builder.WriteString("\t MISMAT -> ")
	ints = mapping.MismatchesRead
	strs = make([]string, len(ints))
	for i, v := range ints {
		strs[i] = strconv.Itoa(v)
	}
	builder.WriteString(strings.Join(strs, ","))
	builder.WriteString("\n")
	return builder.String()
}
