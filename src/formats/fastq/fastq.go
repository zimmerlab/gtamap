package fastq

import (
	"bufio"
	"github.com/sirupsen/logrus"
	"os"
	"time"
)

type Read struct {
	Header   string
	Sequence *[]byte
	Quality  *[]byte
}

type ReadPair struct {
	ReadR1 *Read
	ReadR2 *Read
}

type Reader struct {
	scannerR1 *bufio.Scanner
	scannerR2 *bufio.Scanner
	Duration  *time.Duration
}

func InitFromPaths(pathR1Reads *string, pathR2Reads *string) *Reader {

	logrus.WithFields(logrus.Fields{
		"R1": pathR1Reads,
		"R2": pathR2Reads,
	}).Debug("Initializing reader")

	var fileR1Reads *os.File = nil
	var errR1 error = nil
	var scannerR1 *bufio.Scanner = nil

	var fileR2Reads *os.File = nil
	var errR2 error = nil
	var scannerR2 *bufio.Scanner = nil

	fileR1Reads, errR1 = os.Open(*pathR1Reads)
	if errR1 != nil {
		logrus.Fatal("Error reading R1 reads: ", errR1)
	}
	scannerR1 = bufio.NewScanner(fileR1Reads)

	logrus.Debug("Initialized reader with R1 reads")

	// if R2 reads are given, open them as well
	if pathR2Reads != nil {
		logrus.Debug("Reverse reads given")
		fileR2Reads, errR2 = os.Open(*pathR2Reads)
		if errR2 != nil {
			logrus.Fatal("Error reading R2 reads: ", errR2)
		}
		scannerR2 = bufio.NewScanner(fileR2Reads)

		logrus.Debug("Initialized reader with R2 reads")
	}

	d := time.Duration(0)

	return &Reader{
		scannerR1: scannerR1,
		scannerR2: scannerR2,
		Duration:  &d,
	}
}

func (r Reader) NextRead() *ReadPair {

	timerStart := time.Now()

	if !r.scannerR1.Scan() {
		return nil
	}

	// read the header, sequence and quality of the read R1
	// scanner.Text() returns a copy of the string but .Bytes() returns a reference
	// to the scanners internal buffer which means it has to be copied to a new slice
	// to avoid clashes when reading the next line
	headerR1 := r.scannerR1.Text()
	r.scannerR1.Scan()
	sequenceR1Bytes := r.scannerR1.Bytes()
	sequenceR1 := make([]byte, len(sequenceR1Bytes))
	copy(sequenceR1, sequenceR1Bytes)
	r.scannerR1.Scan()
	r.scannerR1.Scan()
	qualityR1Bytes := r.scannerR1.Bytes()
	qualityR1 := make([]byte, len(qualityR1Bytes))
	copy(qualityR1, qualityR1Bytes)

	var fwRead *Read = &Read{
		Header:   headerR1,
		Sequence: &sequenceR1,
		Quality:  &qualityR1,
	}
	var rvRead *Read = nil

	if r.scannerR2 != nil {
		if !r.scannerR2.Scan() {
			return nil
		}

		headerR2 := r.scannerR2.Text()
		r.scannerR2.Scan()
		sequenceR2Bytes := r.scannerR2.Bytes()
		sequenceR2 := make([]byte, len(sequenceR2Bytes))
		copy(sequenceR2, sequenceR2Bytes)
		r.scannerR2.Scan()
		r.scannerR2.Scan()
		qualityR2Bytes := r.scannerR2.Bytes()
		qualityR2 := make([]byte, len(qualityR2Bytes))
		copy(qualityR2, qualityR2Bytes)

		rvRead = &Read{
			Header:   headerR2,
			Sequence: &sequenceR2,
			Quality:  &qualityR2,
		}
	}

	*r.Duration += time.Since(timerStart)

	return &ReadPair{
		ReadR1: fwRead,
		ReadR2: rvRead,
	}
}
