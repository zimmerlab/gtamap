package fastq

import (
	"bufio"
	"fmt"
	"github.com/sirupsen/logrus"
	"os"
)

type Read struct {
	Header   string
	Sequence string
	Quality  string
}

type ReadPair struct {
	ReadR1 *Read
	ReadR2 *Read
}

type Reader struct {
	scannerR1 *bufio.Scanner
	scannerR2 *bufio.Scanner
}

func InitFromPaths(pathR1Reads string, pathR2Reads string) *Reader {

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

	fileR1Reads, errR1 = os.Open(pathR1Reads)
	if errR1 != nil {
		logrus.Fatal("Error reading R1 reads", errR1)
	}
	scannerR1 = bufio.NewScanner(fileR1Reads)

	logrus.Debug("Initialized reader with R1 reads")

	// if R2 reads are given, open them as well
	if len(pathR2Reads) > 0 {
		fmt.Println("Reverse reads given")
		fileR2Reads, errR2 = os.Open(pathR2Reads)
		if errR2 != nil {
			logrus.Fatal("Error reading R2 reads", errR2)
		}
		scannerR2 = bufio.NewScanner(fileR2Reads)

		logrus.Debug("Initialized reader with R2 reads")
	}

	return &Reader{
		scannerR1: scannerR1,
		scannerR2: scannerR2,
	}
}

func (r Reader) NextRead() *ReadPair {

	if !r.scannerR1.Scan() {
		return nil
	}

	header := r.scannerR1.Text()
	r.scannerR1.Scan()
	sequence := r.scannerR1.Text()
	r.scannerR1.Scan()
	r.scannerR1.Scan()
	quality := r.scannerR1.Text()

	var fwRead *Read = &Read{
		Header:   header,
		Sequence: sequence,
		Quality:  quality,
	}
	var rvRead *Read = nil

	if r.scannerR2 != nil {
		if !r.scannerR2.Scan() {
			return nil
		}

		header := r.scannerR2.Text()
		r.scannerR2.Scan()
		sequence := r.scannerR2.Text()
		r.scannerR2.Scan()
		r.scannerR2.Scan()
		quality := r.scannerR2.Text()

		rvRead = &Read{
			Header:   header,
			Sequence: sequence,
			Quality:  quality,
		}
	}

	return &ReadPair{
		ReadR1: fwRead,
		ReadR2: rvRead,
	}
}
