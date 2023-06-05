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
	FwRead *Read
	RvRead *Read
}

type Reader struct {
	scannerFw *bufio.Scanner
	scannerRv *bufio.Scanner
}

func InitFromPaths(pathFwReads string, pathRvReads string) *Reader {

	logrus.WithFields(logrus.Fields{
		"forward": pathFwReads,
		"reverse": pathRvReads,
	}).Debug("Initializing reader")

	var fileFwReads *os.File = nil
	var errFw error = nil
	var scannerFw *bufio.Scanner = nil

	var fileRvReads *os.File = nil
	var errRv error = nil
	var scannerRv *bufio.Scanner = nil

	fileFwReads, errFw = os.Open(pathFwReads)
	if errFw != nil {
		logrus.Fatal("Error reading forward reads", errFw)
	}
	scannerFw = bufio.NewScanner(fileFwReads)

	logrus.Debug("Initialized reader with forward reads")

	// If reverse reads are given, open them as well
	if len(pathRvReads) > 0 {
		fmt.Println("Reverse reads given")
		fileRvReads, errRv = os.Open(pathRvReads)
		if errRv != nil {
			logrus.Fatal("Error reading reverse reads", errRv)
		}
		scannerRv = bufio.NewScanner(fileRvReads)

		logrus.Debug("Initialized reader with reverse reads")
	}

	return &Reader{
		scannerFw: scannerFw,
		scannerRv: scannerRv,
	}
}

func (r Reader) NextRead() *ReadPair {

	fmt.Println("Next read")

	if !r.scannerFw.Scan() {
		return nil
	}

	header := r.scannerFw.Text()
	r.scannerFw.Scan()
	sequence := r.scannerFw.Text()
	r.scannerFw.Scan()
	r.scannerFw.Scan()
	quality := r.scannerFw.Text()

	var fwRead *Read = &Read{
		Header:   header,
		Sequence: sequence,
		Quality:  quality,
	}
	var rvRead *Read = nil

	if r.scannerRv != nil {
		if !r.scannerRv.Scan() {
			return nil
		}

		header := r.scannerRv.Text()
		r.scannerRv.Scan()
		sequence := r.scannerRv.Text()
		r.scannerRv.Scan()
		r.scannerRv.Scan()
		quality := r.scannerRv.Text()

		rvRead = &Read{
			Header:   header,
			Sequence: sequence,
			Quality:  quality,
		}
	}

	return &ReadPair{
		FwRead: fwRead,
		RvRead: rvRead,
	}
}
