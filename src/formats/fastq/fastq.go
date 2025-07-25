package fastq

import (
	"bufio"
	"compress/gzip"
	"github.com/sirupsen/logrus"
	"io"
	"os"
	"path/filepath"
	"strings"
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
	scannerR1        *bufio.Scanner
	scannerR2        *bufio.Scanner
	Duration         *time.Duration
	ProgressReaderR1 *ProgressReader
	ProgressReaderR2 *ProgressReader
}

type ProgressReader struct {
	r          io.Reader
	BytesRead  int64
	TotalBytes int64
}

func (pr *ProgressReader) Read(p []byte) (int, error) {
	n, err := pr.r.Read(p)
	pr.BytesRead += int64(n)
	return n, err
}

func InitFromPaths(pathR1Reads *string, pathR2Reads *string) (*Reader, error) {

	var fileR1Reads *os.File = nil
	var errR1 error = nil

	var fileR2Reads *os.File = nil
	var errR2 error = nil

	fileR1Reads, errR1 = os.Open(*pathR1Reads)
	if errR1 != nil {
		logrus.Fatal("Error reading R1 reads: ", errR1)
	}

	// if R2 reads are given, open them as well
	if pathR2Reads != nil {
		fileR2Reads, errR2 = os.Open(*pathR2Reads)
		if errR2 != nil {
			logrus.Fatal("Error reading R2 reads: ", errR2)
		}
	}

	return InitFromFiles(fileR1Reads, fileR2Reads)
}

func InitFromFiles(fastqR1File *os.File, fastqR2File *os.File) (*Reader, error) {

	logrus.WithFields(logrus.Fields{
		"R1": fastqR1File.Name(),
		"R2": fastqR2File.Name(),
	}).Debug("Initializing reader")

	var progressReaderR1 *ProgressReader = nil
	var scannerR1 *bufio.Scanner = nil
	var progressReaderR2 *ProgressReader = nil
	var scannerR2 *bufio.Scanner = nil

	statR1, err := fastqR1File.Stat()
	if err != nil {
		return nil, err
	}

	progressReaderR1 = &ProgressReader{
		r:          fastqR1File,
		TotalBytes: statR1.Size(),
	}

	if strings.ToLower(filepath.Ext(fastqR1File.Name())) == ".gz" {
		gzipReader, errGzip := gzip.NewReader(progressReaderR1)
		if errGzip != nil {
			logrus.Fatal("Error reading gzipped R1 reads: ", errGzip)
		}
		scannerR1 = bufio.NewScanner(gzipReader)
	} else {
		scannerR1 = bufio.NewScanner(progressReaderR1)
	}

	logrus.Debug("Initialized reader with R1 reads")

	// if R2 reads are given, open them as well
	if fastqR2File != nil {
		logrus.Debug("Reverse reads given")

		statR2, err := fastqR2File.Stat()
		if err != nil {
			return nil, err
		}

		progressReaderR2 = &ProgressReader{
			r:          fastqR2File,
			TotalBytes: statR2.Size(),
		}

		if strings.ToLower(filepath.Ext(fastqR2File.Name())) == ".gz" {
			gzipReader, errGzip := gzip.NewReader(progressReaderR2)
			if errGzip != nil {
				logrus.Fatal("Error reading gzipped R2 reads: ", errGzip)
			}
			scannerR2 = bufio.NewScanner(gzipReader)
		} else {
			scannerR2 = bufio.NewScanner(progressReaderR2)
		}
		logrus.Debug("Initialized reader with R2 reads")
	}

	d := time.Duration(0)

	return &Reader{
		scannerR1:        scannerR1,
		scannerR2:        scannerR2,
		Duration:         &d,
		ProgressReaderR1: progressReaderR1,
		ProgressReaderR2: progressReaderR2,
	}, nil
}

func (r Reader) NextRead() (*ReadPair, uint64) {

	timerStart := time.Now()

	if !r.scannerR1.Scan() {
		return nil, 0
	}

	var bytesRead uint64 = 0

	// read the header, sequence and quality of the read R1
	// scanner.Text() returns a copy of the string but .Bytes() returns a reference
	// to the scanners internal buffer which means it has to be copied to a new slice
	// to avoid clashes when reading the next line

	// remove the leading '@' character from the header
	headerBytesR1 := r.scannerR1.Bytes()
	bytesRead += uint64(len(headerBytesR1))
	headerR1 := string(headerBytesR1)[1:]

	r.scannerR1.Scan()
	sequenceR1Bytes := r.scannerR1.Bytes()
	bytesRead += uint64(len(sequenceR1Bytes))
	sequenceR1 := make([]byte, len(sequenceR1Bytes))
	copy(sequenceR1, sequenceR1Bytes)

	r.scannerR1.Scan()
	bytesRead += uint64(len(r.scannerR1.Bytes()))

	r.scannerR1.Scan()
	qualityR1Bytes := r.scannerR1.Bytes()
	bytesRead += uint64(len(qualityR1Bytes))
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
			return nil, 0
		}

		// remove the leading '@' character from the header
		headerBytesR2 := r.scannerR2.Bytes()
		bytesRead += uint64(len(headerBytesR2))
		headerR2 := string(headerBytesR2)[1:]

		r.scannerR2.Scan()
		sequenceR2Bytes := r.scannerR2.Bytes()
		bytesRead += uint64(len(sequenceR2Bytes))
		sequenceR2 := make([]byte, len(sequenceR2Bytes))
		copy(sequenceR2, sequenceR2Bytes)

		r.scannerR2.Scan()
		bytesRead += uint64(len(r.scannerR2.Bytes()))

		r.scannerR2.Scan()
		qualityR2Bytes := r.scannerR2.Bytes()
		bytesRead += uint64(len(qualityR2Bytes))
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
	}, bytesRead
}
