package datawriter

import (
	"os"
	"path/filepath"

	"github.com/sirupsen/logrus"
)

type Writer struct {
	file *os.File
}

func InitFromPath(outputFilePath string) *Writer {
	err := os.MkdirAll(filepath.Dir(*&outputFilePath), os.ModePerm)
	if err != nil {
		logrus.Fatalf("Error creating output SAM directory: %s", err)
	}

	file, err := os.Create(outputFilePath)
	if err != nil {
		logrus.Fatal("Error creating file", err)
	}

	return InitFromFile(file)
}

func InitFromFile(outputFile *os.File) *Writer {
	return &Writer{
		file: outputFile,
	}
}

func (writer Writer) Write(s string) {
	// fmt.Println("write to file")
	// fmt.Println(s)

	_, err := writer.file.WriteString(s)
	if err != nil {
		logrus.Error(err)
		return
	}

	// fmt.Println("write to file done")
}

func (writer Writer) Close() {
	writer.file.Close()
}
