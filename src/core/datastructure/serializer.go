package datastructure

import (
	"encoding/gob"
	"github.com/sirupsen/logrus"
	"log"
	"os"
	"time"
)

func SerializeSuffixTreeFromFile(tree *SuffixTree, filePath string) {
	file, err := os.Create(filePath)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	SerializeSuffixTree(tree, file)
}

func SerializeSuffixTree(tree *SuffixTree, outputFile *os.File) {

	timerStart := time.Now()

	enc := gob.NewEncoder(outputFile)
	if err := enc.Encode(tree); err != nil {
		log.Fatal("encode error:", err)
	}

	logrus.WithFields(logrus.Fields{
		"duration": time.Since(timerStart),
		"output":   outputFile.Name(),
	}).Info("Serialized suffix tree.")
}

func DezerializeSuffixTreeFromFile(indexFilePath string) *SuffixTree {
	// Open the file for reading
	file, err := os.Open(indexFilePath)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	return DezerializeSuffixTree(file)
}

func DezerializeSuffixTree(indexFile *os.File) *SuffixTree {

	timerStart := time.Now()

	// Create a decoder and deserialize the person struct from the file
	decoder := gob.NewDecoder(indexFile)
	var tree SuffixTree
	if err := decoder.Decode(&tree); err != nil {
		panic(err)
	}

	logrus.WithFields(logrus.Fields{
		"duration": time.Since(timerStart),
	}).Info("Deserialized suffix tree")

	return &tree
}
