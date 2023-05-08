package datastructure

import (
	"encoding/gob"
	"github.com/sirupsen/logrus"
	"log"
	"os"
	"time"
)

func SerializeSuffixTree(tree *SuffixTree, filePath string) {

	timerStart := time.Now()

	file, err := os.Create(filePath)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	enc := gob.NewEncoder(file)
	if err := enc.Encode(tree); err != nil {
		log.Fatal("encode error:", err)
	}

	logrus.WithFields(logrus.Fields{
		"duration": time.Since(timerStart),
	}).Info("serialized suffix tree")
}

func DezerializeSuffixTree(filePath string) *SuffixTree {

	timerStart := time.Now()

	// Open the file for reading
	file, err := os.Open(filePath)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	// Create a decoder and deserialize the person struct from the file
	decoder := gob.NewDecoder(file)
	var tree SuffixTree
	if err := decoder.Decode(&tree); err != nil {
		panic(err)
	}

	logrus.WithFields(logrus.Fields{
		"duration": time.Since(timerStart),
	}).Info("deserialized suffix tree")

	return &tree
}
