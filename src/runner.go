package main

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"hash/fnv"
	"math"
	"os"
	"runtime"
	"time"

	"github.com/KleinSamuel/gtamap/src/dataloader"
)

func printMemUsage() {
	var m runtime.MemStats
	runtime.ReadMemStats(&m)

	fmt.Printf("Alloc = %v MB\n", bToMb(m.Alloc))
	fmt.Printf("TotalAlloc = %v MB\n", bToMb(m.TotalAlloc))
	fmt.Printf("Sys = %v MB\n", bToMb(m.Sys))
	fmt.Printf("NumGC = %v\n", m.NumGC)
}

func bToMb(b uint64) uint64 {
	return b / 1024 / 1024
}

// BloomFilter represents the bit array and hash functions
type BloomFilter struct {
	bitset []bool
	size   uint
	hashes uint
}

// NewBloomFilter creates a new Bloom filter
func NewBloomFilter(n uint, falsePositiveRate float64) *BloomFilter {
	// Optimal size of bit array m
	m := uint(math.Ceil(-1 * float64(n) * math.Log(falsePositiveRate) / (math.Ln2 * math.Ln2)))
	// Optimal number of hash functions k
	k := uint(math.Ceil((float64(m) / float64(n)) * math.Ln2))

	return &BloomFilter{
		bitset: make([]bool, m),
		size:   m,
		hashes: k,
	}
}

// Add inserts a k-mer (as []byte) into the filter
func (bf *BloomFilter) Add(item []byte) {
	for i := uint(0); i < bf.hashes; i++ {
		idx := bf.hash(item, i) % bf.size
		bf.bitset[idx] = true
	}
}

// Check returns true if the k-mer may exist (false positive possible)
func (bf *BloomFilter) Check(item []byte) bool {
	for i := uint(0); i < bf.hashes; i++ {
		idx := bf.hash(item, i) % bf.size
		if !bf.bitset[idx] {
			return false
		}
	}
	return true
}

// hash function combining FNV with iteration index
func (bf *BloomFilter) hash(item []byte, i uint) uint {
	h := fnv.New64a()
	h.Write(item)
	sum := h.Sum64()
	// Combine with index for multiple hash functions
	return uint(sum + uint64(i*0x9e3779b97f4a7c13))
}

const (
	valueBits = 29
	valueMask = (uint64(1) << valueBits) - 1
)

func pack(index uint16, value uint32) uint64 {
	// Optional safety checks (recommended during development)
	if index >= 500 {
		panic("index overflow")
	}
	if value >= 300_000_000 {
		panic("value overflow")
	}

	return (uint64(index) << valueBits) | uint64(value)
}

func unpackIndex(packed uint64) uint16 {
	return uint16(packed >> valueBits)
}

func unpackValue(packed uint64) uint32 {
	return uint32(packed & valueMask)
}

func writePosToCluster(posToCluster map[int][]uint64, filepath string) error {
	file, err := os.Create(filepath)
	if err != nil {
		return fmt.Errorf("failed to create file: %w", err)
	}
	defer file.Close()

	// gzWriter := gzip.NewWriter(file)
	writer := bufio.NewWriterSize(file, 1<<20) // 1MB buffer
	defer writer.Flush()

	// Write map length
	binary.Write(writer, binary.LittleEndian, uint64(len(posToCluster)))

	// Write each entry
	for key, slice := range posToCluster {
		binary.Write(writer, binary.LittleEndian, int64(key))
		binary.Write(writer, binary.LittleEndian, uint64(len(slice)))
		binary.Write(writer, binary.LittleEndian, slice)
	}

	return nil
}

func writePosToCount(posToCount map[int]uint32, filepath string) error {
	file, err := os.Create(filepath)
	if err != nil {
		return fmt.Errorf("failed to create file: %w", err)
	}
	defer file.Close()

	// gzWriter := gzip.NewWriter(file)
	writer := bufio.NewWriterSize(file, 1<<20) // 1MB buffer
	defer writer.Flush()

	// Write map length
	binary.Write(writer, binary.LittleEndian, uint64(len(posToCount)))

	// Write each entry
	for key, count := range posToCount {
		binary.Write(writer, binary.LittleEndian, int64(key))
		binary.Write(writer, binary.LittleEndian, count)
	}

	return nil
}

func ReadPosToCluster(filepath string) (map[int][]uint64, error) {
	file, err := os.Open(filepath)
	if err != nil {
		return nil, fmt.Errorf("failed to open file: %w", err)
	}
	defer file.Close()

	reader := bufio.NewReaderSize(file, 1<<20)

	var mapLen uint64
	binary.Read(reader, binary.LittleEndian, &mapLen)

	posToCluster := make(map[int][]uint64, mapLen)

	for i := uint64(0); i < mapLen; i++ {
		var key int64
		var sliceLen uint64
		binary.Read(reader, binary.LittleEndian, &key)
		binary.Read(reader, binary.LittleEndian, &sliceLen)

		slice := make([]uint64, sliceLen)
		binary.Read(reader, binary.LittleEndian, slice)
		posToCluster[int(key)] = slice
	}

	return posToCluster, nil
}

func ReadPosToCount(filepath string) (map[int]uint32, error) {
	file, err := os.Open(filepath)
	if err != nil {
		return nil, fmt.Errorf("failed to open file: %w", err)
	}
	defer file.Close()

	reader := bufio.NewReaderSize(file, 1<<20)

	var mapLen uint64
	binary.Read(reader, binary.LittleEndian, &mapLen)

	posToCount := make(map[int]uint32, mapLen)

	for i := uint64(0); i < mapLen; i++ {
		var key int64
		var count uint32
		binary.Read(reader, binary.LittleEndian, &key)
		binary.Read(reader, binary.LittleEndian, &count)
		posToCount[int(key)] = count
	}

	return posToCount, nil
}

func test1() {
	geneFilePath := "/home/sam/Data/ENSG00000198947.fa"
	geneFile, errGeneFile := os.Open(geneFilePath)
	if errGeneFile != nil {
		fmt.Println(errGeneFile)
	}

	// fastaFilePath := "/home/sam/Data/reference-genomes/ensembl/115/homo_sapiens/Homo_sapiens.GRCh38.dna.chromosome.20.fa"
	fastaFilePath := "/home/sam/Data/reference-genomes/ensembl/115/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
	fastaFile, errFile := os.Open(fastaFilePath)
	if errFile != nil {
		fmt.Println(errFile)
	}

	// outFilePath := "/home/sam/Data/test.gtai"
	// outFile, errOutFile := os.Create(outFilePath)
	// if errOutFile != nil {
	// 	fmt.Println(errOutFile)
	// }

	geneEntries, errGeneEntries := dataloader.ReadFasta(geneFile)
	if errGeneEntries != nil {
		fmt.Println(errGeneEntries)
	}

	fastaEntries, err := dataloader.ReadFasta(fastaFile)
	if err != nil {
		fmt.Println(err)
	}

	fmt.Println(len(fastaEntries))
	fmt.Println(len(fastaEntries[0].Sequence))
	fmt.Println(fastaEntries[0].Header)

	geneKmers := make(map[[20]byte][]int)

	for kStart := 0; kStart <= len(geneEntries[0].Sequence)-20; kStart++ {
		seqBytes := *(*[20]byte)((geneEntries[0].Sequence)[kStart : kStart+20])
		geneKmers[seqBytes] = make([]int, 0)
	}

	fmt.Println("gene kmers:", len(geneKmers))

	timerStart := time.Now()

	genomeKmers := make(map[[20]byte]map[uint64]struct{}, len(geneKmers))

	for i := 0; i < len(fastaEntries); i++ {

		fmt.Println("processing chromosome", i, "of", len(fastaEntries), fastaEntries[i].Header)
		printMemUsage()

		for kStart := 0; kStart <= len(fastaEntries[i].Sequence)-20; kStart++ {
			seqBytes := *(*[20]byte)((fastaEntries[i].Sequence)[kStart : kStart+20])

			if _, found := geneKmers[seqBytes]; !found {
				continue
			}

			// var posIdx int = kStart / 500
			posIdx := pack(uint16(i), uint32(kStart/500))

			if _, found := genomeKmers[seqBytes]; !found {
				genomeKmers[seqBytes] = make(map[uint64]struct{})
			}

			genomeKmers[seqBytes][posIdx] = struct{}{}
		}

		if i >= 24 {
			break
		}
	}

	d := time.Since(timerStart)

	fmt.Println("genome kmers:", len(genomeKmers), "in", d)

	// posToCluster := make(map[int][]uint64)
	posToCount := make(map[int]uint32)

	for kStart := 0; kStart < len(geneEntries[0].Sequence)-20; kStart++ {
		seqBytes := *(*[20]byte)((geneEntries[0].Sequence)[kStart : kStart+20])

		positions, found := genomeKmers[seqBytes]
		if !found {
			continue
		}

		// posToCluster[kStart] = make([]uint64, 0)
		posToCount[kStart] = uint32(len(positions))

		// for posIdx := range positions {
		// 	posToCluster[kStart] = append(posToCluster[kStart], posIdx)
		// }

		// for posIdx := range posToCluster[kStart-1] {
		// 	if _, found := posToCluster[kStart][posIdx]; !found {
		// 		delete(posToCluster[kStart-1], posIdx)
		// 	}
		// }

		if kStart%100000 == 0 {
			fmt.Println(kStart, "of", len(geneEntries[0].Sequence))
			printMemUsage()
		}

	}

	// free memory
	genomeKmers = nil

	printMemUsage()

	serPath := "/home/sam/Data/gtamap/evaluation/genes/dmd/test.gob"

	// writePosToCluster(posToCluster, serPath)
	writePosToCount(posToCount, serPath)

	// bf := NewBloomFilter(uint(len(geneKmers)), 0.1)
	//
	// for kmer := range geneKmers {
	// 	bf.Add(kmer[:])
	// }

	// compute clusters
	// numChr := len(fastaEntries)
	// clusterSize := 500

	// fastaEntries[0].Header = "20\tchr20\t+\t1\t64444167"

	// count := 0
	//
	// for kStart := 0; kStart <= len(fastaEntries[0].Sequence)-10; kStart++ {
	// 	seqBytes := *(*[10]byte)((fastaEntries[0].Sequence)[kStart : kStart+10])
	//
	// 	_, found := geneKmers[seqBytes]
	//
	// 	if !found {
	// 		continue
	// 	}
	//
	// 	count++
	//
	// 	// compute cluster index
	// 	var clusterIdx int
	// 	clusterIdx = kStart / clusterSize
	//
	// 	// geneKmers[seqBytes] += 1
	// 	geneKmers[seqBytes] = append(geneKmers[seqBytes], clusterIdx)
	// }
	//
	// fmt.Println(count)

	// genomeIndex := index.BuildGenomeIndex(fastaEntries)
	// index.WriteGenomeIndex(genomeIndex, outFile)
}

func test2() {
	serPath := "/home/sam/Data/gtamap/evaluation/genes/dmd/test.gob"

	posToCluster, err := ReadPosToCluster(serPath)
	if err != nil {
		fmt.Println(err)
	}

	fmt.Println(len(posToCluster))
}

func test3() {
	serPath := "/home/sam/Data/gtamap/evaluation/genes/dmd/genomic-kmers.count.all.gob"

	posToCount, err := ReadPosToCount(serPath)
	if err != nil {
		fmt.Println(err)
	}

	maxVal := 0

	for _, v := range posToCount {
		if int(v) > maxVal {
			maxVal = int(v)
		}
		// fmt.Println(k, v)
		// break
	}

	for i := 0; i <= maxVal; i++ {
		fmt.Println(i, ",", posToCount[i])
	}
}

func main() {
	// test1()
	test3()
}
