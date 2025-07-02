package server

import (
	"encoding/json"
	"fmt"
	"github.com/KleinSamuel/gtamap/src/analysis"
	"github.com/KleinSamuel/gtamap/src/formats/fasta"
	"github.com/gorilla/mux"
	"github.com/rs/cors"
	"log"
	"math"
	"net/http"
	"os"
	"strconv"
	"strings"
	"time"
)

type Server struct {
	Router          *mux.Router
	AnalysisService *analysis.AnalysisService
}

func (s *Server) InitRoutes() {

	api := s.Router.PathPrefix("/api").Subrouter()

	api.HandleFunc("/test", s.test).Methods("GET")

	api.HandleFunc("/health", s.healthCheck).Methods("GET")

	api.HandleFunc("/genomeConfig", getTestGenomeConfig).Methods("GET")

	api.HandleFunc("/readSummaryTable", getReadInfo).Methods("GET")

	api.HandleFunc("/upsetDataRead", s.upsetDataRead).Methods("GET")
	api.HandleFunc("/upsetDataRecordPos", s.upsetDataRecordPos).Methods("GET")

	api.HandleFunc("/readNumMatchesDistribution", s.readNumMatchesDistribution).Methods("GET")
	api.HandleFunc("/readMultimappingDensityData", s.readMultimappingDensityData).Methods("GET")
	api.HandleFunc("/readMultimappingDensityDataHeatmap", s.readMultimappingDensityDataHeatmap).Methods("GET")
	api.HandleFunc("/mapperEmbedding", s.mapperEmbedding).Methods("GET")
	api.HandleFunc("/mapperMultimappingParallel", s.mapperMultimappingParallel).Methods("GET")

	download := s.Router.PathPrefix("/download").Subrouter()

	download.HandleFunc("/genome/fasta", serveGenomeFastaFile).Methods("GET")
	download.HandleFunc("/genome/fastaIndex", serveGenomeFastaIndexFile).Methods("GET")
	download.HandleFunc("/genome/gtf", serveGenomeAnnotationFile).Methods("GET")
	download.HandleFunc("/gtamap/bam", serveBamFile).Methods("GET")
	download.HandleFunc("/gtamap/bamIndex", serveBamIndexFile).Methods("GET")
}

func (s *Server) Start() {
	port := ":8000"

	c := enableCORS()
	handler := c.Handler(s.Router)

	fmt.Printf("Server starting on port %s\n", port)
	log.Fatal(http.ListenAndServe(port, handler))
}

func NewServer() *Server {
	server := &Server{
		Router:          mux.NewRouter(),
		AnalysisService: analysis.NewAnalysisService(),
	}

	server.InitRoutes()

	return server
}

func enableCORS() *cors.Cors {
	return cors.New(cors.Options{
		AllowedOrigins: []string{"http://localhost:8080", "http://localhost:5173"},
		AllowedMethods: []string{"GET", "POST", "PUT", "DELETE", "OPTIONS"},
		AllowedHeaders: []string{"*"},
	})
}

func (s *Server) healthCheck(w http.ResponseWriter, r *http.Request) {
	response := map[string]string{
		"status":    "healthy",
		"timestamp": time.Now().Format(time.RFC3339),
	}
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(response)
}

func (s *Server) test(w http.ResponseWriter, r *http.Request) {

	response := map[string]string{
		"message":  "This is a test response from the server",
		"gtamap":   strconv.Itoa(len(s.AnalysisService.MapperInfos["gtamap"].ParsedFile.Records)),
		"minimap2": strconv.Itoa(len(s.AnalysisService.MapperInfos["minimap2"].ParsedFile.Records)),
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(response)
}

type UpsetElement struct {
	Name string   `json:"name"`
	Sets []string `json:"sets"`
}

func (s *Server) upsetDataRead(w http.ResponseWriter, r *http.Request) {

	data := make([]UpsetElement, 0)

	qnames := make([]string, 0, len(s.AnalysisService.MapperInfos["gtamap"].RecordsByQname))
	for qname := range s.AnalysisService.MapperInfos["gtamap"].RecordsByQname {
		qnames = append(qnames, qname)
	}

	for _, qname := range qnames {

		sets := make([]string, 0)
		for mapperName, mapperInfo := range s.AnalysisService.MapperInfos {
			if _, exists := mapperInfo.RecordsByQname[qname]; exists {
				sets = append(sets, mapperName)
			}
		}

		elem := UpsetElement{
			Name: qname,
			Sets: sets,
		}

		data = append(data, elem)
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(data)
}

func (s *Server) upsetDataRecordPos(w http.ResponseWriter, r *http.Request) {

	data := make([]UpsetElement, 0)

	recordPosMap := make(map[string]map[string]bool)

	for mapperName, mapperInfo := range s.AnalysisService.MapperInfos {

		for _, record := range mapperInfo.ParsedFile.Records {

			id := record.Qname + "X" + strconv.Itoa(record.Pos)

			if _, exists := recordPosMap[id]; !exists {
				//recordPosMap[id] = make([]string, 0)
				recordPosMap[id] = make(map[string]bool)
			}

			//recordPosMap[id] = append(recordPosMap[id], mapperName)
			recordPosMap[id][mapperName] = true
		}
	}

	for id, mappers := range recordPosMap {

		mapperList := make([]string, 0, len(mappers))
		for mapperName := range mappers {
			mapperList = append(mapperList, mapperName)
		}

		elem := UpsetElement{
			Name: id,
			Sets: mapperList,
		}

		data = append(data, elem)
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(data)
}

func (s *Server) readNumMatchesDistribution(w http.ResponseWriter, r *http.Request) {

	// TODO: currently only gtamap

	data := s.AnalysisService.MapperInfos["gtamap"].RecordsByQname

	recordInfoList := make([]RecordInfo, 0, len(data))

	for qname, records := range data {
		recordIds := make([]int, 0, len(records))
		for _, record := range records {
			recordIds = append(recordIds, record.Id)
		}

		recordInfo := RecordInfo{
			Qname:     qname,
			RecordIds: recordIds,
		}

		recordInfoList = append(recordInfoList, recordInfo)
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(recordInfoList)
}

type MultimappingDensity struct {
	MaxXValue int                       `json:"maxXValue"`
	Data      []MultimappingDensityInfo `json:"data"`
}

type MultimappingDensityInfo struct {
	Mapper          string `json:"mapper"`
	ReadAssignments []int  `json:"readAssignments"`
}

type MultimappingDensityInfoHeatmap struct {
	Mapper      string `json:"mapper"`
	NumMappings int    `json:"numMappings"`
	Count       int    `json:"count"`
}

func (s *Server) readMultimappingDensityDataHeatmap(w http.ResponseWriter, r *http.Request) {

	data := make([]MultimappingDensityInfoHeatmap, 0, len(s.AnalysisService.MapperNames))

	maxValue := 0

	mapperToMap := make(map[string]map[int]int)

	for _, mapperName := range s.AnalysisService.MapperNames {

		numMappingsToCount := make(map[int]int)

		for _, records := range s.AnalysisService.MapperInfos[mapperName].RecordsByQname {

			numMappings := len(records)

			if numMappings > maxValue {
				maxValue = numMappings
			}

			if _, exists := numMappingsToCount[numMappings]; !exists {
				numMappingsToCount[numMappings] = 0
			}

			numMappingsToCount[numMappings]++
		}

		if _, exists := mapperToMap[mapperName]; !exists {
			mapperToMap[mapperName] = make(map[int]int)
		}
		mapperToMap[mapperName] = numMappingsToCount
	}

	for mapperName, numMappingsToCount := range mapperToMap {
		for i := 1; i <= maxValue; i++ {
			count, exists := numMappingsToCount[i]
			if !exists {
				count = 0
			}

			multimappingDensityInfo := MultimappingDensityInfoHeatmap{
				Mapper:      mapperName,
				NumMappings: i,
				Count:       count,
			}

			data = append(data, multimappingDensityInfo)
		}
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(data)
}

func (s *Server) readMultimappingDensityData(w http.ResponseWriter, r *http.Request) {

	maxXValue := 0

	data := make([]MultimappingDensityInfo, 0)

	for mapperName, mapperInfo := range s.AnalysisService.MapperInfos {

		readAssignments := make([]int, 0, len(mapperInfo.RecordsByQname))

		for _, records := range mapperInfo.RecordsByQname {
			readAssignments = append(readAssignments, len(records))

			if len(records) > int(maxXValue) {
				maxXValue = len(records)
			}
		}

		multimappingDensityInfo := MultimappingDensityInfo{
			Mapper:          mapperName,
			ReadAssignments: readAssignments,
		}

		data = append(data, multimappingDensityInfo)
	}

	multimappingDensity := MultimappingDensity{
		MaxXValue: maxXValue,
		Data:      data,
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(multimappingDensity)
}

type RecordInfo struct {
	Qname     string `json:"qname"`
	RecordIds []int  `json:"recordIds"`
}

type MapperEmbedding struct {
	Mapper    string      `json:"mapper"`
	Embedding [][]float32 `json:"embedding"`
}

type MapperDistanceMatrix struct {
	Mappers   []string    `json:"mappers"`
	Distances [][]float32 `json:"distances"`
}

func (s *Server) mapperEmbedding(w http.ResponseWriter, r *http.Request) {

	numQnames := len(s.AnalysisService.MapperInfos["gtamap"].RecordsByQname)
	qnameEmbedding := make(map[string]int, numQnames)
	for qname, _ := range s.AnalysisService.MapperInfos["gtamap"].RecordsByQname {
		qnameEmbedding[qname] = len(qnameEmbedding)
	}

	fastaIndexPath := "/home/sam/Projects/gtamap-paper/pipeline/input/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"

	fastaIndex := fasta.ReadFastaIndexUsingPath(fastaIndexPath)

	numContigs := len(fastaIndex.Entries)
	contigEmbedding := make(map[string]int, numContigs)

	for contig, _ := range fastaIndex.Entries {
		contigEmbedding[contig] = len(contigEmbedding)
	}

	embeddings := make([]MapperEmbedding, 0, len(s.AnalysisService.MapperNames))

	for _, mapperName := range s.AnalysisService.MapperNames {

		mapperEmbedding := make([][]float32, 0)

		for qname, records := range s.AnalysisService.MapperInfos[mapperName].RecordsByQname {

			for _, record := range records {

				contigName := strings.Split(record.Rname, " ")[0]

				oneHotEncoding := make([]float32, 2)

				indexQname := qnameEmbedding[qname]
				indexContig := contigEmbedding[contigName]

				//oneHotEncoding[indexQname] = float32(1)
				//oneHotEncoding[numQnames+indexContig] = float32(1)
				oneHotEncoding[0] = float32(indexQname)
				oneHotEncoding[1] = float32(indexContig)

				posValue := float32(record.Pos) / float32(fastaIndex.Entries[contigName].Length)

				oneHotEncoding = append(oneHotEncoding, posValue)

				mapperEmbedding = append(mapperEmbedding, oneHotEncoding)
			}
		}

		embeddings = append(embeddings, MapperEmbedding{
			Mapper:    mapperName,
			Embedding: mapperEmbedding,
		})
	}
	//
	//distanceMatrix := make([][]float32, len(embeddings))
	//
	//for i := 0; i < len(embeddings); i++ {
	//	distanceMatrix[i] = make([]float32, len(embeddings))
	//}
	//
	//for i := 0; i < len(embeddings); i++ {
	//
	//	for j := i + 1; j < len(embeddings); j++ {
	//		embeddingA := embeddings[i]
	//		embeddingB := embeddings[j]
	//
	//		chamferDist := ChamferDistance(embeddingA.Embedding, embeddingB.Embedding)
	//		distanceMatrix[i][j] = chamferDist
	//		distanceMatrix[j][i] = chamferDist
	//	}
	//
	//	fmt.Printf("Done with embedding %d of %d\n", i+1, len(embeddings))
	//}
	//
	//mapperDistanceMatrix := MapperDistanceMatrix{
	//	Mappers:   make([]string, 0, len(embeddings)),
	//	Distances: distanceMatrix,
	//}
	//for _, embedding := range embeddings {
	//	mapperDistanceMatrix.Mappers = append(mapperDistanceMatrix.Mappers, embedding.Mapper)
	//}

	//// write the distance matrix to file
	//file, err := os.Create("distance_matrix.json")
	//if err != nil {
	//	http.Error(w, "Error creating distance matrix file", http.StatusInternalServerError)
	//	return
	//}
	//defer file.Close()
	//
	//err = json.NewEncoder(file).Encode(distanceMatix)
	//if err != nil {
	//	http.Error(w, "Error writing distance matrix to file", http.StatusInternalServerError)
	//	return
	//}

	//// read the distance matrix from file
	//file, err := os.Open("distance_matrix.json")
	//if err != nil {
	//	http.Error(w, "Error reading distance matrix file", http.StatusInternalServerError)
	//	return
	//}
	//defer file.Close()
	//
	//var mapperDistanceMatrix MapperDistanceMatrix
	//
	//err = json.NewDecoder(file).Decode(&mapperDistanceMatrix)
	//if err != nil {
	//	http.Error(w, "Error decoding distance matrix", http.StatusInternalServerError)
	//	return
	//}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(embeddings)
}

type MapperMultimmapingParallel struct {
	Qname   string         `json:"qname"`
	Mappers map[string]int `json:"mappers"`
}

func (s *Server) mapperMultimappingParallel(w http.ResponseWriter, r *http.Request) {

	qnameMap := make(map[string]MapperMultimmapingParallel)

	for qname, _ := range s.AnalysisService.MapperInfos["gtamap"].RecordsByQname {

		mapperMap := make(map[string]int, len(s.AnalysisService.MapperNames))
		for _, mapperName := range s.AnalysisService.MapperNames {
			mapperMap[mapperName] = 0
		}

		qnameMap[qname] = MapperMultimmapingParallel{
			Qname:   qname,
			Mappers: mapperMap,
		}
	}

	for mapperName, mapperInfo := range s.AnalysisService.MapperInfos {
		for qname, records := range mapperInfo.RecordsByQname {
			qnameMap[qname].Mappers[mapperName] = len(records)
		}
	}

	result := make([]MapperMultimmapingParallel, 0, len(qnameMap))

	for _, value := range qnameMap {
		result = append(result, value)
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(result)
}

// EuclideanDistance computes the Euclidean distance between two points.
func EuclideanDistance(a, b []float32) float32 {
	if len(a) != len(b) {
		panic("points must have the same dimension")
	}
	sum := float32(0.0)
	for i := range a {
		diff := a[i] - b[i]
		sum += diff * diff
	}
	return float32(math.Sqrt(float64(sum)))
}

// avgMinDist computes the average of minimum distances
// from each point in fromSet to the closest point in toSet.
func avgMinDist(fromSet, toSet [][]float32) float32 {
	total := float32(0.0)
	for _, p := range fromSet {
		minDist := float32(math.Inf(1))
		for _, q := range toSet {
			dist := EuclideanDistance(p, q)
			if dist < minDist {
				minDist = dist
			}
		}
		total += minDist
	}
	return total / float32(len(fromSet))
}

// ChamferDistance computes the Chamfer distance between two point sets.
func ChamferDistance(setA, setB [][]float32) float32 {
	return avgMinDist(setA, setB) + avgMinDist(setB, setA)
}

func serveGenomeFastaFile(w http.ResponseWriter, r *http.Request) {

	filePath := "/home/sam/Projects/gtamap-paper/pipeline/input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

	if _, err := os.Stat(filePath); os.IsNotExist(err) {
		http.Error(w, "File not found", http.StatusNotFound)
		return
	}

	http.ServeFile(w, r, filePath)
}

func serveGenomeFastaIndexFile(w http.ResponseWriter, r *http.Request) {

	filePath := "/home/sam/Projects/gtamap-paper/pipeline/input/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"

	if _, err := os.Stat(filePath); os.IsNotExist(err) {
		http.Error(w, "File not found", http.StatusNotFound)
		return
	}

	http.ServeFile(w, r, filePath)
}

func serveGenomeAnnotationFile(w http.ResponseWriter, r *http.Request) {

	filePath := "/home/sam/Projects/gtamap-paper/pipeline/input/Homo_sapiens.GRCh38.113.gtf"

	if _, err := os.Stat(filePath); os.IsNotExist(err) {
		http.Error(w, "File not found", http.StatusNotFound)
		return
	}

	http.ServeFile(w, r, filePath)
}

func serveBamFile(w http.ResponseWriter, r *http.Request) {

	filePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.gtamap.target.bam"

	if _, err := os.Stat(filePath); os.IsNotExist(err) {
		http.Error(w, "File not found", http.StatusNotFound)
		return
	}

	http.ServeFile(w, r, filePath)
}

func serveBamIndexFile(w http.ResponseWriter, r *http.Request) {

	filePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000173585/ENSG00000173585.gtamap.target.bam.bai"

	if _, err := os.Stat(filePath); os.IsNotExist(err) {
		http.Error(w, "File not found", http.StatusNotFound)
		return
	}

	http.ServeFile(w, r, filePath)
}

func getReadInfo(w http.ResponseWriter, r *http.Request) {

	response := analysis.ReadInfoTable()

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(response)
}

type IgvGenomeConfig struct {
	Id       string `json:"id"`
	Label    string `json:"label"`
	FastaUrl string `json:"fastaURL"`
	IndexUrl string `json:"indexURL"`
}

type IgvTrackConfig struct {
	Name        string `json:"name"`
	Format      string `json:"format"`
	DisplayMode string `json:"displayMode"`
	Url         string `json:"url"`
	IndexUrl    string `json:"indexURL"`
	Type        string `json:"type"`
}

func getTestGenomeConfig(w http.ResponseWriter, r *http.Request) {

	genomeConfig := IgvGenomeConfig{
		Id:       "BASE_GENOME",
		Label:    "BASE_GENOME_LABEL",
		FastaUrl: "http://localhost:8000/download/genome/fasta",
		IndexUrl: "http://localhost:8000/download/genome/fastaIndex",
	}

	track := IgvTrackConfig{
		Name:        "TEST_TRACK",
		Format:      "bam",
		DisplayMode: "EXPANDED",
		Url:         "http://localhost:8000/download/gtamap/bam",
		IndexUrl:    "http://localhost:8000/download/gtamap/bamIndex",
		Type:        "alignment",
	}

	tracks := make([]IgvTrackConfig, 0)
	tracks = append(tracks, track)

	response := map[string]interface{}{
		"genomeConfig": genomeConfig,
		"tracks":       tracks,
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(response)
}
