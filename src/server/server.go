package server

import (
	"encoding/json"
	"fmt"
	"github.com/KleinSamuel/gtamap/src/analysis"
	"github.com/KleinSamuel/gtamap/src/formats/fasta"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
	"github.com/gorilla/mux"
	"github.com/rs/cors"
	"github.com/sirupsen/logrus"
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

	//api.HandleFunc("/genomeConfig", getTestGenomeConfig).Methods("GET")
	api.HandleFunc("/genomeConfig", getTargetRegionIgvConfig).Methods("GET")

	api.HandleFunc("/readSummaryTable", s.getReadSummaryTable).Methods("GET")

	api.HandleFunc("/upsetDataRead", s.upsetDataRead).Methods("GET")
	api.HandleFunc("/upsetDataRecordPos", s.upsetDataRecordPos).Methods("GET")
	api.HandleFunc("/upsetDataRecordPosCigar", s.upsetDataRecordPosCigar).Methods("GET")

	api.HandleFunc("/getRecordXMapperByUpsetElementName", s.getRecordXMapperByUpsetElementName).Methods("GET")
	api.HandleFunc("/getRecordXMapperByUpsetElementNames", s.getRecordXMapperByUpsetElementNames).Methods("POST")

	//api.HandleFunc("/readNumMatchesDistribution", s.readNumMatchesDistribution).Methods("GET")
	api.HandleFunc("/readMultimappingDensityData", s.readMultimappingDensityData).Methods("GET")
	api.HandleFunc("/readMultimappingDensityDataHeatmap", s.readMultimappingDensityDataHeatmap).Methods("GET")
	api.HandleFunc("/mapperEmbedding", s.mapperEmbedding).Methods("GET")
	api.HandleFunc("/mapperMultimappingParallel", s.mapperMultimappingParallel).Methods("GET")
	api.HandleFunc("/mapperDistances", s.mapperDistances).Methods("GET")

	download := s.Router.PathPrefix("/download").Subrouter()

	download.HandleFunc("/genome/fasta", serveGenomeFastaFile).Methods("GET")
	download.HandleFunc("/genome/fastaIndex", serveGenomeFastaIndexFile).Methods("GET")

	download.HandleFunc("/genome/gtf", serveGenomeAnnotationFile).Methods("GET")
	download.HandleFunc("/gtamap/bam", serveBamFile).Methods("GET")
	download.HandleFunc("/gtamap/bamIndex", serveBamIndexFile).Methods("GET")

	download.HandleFunc("/targetRegion/fasta", serveTargetRegionFastaFile).Methods("GET")
	download.HandleFunc("/targetRegion/fastaIndex", serveTargetRegionFastaIndexFile).Methods("GET")
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
				recordPosMap[id] = make(map[string]bool)
			}

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

func (s *Server) upsetDataRecordPosCigar(w http.ResponseWriter, r *http.Request) {

	data := make([]UpsetElement, 0)

	recordPosMap := make(map[string]map[string]bool)

	for mapperName, mapperInfo := range s.AnalysisService.MapperInfos {

		for _, record := range mapperInfo.ParsedFile.Records {

			id := record.Qname + "::" + strconv.Itoa(record.Pos) + "::" + record.UniformCigar()

			if _, exists := recordPosMap[id]; !exists {
				recordPosMap[id] = make(map[string]bool)
			}

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

type UpsetElementNamesDto struct {
	Names []string `json:"names"`
}

func (s *Server) getRecordXMapperByUpsetElementNames(w http.ResponseWriter, r *http.Request) {

	var dto UpsetElementNamesDto
	err := json.NewDecoder(r.Body).Decode(&dto)
	if err != nil {
		http.Error(w, "Invalid request body", http.StatusBadRequest)
		return
	}

	fmt.Println("len names:", len(dto.Names))
}

func testRecords(names []string) {

}

func (s *Server) getRecordXMapperByUpsetElementName(w http.ResponseWriter, r *http.Request) {

	name := r.URL.Query().Get("name")

	fmt.Println("upset element name:", name)

	parts := strings.Split(name, "::")

	if len(parts) < 1 || len(parts) > 3 {
		http.Error(w, "Invalid name format", http.StatusBadRequest)
		return
	}

	qname := parts[0]
	pos, err := strconv.Atoi(parts[1])
	if err != nil {
		http.Error(w, "Invalid position in name", http.StatusBadRequest)
		return
	}
	cigar := ""
	if len(parts) == 3 {
		cigar = parts[2]
	}

	mapperToRecord := make(map[string][]sam.Record)

	for mapperName, mapperInfo := range s.AnalysisService.MapperInfos {

		for _, record := range mapperInfo.ParsedFile.Records {

			qnameMatches := record.Qname == qname
			posMatches := record.Pos == pos
			cigarMatches := true
			if len(parts) == 3 {
				cigarMatches = record.UniformCigar() == cigar
			}

			if qnameMatches && posMatches && cigarMatches {
				if _, exists := mapperToRecord[mapperName]; !exists {
					mapperToRecord[mapperName] = make([]sam.Record, 0)
				}
				mapperToRecord[mapperName] = append(mapperToRecord[mapperName], record)
			}
		}
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(mapperToRecord)
}

//func (s *Server) readNumMatchesDistribution(w http.ResponseWriter, r *http.Request) {
//
//	// TODO: currently only gtamap
//
//	data := s.AnalysisService.MapperInfos["gtamap"].RecordsByQname
//
//	recordInfoList := make([]RecordInfo, 0, len(data))
//
//	for qname, records := range data {
//		recordIds := make([]int, 0, len(records))
//		for _, record := range records {
//			recordIds = append(recordIds, record.Id)
//		}
//
//		recordInfo := RecordInfo{
//			Qname:     qname,
//			RecordIds: recordIds,
//		}
//
//		recordInfoList = append(recordInfoList, recordInfo)
//	}
//
//	w.Header().Set("Content-Type", "application/json")
//	json.NewEncoder(w).Encode(recordInfoList)
//}

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

type MapperDistanceMatrix64 struct {
	Mappers   []string    `json:"mappers"`
	Distances [][]float64 `json:"distances"`
}

func (s *Server) mapperDistances(w http.ResponseWriter, r *http.Request) {

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

	distanceMatrix := make([][]float64, len(s.AnalysisService.MapperNames))
	for i := 0; i < len(s.AnalysisService.MapperNames); i++ {
		distanceMatrix[i] = make([]float64, len(s.AnalysisService.MapperNames))
	}

	for i := 0; i < len(s.AnalysisService.MapperNames); i++ {

		mapperNameA := s.AnalysisService.MapperNames[i]
		recordsA := s.AnalysisService.MapperInfos[mapperNameA].RecordsByQname

		for j := i + 1; j < len(s.AnalysisService.MapperNames); j++ {

			numSame := 0
			numOther := 0

			mapperNameB := s.AnalysisService.MapperNames[j]

			for qnameB, recordsB := range s.AnalysisService.MapperInfos[mapperNameB].RecordsByQname {

				for _, recordB := range recordsB {

					found := false

					for _, recordA := range recordsA[qnameB] {

						if recordA.Rname == recordB.Rname && recordA.Pos == recordB.Pos {
							found = true
							break
						}
					}

					if found {
						numSame++
					} else {
						numOther++
					}
				}
			}

			distance := float64(numSame)
			distanceMatrix[i][j] = distance
			distanceMatrix[j][i] = distance
		}
	}

	result := MapperDistanceMatrix64{
		Mappers:   s.AnalysisService.MapperNames,
		Distances: distanceMatrix,
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(result)
}

func serveTargetRegionFastaFile(w http.ResponseWriter, r *http.Request) {
	if _, err := os.Stat(GetTargetFasta()); os.IsNotExist(err) {
		http.Error(w, "File not found", http.StatusNotFound)
		return
	}
	http.ServeFile(w, r, GetTargetFasta())
}

func serveTargetRegionFastaIndexFile(w http.ResponseWriter, r *http.Request) {
	if _, err := os.Stat(GetTargetFastaIndex()); os.IsNotExist(err) {
		http.Error(w, "File not found", http.StatusNotFound)
		return
	}
	http.ServeFile(w, r, GetTargetFastaIndex())
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

	filePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000130305/ENSG00000130305.gtamap.target.bam"

	if _, err := os.Stat(filePath); os.IsNotExist(err) {
		http.Error(w, "File not found", http.StatusNotFound)
		return
	}

	http.ServeFile(w, r, filePath)
}

func serveBamIndexFile(w http.ResponseWriter, r *http.Request) {

	filePath := "/home/sam/Projects/gtamap-paper/pipeline/output/ENSG00000130305/ENSG00000130305.gtamap.target.bam.bai"

	if _, err := os.Stat(filePath); os.IsNotExist(err) {
		http.Error(w, "File not found", http.StatusNotFound)
		return
	}

	http.ServeFile(w, r, filePath)
}

func (s *Server) getReadSummaryTable(w http.ResponseWriter, r *http.Request) {

	response := s.ReadSummaryTableData()

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

func getTargetRegionIgvConfig(w http.ResponseWriter, r *http.Request) {

	genomeConfig := IgvGenomeConfig{
		Id:       "Target Region (chr3:123456-789012)",
		Label:    "BASE_GENOME_LABEL",
		FastaUrl: "http://localhost:8000/download/targetRegion/fasta",
		IndexUrl: "http://localhost:8000/download/targetRegion/fastaIndex",
	}

	//track := IgvTrackConfig{
	//	Name:        "TEST_TRACK",
	//	Format:      "bam",
	//	DisplayMode: "EXPANDED",
	//	Url:         "http://localhost:8000/download/gtamap/bam",
	//	IndexUrl:    "http://localhost:8000/download/gtamap/bamIndex",
	//	Type:        "alignment",
	//}

	tracks := make([]IgvTrackConfig, 0)
	//tracks = append(tracks, track)

	// convert global location to local location
	targetRegionStr := GetTargetRegion()
	startEnd := strings.Split(strings.Split(targetRegionStr, ":")[1], "-")
	start, err := strconv.Atoi(startEnd[0])
	if err != nil {
		http.Error(w, "Invalid start position in target region", http.StatusBadRequest)
	}
	end, err := strconv.Atoi(startEnd[1])
	if err != nil {
		http.Error(w, "Invalid end position in target region", http.StatusBadRequest)
	}
	targetRegionStrLocal := fmt.Sprintf("%s:%d-%d", GetTargetName(), 1, end-start)

	response := map[string]interface{}{
		"genomeConfig": genomeConfig,
		"tracks":       tracks,
		"location":     targetRegionStrLocal,
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(response)

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

type MapperResult struct {
	Name   string `json:"name"`
	Path   string `json:"path"`
	Target string `json:"target"`
}

// GetMapperResults retrieves the results of mappers from the output directory.
// It assumes that the output directory contains files named in the format "mapperName.target.sam".
func (s *Server) GetMapperResults() []MapperResult {

	results := make([]MapperResult, 0)

	files, err := os.ReadDir(GetRunDir())
	if err != nil {
		logrus.Error("Error reading output directory: ", err)
		return results
	}

	for _, file := range files {
		if strings.HasSuffix(file.Name(), ".sam") {
			parts := strings.Split(file.Name(), ".")

			if len(parts) < 3 {
				logrus.Warn("Skipping mapper file: ", file.Name())
				continue
			}

			result := MapperResult{
				Name:   parts[1],
				Path:   GetRunDir() + "/" + file.Name(),
				Target: parts[2],
			}

			results = append(results, result)
		}
	}

	return results
}

func (s *Server) LoadMapperResults() {

	results := s.GetMapperResults()

	for _, result := range results {

		if err := s.AnalysisService.AddMapperInfo(result.Name, result.Path, result.Target); err != nil {
			log.Fatalf("Error initializing analysis service: %v", err)
		}

		logrus.WithFields(logrus.Fields{
			"mapperName":        result.Name,
			"mapperSamFilePath": result.Path,
			"target":            result.Target,
		}).Debug("added mapper result")
	}
}

type ReadOverviewInfo struct {
	Qname        string              `json:"qname"`
	Length       int                 `json:"readLength"`
	NumMappedBy  int                 `json:"numMappedBy"`
	MappedBy     []string            `json:"mappedBy"`
	NumLocations int                 `json:"numLocations"`
	Locations    []*ReadLocationInfo `json:"locations"`
}

type ReadLocationInfo struct {
	Pair          string   `json:"pairType"` // "first" or "second"
	Contig        string   `json:"contigName"`
	Strand        bool     `json:"isForwardStrand"` // true for forward strand, false for reverse strand
	Position      int      `json:"position"`
	CigarString   string   `json:"cigarString"`
	NumMismatches int      `json:"numMismatches"`
	NumGaps       int      `json:"numGaps"`
	NumMappedBy   int      `json:"numMappedBy"`
	MappedBy      []string `json:"mappedBy"`
	ReadIndices   []int    `json:"readIndices"`
}

func (s *Server) ReadSummaryTableData() []ReadOverviewInfo {

	readSummary := make(map[string]*ReadOverviewInfo)

	for mapperName, info := range s.AnalysisService.MapperInfos {

		for qname, records := range info.RecordsByQname {

			if _, exists := readSummary[qname]; !exists {
				readSummary[qname] = &ReadOverviewInfo{
					Qname:        qname,
					Length:       len(records[0].Seq),
					NumMappedBy:  0,
					MappedBy:     make([]string, 0),
					NumLocations: 0,
					Locations:    make([]*ReadLocationInfo, 0),
				}
			}

			readInfo := readSummary[qname]
			readInfo.MappedBy = append(readInfo.MappedBy, mapperName)

			for _, record := range records {

				found := false
				for _, location := range readInfo.Locations {
					if location.Strand == !record.Flag.IsReverseStrand() &&
						location.Contig == record.Rname &&
						location.Position == record.Pos &&
						location.CigarString == record.UniformCigar() {

						// do not add the same mapper multiple times (could be because the same location is used
						// for multiple read pairs, especially in the current 0.2 gtamap version)
						for _, mappedBy := range location.MappedBy {
							if mappedBy == mapperName {
								found = true
								break
							}
						}

						if found {
							break
						}

						location.MappedBy = append(location.MappedBy, mapperName)
						location.ReadIndices = append(location.ReadIndices, record.IndexInSam)
						location.NumMappedBy++
						found = true
						break
					}
				}

				if found {
					continue
				}

				pair := "none"
				if record.Flag.IsFirstInPair() {
					pair = "first"
				} else if record.Flag.IsSecondInPair() {
					pair = "second"
				}

				locationInfo := ReadLocationInfo{
					Pair:          pair,
					Contig:        record.Rname,
					Strand:        !record.Flag.IsReverseStrand(),
					Position:      record.Pos,
					CigarString:   record.UniformCigar(),
					NumMismatches: -1,
					NumGaps:       -1,
					NumMappedBy:   1,
					MappedBy:      []string{mapperName},
					ReadIndices:   make([]int, 0),
				}

				readInfo.Locations = append(readInfo.Locations, &locationInfo)
				locationInfo.ReadIndices = append(locationInfo.ReadIndices, record.IndexInSam)
			}
		}
	}

	reads := make([]ReadOverviewInfo, 0)

	for _, readInfo := range readSummary {
		readInfo.NumMappedBy = len(readInfo.MappedBy)
		readInfo.NumLocations = len(readInfo.Locations)
		reads = append(reads, *readInfo)
	}

	return reads
}
