package server

import (
	"bytes"
	"encoding/json"
	"fmt"
	"log"
	"math"
	"net/http"
	"os"
	"os/exec"
	"slices"
	"strconv"
	"strings"
	"time"

	"github.com/KleinSamuel/gtamap/src/analysis"
	"github.com/KleinSamuel/gtamap/src/config"
	"github.com/KleinSamuel/gtamap/src/formats/fasta"
	"github.com/KleinSamuel/gtamap/src/formats/sam"
	"github.com/gorilla/mux"
	"github.com/rs/cors"
)

type Server struct {
	Router          *mux.Router
	Handler         *MappingDataHandler
	AnalysisService *analysis.AnalysisService // TODO: change all usages to Handler
}

func (s *Server) InitRoutes() {

	api := s.Router.PathPrefix("/api").Subrouter()

	api.HandleFunc("/health", s.healthCheck).Methods("GET")

	//api.HandleFunc("/genomeConfig", getTestGenomeConfig).Methods("GET")
	api.HandleFunc("/igvConfigTarget", getTargetRegionIgvConfig).Methods("GET")

	api.HandleFunc("/readSummaryTable", s.getReadSummaryTable).Methods("GET")

	api.HandleFunc("/upsetData", s.upsetData).Methods("GET")
	// api.HandleFunc("/upsetDataRecordPos", s.upsetDataRecordPos).Methods("GET")
	// api.HandleFunc("/upsetDataRecordPosCigar", s.upsetDataRecordPosCigar).Methods("GET")

	// api.HandleFunc("/getRecordXMapperByUpsetElementName", s.getRecordXMapperByUpsetElementName).Methods("GET")
	// api.HandleFunc("/getRecordXMapperByUpsetElementNames", s.getRecordXMapperByUpsetElementNames).Methods("POST")

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

	download.HandleFunc("/targetRegion/combinedSam", s.serveTargetRegionCombinedMappingSam).Methods("GET")
	download.HandleFunc("/targetRegion/combinedBam", s.serveTargetRegionCombinedMappingBam).Methods("GET")
	download.HandleFunc("/targetRegion/combinedBamIndex", s.serveTargetRegionCombinedMappingBamIndex).Methods("GET")

	readDetailsRouter := api.PathPrefix("/readDetails").Subrouter()

	readDetailsRouter.HandleFunc("/readInfo", s.getReadDetailsInfo).Methods("GET")
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

func (s *Server) InitMappingDataHandler(fastaFilePath string, fastaIndexFilePath string) {
	handler, err := NewMappingDataHandler(fastaFilePath, fastaIndexFilePath)
	if err != nil {
		panic("could not create mapping data handler")
	}
	s.Handler = handler
}

type UpsetElement struct {
	Name string   `json:"name"`
	Sets []string `json:"sets"`
}

func (s *Server) upsetData(w http.ResponseWriter, r *http.Request) {

	onlyTargetRegion := false
	if val := r.URL.Query().Get("onlyTargetRegion"); val == "true" {
		onlyTargetRegion = true
	}
	usePosition := false
	if val := r.URL.Query().Get("usePosition"); val == "true" {
		usePosition = true
	}
	useCigar := false
	if val := r.URL.Query().Get("useCigar"); val == "true" {
		useCigar = true
	}

	data := make([]UpsetElement, 0)

	for _, qnameClust := range s.Handler.QnameCluster {

		set := make(map[string][]string)

		for _, c1 := range qnameClust.ClusterR1.Records {

			// skip reads that do not map to the target region of flag is set
			if onlyTargetRegion && (c1.Rname != config.GetTargetContig() ||
				!c1.MappedGenome.Overlaps(config.GetTargetStart(), config.GetTargetEnd())) {
				continue
			}

			readKey := getReadKey(c1, "R1", usePosition, useCigar)

			if _, exists := set[readKey]; !exists {
				set[readKey] = make([]string, 0)
			}
			set[readKey] = append(set[readKey], s.Handler.MapperInfos[c1.MapperIndex].MapperName)
		}
		for _, c2 := range qnameClust.ClusterR2.Records {
			// skip reads that do not map to the target region of flag is set
			if onlyTargetRegion && (c2.Rname != config.GetTargetContig() ||
				!c2.MappedGenome.Overlaps(config.GetTargetStart(), config.GetTargetEnd())) {
				continue
			}

			readKey := getReadKey(c2, "R2", usePosition, useCigar)

			if _, exists := set[readKey]; !exists {
				set[readKey] = make([]string, 0)
			}
			set[readKey] = append(set[readKey], s.Handler.MapperInfos[c2.MapperIndex].MapperName)
		}

		for qname, mapperNames := range set {
			data = append(data, UpsetElement{
				Name: qname,
				Sets: mapperNames,
			})
		}
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(data)
}

func getReadKey(record *EnhancedRecord, prefix string, usePosition bool, useCigar bool) string {
	readKey := prefix + "::" + record.Qname
	if usePosition {
		readKey += "::" + strconv.Itoa(record.Pos)
	}
	if useCigar {
		readKey += "::" + record.UniformCigar()
	}
	return readKey
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

// func (s *Server) getRecordXMapperByUpsetElementName(w http.ResponseWriter, r *http.Request) {
//
// 	name := r.URL.Query().Get("name")
//
// 	fmt.Println("upset element name:", name)
//
// 	parts := strings.Split(name, "::")
//
// 	if len(parts) < 1 || len(parts) > 3 {
// 		http.Error(w, "Invalid name format", http.StatusBadRequest)
// 		return
// 	}
//
// 	qname := parts[0]
// 	pos, err := strconv.Atoi(parts[1])
// 	if err != nil {
// 		http.Error(w, "Invalid position in name", http.StatusBadRequest)
// 		return
// 	}
// 	cigar := ""
// 	if len(parts) == 3 {
// 		cigar = parts[2]
// 	}
//
// 	mapperToRecord := make(map[string][]sam.Record)
//
// 	for mapperName, mapperInfo := range s.AnalysisService.MapperInfos {
//
// 		for _, record := range mapperInfo.ParsedFile.Records {
//
// 			qnameMatches := record.Qname == qname
// 			posMatches := record.Pos == pos
// 			cigarMatches := true
// 			if len(parts) == 3 {
// 				cigarMatches = record.UniformCigar() == cigar
// 			}
//
// 			if qnameMatches && posMatches && cigarMatches {
// 				if _, exists := mapperToRecord[mapperName]; !exists {
// 					mapperToRecord[mapperName] = make([]sam.Record, 0)
// 				}
// 				mapperToRecord[mapperName] = append(mapperToRecord[mapperName], record)
// 			}
// 		}
// 	}
//
// 	w.Header().Set("Content-Type", "application/json")
// 	json.NewEncoder(w).Encode(mapperToRecord)
// }

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
	if _, err := os.Stat(config.GetTargetFasta()); os.IsNotExist(err) {
		http.Error(w, "File not found", http.StatusNotFound)
		return
	}
	http.ServeFile(w, r, config.GetTargetFasta())
}

func serveTargetRegionFastaIndexFile(w http.ResponseWriter, r *http.Request) {
	if _, err := os.Stat(config.GetTargetFastaIndex()); os.IsNotExist(err) {
		http.Error(w, "File not found", http.StatusNotFound)
		return
	}
	http.ServeFile(w, r, config.GetTargetFastaIndex())
}

func (s *Server) getCombinedSam() string {
	var builder strings.Builder

	builder.WriteString("@HD\tVN:1.6\tSO:unsorted\n")
	builder.WriteString(fmt.Sprintf("@SQ\tSN:%s\tLN:%d\n", config.GetTargetContig(), config.GetTargetEnd()-config.GetTargetStart()+1))

	for mIndex, mapperInfo := range s.Handler.MapperInfos {
		for _, record := range s.Handler.RecordsByMapper[mIndex] {
			// only include records that map to the target region
			isOnTargetContig := record.Rname == config.GetTargetContig()
			isInTargetRegion := record.Pos >= config.GetTargetStart() && record.Pos <= config.GetTargetEnd()
			if !isOnTargetContig || !isInTargetRegion {
				continue
			}

			builder.WriteString(record.StringWithMapperInfo(config.GetTargetStart(), mapperInfo.MapperName))
			builder.WriteString("\n")
		}
	}

	return builder.String()
}

func (s *Server) serveTargetRegionCombinedMappingSam(w http.ResponseWriter, r *http.Request) {

	samString := s.getCombinedSam()

	w.Header().Set("Content-Type", "plain/text")
	w.Write([]byte(samString))
}

func (s *Server) serveTargetRegionCombinedMappingBam(w http.ResponseWriter, r *http.Request) {

	samString := s.getCombinedSam()

	bamBytes, err := samToSortedBam(samString)
	if err != nil {
		http.Error(w, fmt.Sprintf("Error converting SAM to BAM: %v", err), http.StatusInternalServerError)
		return
	}

	w.Header().Set("Content-Type", "application/octet-stream")
	w.Write(bamBytes)
}

func (s *Server) serveTargetRegionCombinedMappingBamIndex(w http.ResponseWriter, r *http.Request) {

	samString := s.getCombinedSam()

	bamBytes, err := samToSortedBam(samString)
	if err != nil {
		http.Error(w, fmt.Sprintf("Error converting SAM to BAM: %v", err), http.StatusInternalServerError)
		return
	}

	bamIndex, err := bamToBamIndex(bamBytes)
	if err != nil {
		http.Error(w, fmt.Sprintf("Error creating BAM index: %v", err), http.StatusInternalServerError)
		return
	}

	w.Header().Set("Content-Type", "plain/text")
	w.Write([]byte(bamIndex))
}

func samToSortedBam(samString string) ([]byte, error) {

	cmd := exec.Command("samtools", "sort", "-o", "-", "-")

	var stderrBuf bytes.Buffer
	cmd.Stderr = &stderrBuf

	stdin, err := cmd.StdinPipe()
	if err != nil {
		return nil, fmt.Errorf("error creating stdin pipe: %w, stderr: %s", err, stderrBuf.String())
	}
	var bamBuf bytes.Buffer
	cmd.Stdout = &bamBuf

	if err := cmd.Start(); err != nil {
		return nil, fmt.Errorf("error starting samtools: %w, stderr: %s", err, stderrBuf.String())
	}

	_, err = stdin.Write([]byte(samString))
	if err != nil {
		stdin.Close()
		return nil, fmt.Errorf("error writing to stdin: %w, stderr: %s", err, stderrBuf.String())
	}
	stdin.Close()

	if err := cmd.Wait(); err != nil {
		return nil, fmt.Errorf("samtools failed: %w, stderr: %s", err, stderrBuf.String())
	}

	if bamBuf.Len() == 0 {
		return nil, fmt.Errorf("no output from samtools")
	}

	return bamBuf.Bytes(), nil
}

func bamToBamIndex(bamBytes []byte) (string, error) {

	cmd := exec.Command("samtools", "index", "-o", "-", "-")

	var stderrBuf bytes.Buffer
	cmd.Stderr = &stderrBuf

	stdin, err := cmd.StdinPipe()
	if err != nil {
		return "", fmt.Errorf("error creating stdin pipe: %w, stderr: %s", err, stderrBuf.String())
	}

	var indexBuf bytes.Buffer
	cmd.Stdout = &indexBuf

	if err := cmd.Start(); err != nil {
		return "", fmt.Errorf("error starting samtools: %w, stderr: %s", err, stderrBuf.String())
	}

	_, err = stdin.Write(bamBytes)
	if err != nil {
		stdin.Close()
		return "", fmt.Errorf("error writing to stdin: %w, stderr: %s", err, stderrBuf.String())
	}
	stdin.Close()

	if err := cmd.Wait(); err != nil {
		return "", fmt.Errorf("samtools failed: %w, stderr: %s", err, stderrBuf.String())
	}

	if indexBuf.Len() == 0 {
		return "", fmt.Errorf("no output from samtools")
	}

	return indexBuf.String(), nil
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
	Height      int    `json:"height,omitempty"`
	MaxHeight   int    `json:"maxHeight,omitempty"`
}

func getTargetRegionIgvConfig(w http.ResponseWriter, r *http.Request) {

	genomeConfig := IgvGenomeConfig{
		Id:       "Target Region (chr3:123456-789012)",
		Label:    "BASE_GENOME_LABEL",
		FastaUrl: "http://localhost:8000/download/targetRegion/fasta",
		IndexUrl: "http://localhost:8000/download/targetRegion/fastaIndex",
	}

	track := IgvTrackConfig{
		Name:        "TEST_TRACK",
		Format:      "sam",
		DisplayMode: "EXPANDED",
		Url:         "http://localhost:8000/download/targetRegion/combinedBam",
		IndexUrl:    "http://localhost:8000/download/targetRegion/combinedBamIndex",
		Type:        "alignment",
		Height:      800,
		MaxHeight:   1000,
	}

	tracks := make([]IgvTrackConfig, 0)
	tracks = append(tracks, track)

	// convert global location to local location
	targetRegionStrLocal := fmt.Sprintf("%s:%d-%d", config.GetTargetContig(), 1, config.GetTargetEnd()-config.GetTargetStart()+1)

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

type ReadOverviewInfo struct {
	Qname           string              `json:"qname"`
	ReadLengthR1    int                 `json:"readLengthR1"`
	ReadLengthR2    int                 `json:"readLengthR2"`
	NumMappedBy     int                 `json:"numMappedBy"`
	MappedBy        []string            `json:"mappedBy"`
	NumLocations    int                 `json:"numLocations"`
	Locations       []*ReadLocationInfo `json:"locations"`
	DistanceScore   float64             `json:"distanceScore,omitempty"`   // optional, for distance score
	ConfidenceLevel int                 `json:"confidenceLevel,omitempty"` // optional, for confidence level
}

type ReadLocationInfo struct {
	Pair          string   `json:"pairType"` // "first" or "second"
	Contig        string   `json:"contigName"`
	Strand        bool     `json:"isForwardStrand"` // true for forward strand, false for reverse strand
	Position      int      `json:"position"`
	Cigar         string   `json:"cigar"`
	CigarDetailed string   `json:"cigarDetailed,omitempty"` // optional, for detailed CIGAR representation
	NumMismatches int      `json:"numMismatches"`
	NumGaps       int      `json:"numGaps"`
	NumMappedBy   int      `json:"numMappedBy"`
	MappedBy      []string `json:"mappedBy"`
	ReadIndices   []int    `json:"readIndices"`
}

func (s *Server) ReadSummaryTableData() []ReadOverviewInfo {

	reads := make([]ReadOverviewInfo, 0, len(s.Handler.ReadInfos))

	for _, readInfo := range s.Handler.ReadInfos {

		qnameCluster := s.Handler.QnameCluster[readInfo.Qname]

		numClusterR1 := len(qnameCluster.ClusterR1.SimilarRecords)
		numClusterR2 := len(qnameCluster.ClusterR2.SimilarRecords)

		locations := make([]*ReadLocationInfo, 0, numClusterR1+numClusterR2)

		mappedBy := make(map[string]bool)

		for _, cR1 := range qnameCluster.ClusterR1.SimilarRecords {

			mappedByCluster := make(map[string]bool)
			for _, r := range cR1 {
				mappedBy[s.Handler.MapperInfos[r.MapperIndex].MapperName] = true
				mappedByCluster[s.Handler.MapperInfos[r.MapperIndex].MapperName] = true
			}
			mappedByClusterList := make([]string, 0, len(mappedByCluster))
			for mapperName := range mappedByCluster {
				mappedByClusterList = append(mappedByClusterList, mapperName)
			}
			slices.Sort(mappedByClusterList)

			locationInfo := ReadLocationInfo{
				Pair:          "first",
				Contig:        cR1[0].Rname,
				Strand:        !cR1[0].Flag.IsReverseStrand(),
				Position:      cR1[0].Pos,
				Cigar:         cR1[0].CigarObj.StringUniform(),
				CigarDetailed: cR1[0].CigarObj.String(),
				NumMismatches: cR1[0].CigarObj.GetNumMismatches(),
				NumGaps:       cR1[0].CigarObj.GetNumGaps(),
				NumMappedBy:   len(mappedByCluster),
				MappedBy:      mappedByClusterList,
				ReadIndices:   make([]int, 0),
			}

			for _, r := range cR1 {
				locationInfo.ReadIndices = append(locationInfo.ReadIndices, r.IndexInSam)
			}

			locations = append(locations, &locationInfo)
		}

		for _, cR2 := range qnameCluster.ClusterR2.SimilarRecords {

			mappedByCluster := make(map[string]bool)
			for _, r := range cR2 {
				mappedBy[s.Handler.MapperInfos[r.MapperIndex].MapperName] = true
				mappedByCluster[s.Handler.MapperInfos[r.MapperIndex].MapperName] = true
			}
			mappedByClusterList := make([]string, 0, len(mappedByCluster))
			for mapperName := range mappedByCluster {
				mappedByClusterList = append(mappedByClusterList, mapperName)
			}
			slices.Sort(mappedByClusterList)

			locationInfo := ReadLocationInfo{
				Pair:          "second",
				Contig:        cR2[0].Rname,
				Strand:        !cR2[0].Flag.IsReverseStrand(),
				Position:      cR2[0].Pos,
				Cigar:         cR2[0].CigarObj.StringUniform(),
				CigarDetailed: cR2[0].CigarObj.String(),
				NumMismatches: cR2[0].CigarObj.GetNumMismatches(),
				NumGaps:       cR2[0].CigarObj.GetNumGaps(),
				NumMappedBy:   len(mappedByCluster),
				MappedBy:      mappedByClusterList,
				ReadIndices:   make([]int, 0),
			}

			for _, r := range cR2 {
				locationInfo.ReadIndices = append(locationInfo.ReadIndices, r.IndexInSam)
			}

			locations = append(locations, &locationInfo)
		}

		mappedByList := make([]string, 0, len(mappedBy))
		for mapperName := range mappedBy {
			mappedByList = append(mappedByList, mapperName)
		}
		slices.Sort(mappedByList)

		distanceMean := (qnameCluster.ClusterR1.RecordDistanceMean + qnameCluster.ClusterR2.RecordDistanceMean) / 2.0

		readOverviewInfo := ReadOverviewInfo{
			Qname:           readInfo.Qname,
			ReadLengthR1:    len(readInfo.SequenceFirstOfPair),
			ReadLengthR2:    len(readInfo.SequenceSecondOfPair),
			NumMappedBy:     len(mappedBy),
			MappedBy:        mappedByList,
			NumLocations:    len(locations),
			Locations:       locations,
			DistanceScore:   distanceMean,
			ConfidenceLevel: qnameCluster.ConfidenceLevel,
		}

		reads = append(reads, readOverviewInfo)
	}

	return reads
}

type RecordWithMapperInfo struct {
	Record *sam.Record `json:"record"`
	Mapper string      `json:"mapper"`
}

type Interval struct {
	Contig string `json:"contig"`
	Start  int    `json:"start"`
	End    int    `json:"end"`
}

type Cluster struct {
	Interval      Interval `json:"interval"`
	RecordIndices []int    `json:"recordIndices"`
}

type ReadDetailsInfo struct {
	Records []RecordWithMapperInfo `json:"records"`
	Cluster []Cluster              `json:"cluster"`
}

func (s *Server) getReadDetailsInfo(w http.ResponseWriter, r *http.Request) {

	qname := r.URL.Query().Get("qname")

	if qname == "" {
		http.Error(w, "Missing qname parameter", http.StatusBadRequest)
		return
	}

	readDetails := ReadDetailsInfo{
		Records: make([]RecordWithMapperInfo, 0),
		Cluster: make([]Cluster, 0),
	}

	for _, mapperName := range s.AnalysisService.MapperNames {

		for _, record := range s.AnalysisService.MapperInfos[mapperName].RecordsByQname[qname] {

			recordWithMapperInfo := RecordWithMapperInfo{
				Record: record,
				Mapper: mapperName,
			}

			readDetails.Records = append(readDetails.Records, recordWithMapperInfo)

			foundCluster := false
			for _, cluster := range readDetails.Cluster {

				if cluster.Interval.Contig == record.Rname && cluster.Interval.Start <= record.Pos && cluster.Interval.End >= record.Pos {
					foundCluster = true

					readDetails.Records = append(readDetails.Records, recordWithMapperInfo)
					// add the record index to the cluster
					cluster.RecordIndices = append(cluster.RecordIndices, len(readDetails.Records)-1)
				}
			}

			if !foundCluster {
				// create a new cluster for this record
				newCluster := Cluster{
					Interval: Interval{
						Contig: record.Rname,
						Start:  record.Pos,
						End:    record.Pos + len(record.Seq) - 1,
					},
					RecordIndices: []int{len(readDetails.Records) - 1},
				}
				readDetails.Cluster = append(readDetails.Cluster, newCluster)
			}

		}

	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(readDetails)
}
