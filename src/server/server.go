package server

import (
	"encoding/json"
	"fmt"
	"github.com/KleinSamuel/gtamap/src/analysis"
	"github.com/gorilla/mux"
	"github.com/rs/cors"
	"log"
	"net/http"
	"os"
	"strconv"
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
