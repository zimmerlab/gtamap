package server

import (
	"encoding/json"
	"fmt"
	"github.com/gorilla/mux"
	"github.com/rs/cors"
	"log"
	"net/http"
	"time"
)

func enableCORS() *cors.Cors {
	return cors.New(cors.Options{
		AllowedOrigins: []string{"http://localhost:8080", "http://localhost:3000"},
		AllowedMethods: []string{"GET", "POST", "PUT", "DELETE", "OPTIONS"},
		AllowedHeaders: []string{"*"},
	})
}

func healthCheck(w http.ResponseWriter, r *http.Request) {
	response := map[string]string{
		"status":    "healthy",
		"timestamp": time.Now().Format(time.RFC3339),
	}
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(response)
}

func StartServer() {

	// Create router
	r := mux.NewRouter()

	// API routes
	api := r.PathPrefix("/api").Subrouter()
	api.HandleFunc("/health", healthCheck).Methods("GET")

	// Setup CORS
	c := enableCORS()
	handler := c.Handler(r)

	// Start server
	port := ":8000"
	fmt.Printf("Server starting on port %s\n", port)
	log.Fatal(http.ListenAndServe(port, handler))

}
