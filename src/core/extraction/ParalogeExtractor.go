package extraction

import (
	"bufio"
	"encoding/json"
	"fmt"
	"github.com/sirupsen/logrus"
	"io"
	"net/http"
	"os"
	"path/filepath"
)

type HomologyResponse struct {
	Data []struct {
		ID         string `json:"id"`
		Homologies []struct {
			Target struct {
				ID      string  `json:"id"`
				Species string  `json:"species"`
				PercID  float64 `json:"perc_id"`
			} `json:"target"`
			Type string `json:"type"`
		} `json:"homologies"`
	} `json:"data"`
}

type ParalogeSeq struct {
	ID      string
	Species string
	Perc    float64
}

func queryTargetId(geneid string, species string) map[string]struct{} {
	url := fmt.Sprintf("https://rest.ensembl.org/homology/id/%s/%s/?content-type=application/json;target_species=%s", species, geneid, species)

	resp, err := http.Get(url)
	if err != nil {
		logrus.Fatal("Failed to make request: %v", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		logrus.Fatal("API request failed with status: %s", resp.Status)
	}

	body, err := io.ReadAll(resp.Body)
	if err != nil {
		logrus.Fatal("Failed to read response body: %v", err)
	}

	var result HomologyResponse
	if err := json.Unmarshal(body, &result); err != nil {
		logrus.Fatal("Failed to parse JSON: %v", err)
	}

	paralogSeqs := make(map[string]struct{})
	for _, gene := range result.Data {
		for _, hom := range gene.Homologies {
			// we could filter for a certain homology percentage threshold if we wanted by using hom.Target.PercID
			// for now I create it in the format which is required by index..ExtractGeneSequenceFromGtfAndFastaForIndex
			paralogSeqs[hom.Target.ID] = struct{}{}
		}
	}
	// I am not using pointer here since the amount of paraloges will not be too high
	return paralogSeqs
}

func GetParaloges(targetGeneIds []string, species string) map[string]map[string]struct{} {
	paralogeSeqs := make(map[string]map[string]struct{})
	// receive paraloge genes per target gene
	for _, targetGene := range targetGeneIds {
		paralogeGenes := queryTargetId(targetGene, species)
		paralogeSeqs[targetGene] = paralogeGenes
	}
	return paralogeSeqs
}

func WriteParalogesPre(filename string, targetMap map[string][]string) error {
	file, err := os.Create(filename) // or os.OpenFile for appending
	if err != nil {
		return err
	}
	defer file.Close()

	writer := bufio.NewWriter(file)

	isFirstLine := true
	for target, paraloges := range targetMap {
		for _, paraloge := range paraloges {
			absPath, err := filepath.Abs(paraloge)
			if err != nil {
				fmt.Println("Error:", err)
			}
			if isFirstLine {
				_, err := writer.WriteString(target + "," + absPath)
				if err != nil {
					return err
				}
				isFirstLine = false
			} else {
				_, err := writer.WriteString("\n" + target + "," + absPath)
				if err != nil {
					return err
				}
			}
		}
	}
	return writer.Flush()
}
