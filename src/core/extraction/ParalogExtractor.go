package extraction

import (
	"bufio"
	"encoding/json"
	"fmt"
	"io"
	"io/fs"
	"net/http"
	"os"
	"path/filepath"
	"strings"

	"github.com/sirupsen/logrus"
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

type ParalogSeq struct {
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
	// I am not using pointer here since the amount of paralogs will not be too high
	return paralogSeqs
}

func GetParalogs(targetGeneIds []string, species string) map[string]map[string]struct{} {
	paralogSeqs := make(map[string]map[string]struct{})
	// receive paralog genes per target gene
	for _, targetGene := range targetGeneIds {
		paralogGenes := queryTargetId(targetGene, species)
		paralogSeqs[targetGene] = paralogGenes
	}
	return paralogSeqs
}

func WriteParalogsPre(filename string, targetMap map[string][]string) {
	file, err := os.Create(filename)
	if err != nil {
		logrus.Fatalf("Error creating: %s : %s", filename, err)
	}
	defer file.Close()

	writer := bufio.NewWriter(file)

	isFirstLine := true
	for target, paralogs := range targetMap {
		for _, paralog := range paralogs {
			absPath, err := filepath.Abs(paralog)
			if err != nil {
				logrus.Fatalf("Could not get abs path of %s: %s", paralog, err)
			}
			if isFirstLine {
				_, err := writer.WriteString(target + "," + absPath)
				if err != nil {
					logrus.Fatalf("Could not write to %s: %s", file.Name(), err)
				}
				isFirstLine = false
			} else {
				_, err := writer.WriteString("\n" + target + "," + absPath)
				if err != nil {
					logrus.Fatalf("Could not write to %s: %s", file.Name(), err)
				}
			}
		}
	}
	err = writer.Flush()
	if err != nil {
		logrus.Fatalf("Could not flush writer: %s", err)
	}
}

func GetAbsPathsPerTarget(targetParalogs map[string]map[string]struct{}, indexDirParalogPre *string) map[string][]string {
	indexPaths := make(map[string][]string)
	err := filepath.WalkDir(*indexDirParalogPre, func(path string, d fs.DirEntry, err error) error {
		if err != nil {
			return err
		}

		if !d.IsDir() {
			gtaiName := strings.Split(d.Name(), ".gtai")[0]
			for t, paralogs := range targetParalogs {
				_, exists := paralogs[gtaiName]
				if exists {
					paralogIndexPath := filepath.Join(*indexDirParalogPre, d.Name())
					indexPaths[t] = append(indexPaths[t], paralogIndexPath)
				}
			}
		}
		return nil
	})
	if err != nil {
		logrus.Fatalf("Error reading %s directory to get abs paths for paralog.csv meta file: %s", *indexDirParalogPre, err)
	}

	return indexPaths
}
