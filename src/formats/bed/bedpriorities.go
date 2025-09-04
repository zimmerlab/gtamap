package bed

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

func ReadPrioritiesFromPath(filePath string) (map[string]int, error) {

	file, err := os.Open(filePath)
	if err != nil {
		return nil, fmt.Errorf("failed to open file %s: %w", filePath, err)
	}
	defer file.Close()

	return ReadPriorities(file)
}

func ReadPriorities(file *os.File) (map[string]int, error) {

	pMap := make(map[string]int)
	pSet := make(map[int]struct{})

	scanner := bufio.NewScanner(file)
	lineNum := 0

	for scanner.Scan() {
		lineNum++
		line := strings.TrimSpace(scanner.Text())

		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}

		parts := strings.Fields(line)
		if len(parts) != 2 {
			return nil, fmt.Errorf("invalid format at line %d: expected 2 fields, got %d", lineNum, len(parts))
		}

		priority, err := strconv.Atoi(parts[1])
		if err != nil {
			return nil, fmt.Errorf("invalid priority value at line %d: %s", lineNum, parts[1])
		}

		if _, exists := pSet[priority]; exists {
			return nil, fmt.Errorf("duplicate priority value at line %d: %d", lineNum, priority)
		}

		pMap[parts[0]] = priority
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading file: %w", err)
	}

	return pMap, nil
}
