package utils

import (
	"fmt"
	"net/http"
	"strconv"
)

func GetRequestParams(r *http.Request) map[string]string {
	params := make(map[string]string)
	for key, values := range r.URL.Query() {
		if len(values) > 0 {
			params[key] = values[0]
		} else {
			params[key] = ""
		}
	}
	return params
}

func GetSpecificRequestParam(r *http.Request, key string) (string, error) {
	values := r.URL.Query()[key]
	if len(values) > 0 {
		return values[0], nil
	}
	return "", fmt.Errorf("parameter '%s' not found in request", key)
}

func GetSpecificRequestParamBool(r *http.Request, key string) (bool, error) {
	value, err := GetSpecificRequestParam(r, key)
	if err != nil {
		return false, err
	}
	if value != "true" && value != "false" {
		return false, fmt.Errorf("parameter '%s' must be 'true' or 'false', got '%s'", key, value)
	}
	return value == "true", nil
}

func GetSpecificRequestParamInt(r *http.Request, key string) (int, error) {
	value, err := GetSpecificRequestParam(r, key)
	if err != nil {
		return -1, err
	}
	if value == "" {
		return -1, fmt.Errorf("parameter '%s' is empty", key)
	}
	intValue, err := strconv.Atoi(value)
	if err != nil {
		return -1, fmt.Errorf("parameter '%s' must be an integer, got '%s'", key, value)
	}
	return intValue, nil
}
