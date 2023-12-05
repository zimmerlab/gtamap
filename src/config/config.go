package config

var env string = "development"

const toolVersion string = "0.1"

var kmerLength int = 8
var numKmers int = 8

func Env() string {
	return env
}

func ToolVersion() string {
	return toolVersion
}

func KmerLength() int {
	return kmerLength
}

func NumKmers() int {
	return numKmers
}
