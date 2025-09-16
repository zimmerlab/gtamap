package config

var env string = "development"

const toolVersion string = "0.4.0"

var kmerLength uint8 = 10

var LogOut string = "gta_log.tsv"

func Env() string {
	return env
}

func ToolVersion() string {
	return toolVersion
}

func KmerLength() uint8 {
	return kmerLength
}
