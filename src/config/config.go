package config

const toolVersion string = "0.4.0"

var kmerLength uint8 = 10

func ToolVersion() string {
	return toolVersion
}

func KmerLength() uint8 {
	return kmerLength
}
