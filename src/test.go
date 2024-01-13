package main

import "github.com/KleinSamuel/gtamap/src/test"

func main() {
	testCase := test.ReadTestcase("../resources/test/abc-cef-ceg-ceh.testcase")

	testCase2 := test.ReadTestcase("../resources/test/abc-cef-ceg-ceh.testcase")

	testCase.IsEqual(testCase2)
}
