package main

import (
	"github.com/KleinSamuel/gtamap/src/core/datastructure"
	"github.com/KleinSamuel/gtamap/src/test"
	"github.com/sirupsen/logrus"
)

func runTestCases() {

	testCaseFiles := test.GetTestcaseFiles()

	failed := 0

	for _, testCaseFile := range testCaseFiles {

		logrus.Info("Running testcase: ", testCaseFile)

		testCase := test.ReadTestcase(testCaseFile)

		tree := datastructure.CreateTree()

		for i, sequence := range testCase.Sequences {
			tree.AddSequence(sequence, i)
		}

		treeCase := tree.ToTestCase()

		isEqual := testCase.IsEqual(treeCase)

		if !isEqual {
			logrus.Error("Testcase failed!")
			failed++
		} else {
			logrus.Info("Testcase passed!")
		}
	}

	logrus.Info()
	logrus.Info("Summary:\n")
	logrus.Info("Testcases passed: ", len(testCaseFiles)-failed)
	logrus.Info("Testcases failed: ", failed)
}

func main() {
	//testCase := test.ReadTestcase("../resources/test/abc-cef-ceg-ceh.testcase")
	//
	//tree := datastructure.CreateTree()
	//
	//for i, sequence := range testCase.Sequences {
	//	tree.AddSequence(sequence, i)
	//}
	//
	//treeCase := tree.ToTestCase()
	//
	//testCase.IsEqual(treeCase)

	runTestCases()
}
