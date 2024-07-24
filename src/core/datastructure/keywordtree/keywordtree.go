package keywordtree

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core"
	"github.com/sirupsen/logrus"
)

type Node struct {
	Id        uint32
	Edges     map[byte]uint32
	Positions []Position
}

type Edge struct {
	To    uint32
	Label byte
}

type Position struct {
	SequenceIndex uint32
	Position      uint32
}

type KeywordTree struct {
	KeywordLength uint8            // the length of every keyword in this tree
	RootId        uint32           // the id of the root node
	Nodes         map[uint32]*Node // all nodes in the tree, key = node id
}

func NewKeywordTree(keywordLength uint8) *KeywordTree {
	rootNode := &Node{
		Id:    0,
		Edges: make(map[byte]uint32),
	}

	return &KeywordTree{
		KeywordLength: keywordLength,
		RootId:        0,
		Nodes:         map[uint32]*Node{0: rootNode},
	}
}

func (t *KeywordTree) AddAllKeywords(word string, sequenceIndex uint32) {
	for i := 0; i < len(word)-int(t.KeywordLength); i++ {
		t.AddKeyword(word[i:i+int(t.KeywordLength)], sequenceIndex, uint32(i))
	}
}

func (t *KeywordTree) AddKeyword(keyword string, sequenceIndex uint32, position uint32) {

	logrus.WithFields(logrus.Fields{
		"keyword":       keyword,
		"sequenceIndex": sequenceIndex,
		"position":      position,
	}).Debug("Adding keyword to tree")

	activeNode := t.Nodes[t.RootId]

	for i := 0; i < len(keyword); i++ {

		if nextNodeId, ok := activeNode.Edges[keyword[i]]; ok {
			// there is an edge with the current character
			// update the active node to the next node
			activeNode = t.Nodes[nextNodeId]
		} else {
			// there is no edge with the current character
			// create a new node and add an edge to it
			newNode := &Node{Id: uint32(len(t.Nodes)), Edges: make(map[byte]uint32)}
			t.Nodes[newNode.Id] = newNode
			activeNode.Edges[keyword[i]] = newNode.Id
			activeNode = newNode
		}
	}

	// add the sequence index of this keyword and the position within the sequence to the last node
	if activeNode.Positions == nil {
		activeNode.Positions = make([]Position, 0)
	}

	activeNode.Positions = append(activeNode.Positions, Position{SequenceIndex: sequenceIndex, Position: position})

	logrus.Debug("Keyword was added to the tree")
}

func (t *KeywordTree) FindKeyword(keyword *string) *core.ExactMatchResult {

	activeNode := t.Nodes[t.RootId]

	for i := 0; i < len(*keyword); i++ {

		if nextNodeId, ok := activeNode.Edges[(*keyword)[i]]; ok {
			// there is an edge with the current character
			// update the active node to the next node
			activeNode = t.Nodes[nextNodeId]
		} else {
			// there is no edge with the current character
			// the keyword is not in the tree
			return nil
		}
	}

	result := &core.ExactMatchResult{
		Matches: make([]*core.SequenceMatch, len(activeNode.Positions)),
	}

	for i, position := range activeNode.Positions {
		result.Matches[i] = &core.SequenceMatch{
			SequenceIndex: int(position.SequenceIndex),
			FromTarget:    int(position.Position),
			ToTarget:      int(position.Position + uint32(len(*keyword))),
		}
	}

	return result
}

func (t *KeywordTree) PrintEdgeList() {

	for nodeId, node := range t.Nodes {
		for label, to := range node.Edges {
			fmt.Println(nodeId, to, string(label))
		}
	}
}
