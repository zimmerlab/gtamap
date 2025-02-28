package keywordtreebyte

import (
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core/mapper/matchutils"
	"github.com/sirupsen/logrus"
)

type Position struct {
	SequenceIndex uint32
	Position      uint32
}

type Node struct {
	Id        uint32
	Edges     map[byte]uint32
	Positions []Position
}

type Edge struct {
	To    uint32
	Label byte
}

type KeywordTree struct {
	KeywordLength uint8            // the length of every keyword in this tree
	NumSequences  uint8            // the number of sequences in the tree
	RootId        uint32           // the id of the root node
	Nodes         map[uint32]*Node // all nodes in the tree, key = node id
}

func NewKeywordTree(keywordLength uint8) *KeywordTree {
	rootNode := &Node{
		Id:    0,
		Edges: make(map[byte]uint32),
	}

	return &KeywordTree{
		NumSequences:  0,
		KeywordLength: keywordLength,
		RootId:        0,
		Nodes:         map[uint32]*Node{0: rootNode},
	}
}

// AddKeyword adds a keyword to the tree.
// The keyword must have exactly the length as specified in the tree.
// For every character in the keyword either the existing edge is followed or a new edge is created.
// The last node on which this keyword ends is returned.
// The position of the keyword in the sequence is not added in this function but must be added afterward.
func (t *KeywordTree) AddKeyword(keyword []byte) *Node {

	if len(keyword) != int(t.KeywordLength) {
		logrus.WithFields(logrus.Fields{
			"keyword": keyword,
		}).Error("Keyword has invalid length")
		panic(fmt.Sprintf("Keyword has invalid length: %s", keyword))
	}

	logrus.WithFields(logrus.Fields{
		"keyword": keyword,
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

	logrus.Debug("Keyword was added to the tree")

	return activeNode
}

func (t *KeywordTree) FindKeyword(keyword *[]byte, posInRead int) []*matchutils.Match {

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

	result := make([]*matchutils.Match, len(activeNode.Positions))

	for i, position := range activeNode.Positions {

		result[i] = &matchutils.Match{
			SequenceIndex: int(position.SequenceIndex),
			FromGenome:    int(position.Position),
			ToGenome:      int(position.Position + uint32(len(*keyword))),
			FromRead:      posInRead,
			ToRead:        posInRead + int(t.KeywordLength),
			StartGenome:   int(position.Position) - posInRead,
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
