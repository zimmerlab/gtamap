package datastructure

import (
	"container/list"
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core"
	"github.com/KleinSamuel/gtamap/src/test"
	"github.com/sirupsen/logrus"
	"strconv"
	"strings"
)

// This generalized suffix tree is built using linear time Ukkonens Algorithm.
// Original Paper by Ukkonen: http://www.cs.helsinki.fi/u/ukkonen/SuffixT1withFigs.pdf
// Inspired in some ways by: https://github.com/abahgat/suffixtree
// There is no terminal symbol added to the strings.
// Each sequence is referenced by its index in the list of sequences.
// Edges contain only the interval of the string represented by the edge label.
// Such an interval refers to the sequence index and the start and end position of the substring in the sequence.
// Each node contains a list of suffix positions which are the positions of the sequence that end at this node.
// Each node contains the list of suffix positions of all its descendants.
// This way, the position of a certain string can be determined without traversing all the way down to the leaf nodes.

type SuffixPosition struct {
	Index int
	Start int
}

func (s SuffixPosition) copy() SuffixPosition {
	return SuffixPosition{
		Index: s.Index,
		Start: s.Start,
	}
}

// Node represents a node in the suffix tree
type Node struct {
	Id        int              // unique identifier of the node
	Edges     map[byte]int     // edges of the node, key is the first character of the edge label
	Link      int              // suffix link pointing to the id of the linked node
	Positions []SuffixPosition // the list of sequence positions that end at this node (include all descendants)
}

// Edge represents an edge in the suffix tree
type Edge struct {
	Label  *SubString // the label/string on this edge represented by the SubString object
	To     int        // the id of the target node of the edge
	isLast bool       // edge connects the node with a leaf node and does not have a label TODO: check if this can be inferred by To == -1 or 0
}

// SubString represents the substring of a sequence identified by the sequence index and the start and end position
// used to efficiently store strings on the edges of the suffix tree
type SubString struct {
	SequenceIndex int // the index of the sequence in the list of sequences
	Start         int // the start position of the substring in the sequence at the specific sequence index
	End           int // the end position of the substring in the sequence at the specific sequence index
}

func (s *SubString) copy() *SubString {
	return &SubString{
		SequenceIndex: s.SequenceIndex,
		Start:         s.Start,
		End:           s.End,
	}
}

// shortenStart returns a new substring with the start position shortened by one (e.g. +1)
func (s *SubString) shortenStart() *SubString {
	if s.length() == 0 {
		return &SubString{
			SequenceIndex: s.SequenceIndex,
			Start:         s.Start,
			End:           s.End,
		}
	}
	return &SubString{
		SequenceIndex: s.SequenceIndex,
		Start:         s.Start + 1,
		End:           s.End,
	}
}

// shortenEnd returns a new substring with the end position shortened by one (e.g. -1)
func (s *SubString) shortenEnd() *SubString {
	if s.length() == 0 {
		return &SubString{
			SequenceIndex: s.SequenceIndex,
			Start:         s.Start,
			End:           s.End,
		}
	}
	return &SubString{
		SequenceIndex: s.SequenceIndex,
		Start:         s.Start,
		End:           s.End - 1,
	}
}

// extendEnd returns a new substring with the end position extended by one (e.g. +1)
func (s *SubString) extendEnd() *SubString {
	return &SubString{
		SequenceIndex: s.SequenceIndex,
		Start:         s.Start,
		End:           s.End + 1,
	}
}

func (s *SubString) length() int {
	return s.End - s.Start
}
func (s *SubString) isEmpty() bool {
	return s.length() <= 0
}

type SuffixTree struct {
	NumNodes             int
	NumEdges             int
	Sequences            []string
	RootId               int
	ActiveLeafId         int
	currentSequenceIndex int
	currentSequenceStart int
	Nodes                map[int]*Node
	Edges                map[int]*Edge
	PositionQueue        []SuffixPosition
	insertedNode         bool
	remainder            int
}

func CreateTree() *SuffixTree {

	tree := &SuffixTree{
		NumNodes:      1, // 1 is root, 0 is nil
		NumEdges:      1, // 1 is root, 0 is nil
		Nodes:         make(map[int]*Node),
		Edges:         make(map[int]*Edge),
		Sequences:     make([]string, 0),
		PositionQueue: make([]SuffixPosition, 0),
	}

	rootNodeId := tree.AddNode(false)
	tree.RootId = rootNodeId
	tree.ActiveLeafId = tree.RootId

	return tree
}

func (t *SuffixTree) GetPositionsCopyOfNode(nodeId int) []SuffixPosition {
	node := t.GetNode(nodeId)
	positionCopies := make([]SuffixPosition, len(node.Positions))
	for i, pos := range node.Positions {
		positionCopies[i] = pos
	}
	return positionCopies
}

func (t *SuffixTree) RemoveAllSuffixLinks() {
	for _, node := range t.Nodes {
		node.Link = 0
	}
}

func (t *SuffixTree) AddNode(isLeaf bool) int {
	nodeId := t.NumNodes
	t.Nodes[nodeId] = &Node{
		Id:        nodeId,
		Edges:     make(map[byte]int),
		Link:      0,
		Positions: make([]SuffixPosition, 0),
	}

	logrus.WithFields(logrus.Fields{
		"node":   nodeId,
		"isLeaf": isLeaf,
	}).Debug("add new node")

	t.NumNodes++
	return nodeId
}

func (t *SuffixTree) AddEdge(label *SubString, fromId int, toId int) int {

	key := t.GetFirstCharOfSubstring(label)

	logrus.WithFields(logrus.Fields{
		"label": t.GetSubstring(label),
		"from":  fromId,
		"to":    toId,
		"key":   string(key),
	}).Debug("add new edge")

	if t.GetFirstCharOfSubstring(label) != key {
		panic("first char of label does not match key when adding edge")
	}

	edgeId := t.NumEdges
	t.Edges[edgeId] = &Edge{
		Label: label,
		To:    toId,
	}
	t.NumEdges++

	t.GetNode(fromId).Edges[key] = edgeId

	return edgeId
}

func (t *SuffixTree) OverwriteEdge(edgeId int, label *SubString, fromId int, toId int) {

	key := t.GetFirstCharOfSubstring(label)

	logrus.WithFields(logrus.Fields{
		"edgeId": edgeId,
		"label":  t.GetSubstring(label),
		"from":   fromId,
		"to":     toId,
		"key":    string(key),
	}).Debug("overwrite edge")

	if t.GetFirstCharOfSubstring(label) != key {
		panic("first char of label does not match key when adding edge")
	}

	t.Edges[edgeId] = &Edge{
		Label: label,
		To:    toId,
	}

	t.GetNode(fromId).Edges[key] = edgeId
}

func (t *SuffixTree) GetRoot() *Node {
	return t.Nodes[1]
}

func (t *SuffixTree) GetNode(id int) *Node {
	return t.Nodes[id]
}

func (t *SuffixTree) GetEdge(id int) *Edge {
	return t.Edges[id]
}

func (t *SuffixTree) GetEdgeLabelCopyByEdgeId(edgeId int) *SubString {
	return t.Edges[edgeId].Label.copy()
}

// GetEdgeLabelStringByEdgeId returns the label string of the edge with the given id
func (t *SuffixTree) GetEdgeLabelStringByEdgeId(edgeId int) string {
	return t.Sequences[t.Edges[edgeId].Label.SequenceIndex][t.Edges[edgeId].Label.Start:t.Edges[edgeId].Label.End]
}

func (t *SuffixTree) GetSubstring(substring *SubString) string {
	return t.Sequences[substring.SequenceIndex][substring.Start:substring.End]
}

func (t *SuffixTree) GetFirstCharOfSubstring(substring *SubString) byte {
	return t.Sequences[substring.SequenceIndex][substring.Start]
}

func (t *SuffixTree) GetCharAt(sequenceIndex int, position int) byte {
	return t.Sequences[sequenceIndex][position]
}

func (t *SuffixTree) AddSequence(sequence string, sequenceIndex int) {

	//sequence = sequence + strconv.Itoa(sequenceIndex)

	logrus.Debug()
	logrus.Debug("--------------------")
	logrus.WithFields(logrus.Fields{
		"sequenceIndex": sequenceIndex,
		"sequence":      sequence,
	}).Debug("adding sequence to suffix tree")

	// add the new sequence to the list of sequences in the tree for character lookup
	// as each edge contains only an interval to represent the edge label
	t.Sequences = append(t.Sequences, sequence)

	t.ActiveLeafId = t.RootId
	t.PositionQueue = make([]SuffixPosition, 0)

	sId := t.RootId

	var text *SubString = &SubString{
		SequenceIndex: sequenceIndex,
		Start:         0,
		End:           0,
	}
	var rest *SubString

	for i := 0; i < len(sequence); i++ {

		logrus.Debug("--------------------")

		rest = &SubString{
			SequenceIndex: sequenceIndex,
			Start:         i,
			End:           len(sequence),
		}

		//t.PrintEdgeList(true)

		logrus.WithFields(logrus.Fields{
			"i":     i,
			"text":  t.GetSubstring(text),
			"char":  string(sequence[i]),
			"rest":  t.GetSubstring(rest),
			"queue": t.PositionQueue,
		}).Debug("adding char to suffix tree")

		// add the current position to the queue of positions that will be added to the leaf nodes
		t.PositionQueue = append(t.PositionQueue, SuffixPosition{
			Index: sequenceIndex, // the index of the sequence in the list of sequences
			Start: i,             // the start position of the suffix in the sequence
		})

		logrus.WithFields(logrus.Fields{
			"new":   t.PositionQueue[len(t.PositionQueue)-1],
			"queue": t.PositionQueue,
		}).Debug("added new position to queue")

		t.insertedNode = false

		sId, text = t.update(sId, text, sequence[i], rest)

		if sId != 1 && text.length() == 0 {

			position := t.PositionQueue[0].copy()

			// add the current suffix start position to the inner node
			t.GetNode(sId).Positions = append(t.GetNode(sId).Positions, position)

			logrus.WithFields(logrus.Fields{
				"nodeId":   sId,
				"position": position,
			}).Debug("added position to inner node")
		}

		// keep track of the length of the suffix which is already contained in the tree
		if t.insertedNode {
			t.remainder = 0
		} else {
			t.remainder++
		}
	}

	if t.ActiveLeafId != 0 && t.ActiveLeafId != t.RootId && t.ActiveLeafId != sId {
		t.GetNode(t.ActiveLeafId).Link = sId

		logrus.WithFields(logrus.Fields{
			"from": t.ActiveLeafId,
			"to":   sId,
		}).Debug("updating suffix link for active leaf")
	}

	logrus.Debug("done adding sequence to suffix tree")

	if !t.insertedNode {
		logrus.Debug("last suffix is already contained in the tree")
		logrus.Debug("adding positions to all of its links")

		// set the active node to the suffix link of the previous node
		// the previous node has its position set within the main loop
		activeNodeId := t.GetNode(sId).Link
		// remove the position of the previous node from the queue
		t.PositionQueue = t.PositionQueue[1:]

		// add the position to all of its links except the root node
		for activeNodeId != t.RootId {

			position := t.PositionQueue[0].copy()
			// remove the position from the queue
			t.PositionQueue = t.PositionQueue[1:]

			t.GetNode(activeNodeId).Positions = append(t.GetNode(activeNodeId).Positions, position)

			logrus.WithFields(logrus.Fields{
				"nodeId":   activeNodeId,
				"position": position,
			}).Debug("added position to node")

			activeNodeId = t.GetNode(activeNodeId).Link
		}
	}
}

func (t *SuffixTree) update(sId int, stringPart *SubString, nextChar byte, rest *SubString) (int, *SubString) {

	logrus.WithFields(logrus.Fields{
		"s":          sId,
		"stringPart": t.GetSubstring(stringPart),
		"nextChar":   string(nextChar),
		"rest":       t.GetSubstring(rest),
	}).Debug("update")

	k := stringPart.extendEnd()

	oldRootId := t.RootId

	endpoint, rId := t.testAndSplit(sId, stringPart, nextChar, rest)

	logrus.WithFields(logrus.Fields{
		"endpoint": endpoint,
		"r":        rId,
		"k":        t.GetSubstring(k),
		"s":        sId,
	}).Debug("in update start after testAndSplit")

	var leafId int

	for !endpoint {

		logrus.Debug("in update !endpoint loop")

		edgeId := t.GetNode(rId).Edges[nextChar]

		if edgeId != 0 {
			logrus.Debug("r has edge starting with nextChar")

			leafId = t.GetEdge(edgeId).To

			logrus.WithFields(logrus.Fields{
				"leaf": leafId,
			}).Debug("updated leaf")

		} else {
			logrus.Debug("r has no edge starting with nextChar")

			leafId = t.AddNode(true)

			t.insertedNode = true

			position := t.PositionQueue[0].copy()
			// add the position to the new leaf node
			t.GetNode(leafId).Positions = append(t.GetNode(leafId).Positions, position)

			t.AddEdge(rest, rId, leafId)

			logrus.WithFields(logrus.Fields{
				"leaf":     leafId,
				"position": position,
			}).Debug("created new leaf")
		}

		if t.ActiveLeafId != t.RootId {
			// update suffix link for new leaf
			t.GetNode(t.ActiveLeafId).Link = leafId

			logrus.WithFields(logrus.Fields{
				"oldActiveLeaf": t.ActiveLeafId,
				"link":          leafId,
			}).Debug("updated suffix link for new leaf")
		}
		logrus.WithFields(logrus.Fields{
			"old": t.ActiveLeafId,
			"new": leafId,
		}).Debug("updated active leaf")
		t.ActiveLeafId = leafId

		if oldRootId != t.RootId {
			t.GetNode(oldRootId).Link = rId

			logrus.WithFields(logrus.Fields{
				"oldRoot": oldRootId,
				"link":    rId,
			}).Debug("updated suffix link for old root")
		}
		logrus.WithFields(logrus.Fields{
			"old": oldRootId,
			"new": rId,
		}).Debug("updated old root")
		oldRootId = rId

		if t.insertedNode {
			// a node was inserted in this step (split edge and maybe a new leaf)
			// either way this indicates that the current suffix is contained within the tree
			// so the current sequence position can be removed from the queue

			position := t.PositionQueue[0].copy()
			// remove the position from the queue
			t.PositionQueue = t.PositionQueue[1:]

			logrus.WithFields(logrus.Fields{
				"removed": position,
				"queue":   t.PositionQueue,
			}).Debug("removed first item from position queue")
		}

		if t.GetNode(sId).Link == 0 {
			// special case
			k = k.shortenStart()

			logrus.WithFields(logrus.Fields{
				"k":         k,
				"kSequence": t.GetSubstring(k),
			}).Debug("s.Link == nil -> special case")

		} else {

			logrus.Debug("s.Link != nil")
			logrus.WithFields(logrus.Fields{
				"suffixLink": t.GetNode(sId).Link,
			}).Debug("canonize from suffix link")

			//if t.GetNode(sId).Link == 1 {
			//	logrus.Debug("followed suffix link is root")
			//	logrus.Debug("current suffix completely contained")
			//
			//	position := t.PositionQueue[0].copy()
			//	t.PositionQueue = t.PositionQueue[1:]
			//
			//	logrus.WithFields(logrus.Fields{
			//		"queue":   t.PositionQueue,
			//		"removed": position,
			//	}).Debug("remove first position from queue")
			//}

			sId, k = t.canonize(t.GetNode(sId).Link, k.shortenEnd())

			if sId != 1 && k.length() == 0 {

				position := t.PositionQueue[0].copy()
				t.GetNode(sId).Positions = append(t.GetNode(sId).Positions, position)

				logrus.WithFields(logrus.Fields{
					"nodeId":   sId,
					"position": position,
				}).Debug("adding position to inner node because it was traversed")
			}

			k = k.extendEnd()

			logrus.WithFields(logrus.Fields{
				"s": sId,
				"k": t.GetSubstring(k),
			}).Debug("updating s and k after canonize")
		}

		endpoint, rId = t.testAndSplit(sId, k.shortenEnd(), nextChar, rest)
	}

	logrus.Debug("in update !endpoint loop finished")

	if oldRootId != t.RootId {
		t.GetNode(oldRootId).Link = rId

		logrus.WithFields(logrus.Fields{
			"oldRoot": oldRootId,
			"link":    rId,
		}).Debug("old root != root -> updating suffix link for old root")
	}

	return t.canonize(sId, k)
}

func (t *SuffixTree) testAndSplit(startNodeId int, searchString *SubString, nextChar byte, remainder *SubString) (bool, int) {

	logrus.WithFields(logrus.Fields{
		"startNode":            startNodeId,
		"searchString":         searchString,
		"searchStringSequence": t.GetSubstring(searchString),
		"nextChar":             string(nextChar),
		"remainder":            t.GetSubstring(remainder),
	}).Debug("testAndSplit")

	// TODO: maybe assert remainder not empty and remainder char at 0 == nextChar

	startNodeId, searchString = t.canonize(startNodeId, searchString)

	logrus.WithFields(logrus.Fields{
		"startNode":    startNodeId,
		"searchString": t.GetSubstring(searchString),
	}).Debug("in testAndSplit start after canonize")

	if !searchString.isEmpty() {

		logrus.Debug("in testAndSplit !searchString.isEmpty")

		edgeId := t.GetNode(startNodeId).Edges[t.GetFirstCharOfSubstring(searchString)]
		edgeSequence := t.GetEdgeLabelStringByEdgeId(edgeId)
		searchSequence := t.GetSubstring(searchString)

		logrus.WithFields(logrus.Fields{
			"edgeSequence":   edgeSequence,
			"searchSequence": searchSequence,
		}).Debug("in testAndSplit !searchString.isEmpty loop")

		// searchString is a substring of the edge label
		if len(edgeSequence) > len(searchSequence) && edgeSequence[len(searchSequence)] == nextChar {

			logrus.Debug("in testAndSplit !searchString.isEmpty loop -> searchString is a substring of the edge label")

			return true, startNodeId
		}

		newNodeId := t.splitNode(startNodeId, edgeId, searchString)

		return false, newNodeId
	}

	remainderSequence := t.GetSubstring(remainder)

	edgeId := t.GetNode(startNodeId).Edges[remainderSequence[0]]

	if edgeId == 0 {

		logrus.Debug("edgeId == 0")

		// there is no t-transition from s
		return false, startNodeId
	}

	edgeSequence := t.GetEdgeLabelStringByEdgeId(edgeId)

	if strings.Index(edgeSequence, remainderSequence) == 0 {

		logrus.Debug("edgeSequence starts with remainderSequence")

		if len(remainderSequence) == len(edgeSequence) {
			logrus.Debug("remainderSequence == edgeSequence")

			return true, startNodeId

		} else {
			logrus.Debug("remainderSequence != edgeSequence")

			t.splitNode(startNodeId, edgeId, remainder)

			return false, startNodeId
		}
	}

	logrus.Debug("edgeSequence does not start with remainderSequence")

	return true, startNodeId
}

func (t *SuffixTree) splitNode(nodeId int, edgeId int, splitFirstPart *SubString) int {

	logrus.WithFields(logrus.Fields{
		"node":           nodeId,
		"splitFirstPart": t.GetSubstring(splitFirstPart),
		"edgeSequence":   t.GetEdgeLabelStringByEdgeId(edgeId),
	}).Debug("splitNode")

	newNodeId := t.AddNode(false)

	// add all positions of the old target node of the edge to be split to the new inner nodes positions
	t.GetNode(newNodeId).Positions = t.GetPositionsCopyOfNode(t.GetEdge(edgeId).To)
	// get the current position of the suffix being added
	position := t.PositionQueue[0].copy()

	logrus.WithFields(logrus.Fields{
		"inherited": t.GetNode(newNodeId).Positions,
		"new":       position,
	}).Debug("add positions to new inner node")

	// add the current positions to the new inner nodes positions
	t.GetNode(newNodeId).Positions = append(t.GetNode(newNodeId).Positions, position)
	// keep track that a new node was inserted in the current iteration
	t.insertedNode = true

	// new edge that connects the given source node with the new inner node
	newEdgeId := t.AddEdge(splitFirstPart, nodeId, newNodeId)

	// update the given edge such that it connects the new inner node with the destination node of the given edge
	newLabel := t.GetEdgeLabelCopyByEdgeId(edgeId)
	newLabel.Start += splitFirstPart.length()

	t.OverwriteEdge(edgeId, newLabel, newNodeId, t.GetEdge(edgeId).To)

	logrus.WithFields(logrus.Fields{
		"sourceNode":          nodeId,
		"innerNode":           newNodeId,
		"destNode":            t.GetEdge(edgeId).To,
		"newEdgeSequence":     t.GetEdgeLabelStringByEdgeId(newEdgeId),
		"updatedEdgeSequence": t.GetEdgeLabelStringByEdgeId(edgeId),
	}).Debug("splitNode finished")

	return newNodeId
}

func (t *SuffixTree) canonize(startNodeId int, input *SubString) (int, *SubString) {

	logrus.WithFields(logrus.Fields{
		"startNode": startNodeId,
		"input":     t.GetSubstring(input),
	}).Debug("canonize")

	nodeId := startNodeId

	var remainder *SubString = input.copy()

	for !remainder.isEmpty() {

		remainderSequence := t.GetSubstring(remainder)

		logrus.WithFields(logrus.Fields{
			"node":               nodeId,
			"remainder":          remainder,
			"remainderSequence":  remainderSequence,
			"firstCharRemainder": string(remainderSequence[0]),
		}).Debug("remainder not empty")

		edgeId := t.GetNode(nodeId).Edges[remainderSequence[0]]

		if edgeId == 0 {
			logrus.Debug("no edge from current node with first char of remainder")
			break
		}

		edgeSequence := t.GetEdgeLabelStringByEdgeId(edgeId)

		logrus.WithFields(logrus.Fields{
			"edgeSequence": edgeSequence,
		}).Debug("edge from current node starting with first char of remainder")

		if strings.Index(remainderSequence, edgeSequence) != 0 {
			logrus.Debug("remainder sequence does not start with edge sequence")
			break
		}

		nodeId = t.GetEdge(edgeId).To
		remainder.Start += t.GetEdgeLabelCopyByEdgeId(edgeId).length()

		logrus.WithFields(logrus.Fields{
			"node":      nodeId,
			"remainder": t.GetSubstring(remainder),
		}).Debug("updated node and remainder")
	}

	logrus.WithFields(logrus.Fields{
		"node":      nodeId,
		"remainder": t.GetSubstring(remainder),
	}).Debug("canonize finished")

	return nodeId, remainder
}

func (t *SuffixTree) FindPatternExact(pattern *string) *core.ExactMatchResult {

	logrus.WithFields(logrus.Fields{
		"pattern": *pattern,
	}).Debug("find pattern exact")

	result := &core.ExactMatchResult{
		Matches: make([]core.SequenceMatch, 0),
	}

	// the active point is used to point to a specific position in the tree
	// the current node is the last node that was visited
	activeNodeId := t.GetRoot().Id
	// the active edge is the current edge that is being traversed
	activeEdgeId := 0
	// the active length is the number of characters that have been traversed on the active edge (offset)
	activeLength := 0

	for indexPattern := 0; indexPattern < len(*pattern); indexPattern++ {

		//logrus.WithFields(logrus.Fields{
		//	"indexPattern": indexPattern,
		//	"activeNode":   activeNodeId,
		//	"activeEdge":   activeEdgeId,
		//	"activeLength": activeLength,
		//	"char":         string((*pattern)[indexPattern]),
		//}).Debug("find pattern exact loop")

		// determine the active edge if there is none
		if activeEdgeId == 0 {

			//logrus.Debug("currently there is no active edge")

			// find the edge that starts with the character of the pattern at position indexPattern
			activeEdgeId = t.GetNode(activeNodeId).Edges[(*pattern)[indexPattern]]

			if activeEdgeId == 0 {
				return result
			}

			//logrus.WithFields(logrus.Fields{
			//	"activeEdge": activeEdgeId,
			//	"edgeLabel":  t.GetEdgeLabelStringByEdgeId(activeEdgeId),
			//	"char":       string((*pattern)[indexPattern]),
			//}).Debug("found new active edge")
		}

		//logrus.WithFields(logrus.Fields{
		//	"activeEdge": activeEdgeId,
		//	"edgeLabel":  t.GetEdgeLabelStringByEdgeId(activeEdgeId),
		//}).Debug("active edge information")

		// the next character on the current edge does not match the current character in the pattern
		if t.GetEdgeLabelStringByEdgeId(activeEdgeId)[activeLength] != (*pattern)[indexPattern] {

			//logrus.Debug("next character on active edge does not match current character in pattern")

			return result
		}

		//logrus.Debug("char on active edge matches char in pattern")

		// there are additional characters on the string of the current edge
		activeLength++

		// the end of the current edge is reached -> update the active node
		if t.GetEdge(activeEdgeId).Label.length() == activeLength {
			// update the active point information
			activeNodeId = t.GetEdge(activeEdgeId).To
			activeLength = 0
			activeEdgeId = 0

			//logrus.WithFields(logrus.Fields{
			//	"activeNode":   activeNodeId,
			//	"activeEdge":   activeEdgeId,
			//	"activeLength": activeLength,
			//}).Debug("end of active edge reached")

			continue
		}

		//logrus.WithFields(logrus.Fields{
		//	"activeEdge":   activeEdgeId,
		//	"activeLength": activeLength,
		//}).Debug("walking on edge")
	}

	// at least one exact match of the pattern was found at this point

	//logrus.Debug("pattern found")

	// the positions of the pattern match(es) are stored in the next node (contains all positions of all descendants)
	// the next node is the active node if there is no active edge or the active length is 0
	// otherwise the next node is the target node of the active edge
	var targetNode *Node
	if activeEdgeId == 0 || activeLength == 0 {
		targetNode = t.GetNode(activeNodeId)
	} else {
		targetNode = t.GetNode(t.GetEdge(activeEdgeId).To)
	}

	// add all positions of the target node to the result
	for _, position := range targetNode.Positions {

		result.Matches = append(result.Matches, core.SequenceMatch{
			SequenceIndex: position.Index,
			FromTarget:    position.Start,
			ToTarget:      position.Start + len(*pattern),
		})

		logrus.WithFields(logrus.Fields{
			"SequenceIndex": position.Index,
			"FromTarget":    position.Start,
			"ToTarget":      position.Start + len(*pattern),
		}).Debug("determined pattern match")
	}

	return result
}

func (t *SuffixTree) findLeafNodesRecursive(nodeId int) []int {

	if len(t.GetNode(nodeId).Edges) == 0 {
		return []int{nodeId}
	}

	var leafNodeIds []int

	for _, edgeId := range t.GetNode(nodeId).Edges {
		leafNodeIds = append(leafNodeIds, t.findLeafNodesRecursive(t.GetEdge(edgeId).To)...)
	}

	return leafNodeIds
}

func (t *SuffixTree) ToTestCase() *test.Testcase {

	testCase := &test.Testcase{
		Sequences: make([]string, 0),
		Edges:     make([]test.Edge, 0),
		Positions: make(map[int][]test.Position),
	}

	// add all sequences to the test case
	for _, sequence := range t.Sequences {
		testCase.Sequences = append(testCase.Sequences, sequence)
	}

	// compute and add the edge list to the test case
	testCase.Edges = t.toTestEdgeList(true)

	// compute and add the position list to the test case
	testCase.Positions = t.toTestPositionList()

	return testCase
}

func (t *SuffixTree) ListAllNodes() {

	queue := list.New()
	queue.PushBack(t.RootId)
	visited := make(map[int]bool)

	for queue.Len() > 0 {

		item := queue.Front()
		queue.Remove(item)
		nodeId := item.Value.(int)
		visited[nodeId] = true

		node := t.GetNode(nodeId)

		fmt.Println(nodeId, node)

		for _, edgeId := range node.Edges {

			edge := t.GetEdge(edgeId)

			if !visited[edge.To] {
				queue.PushBack(edge.To)
			}
		}
	}
}

func (t *SuffixTree) toTestEdgeList(includeSuffixLinks bool) []test.Edge {

	edgeList := make([]test.Edge, 0)

	queue := list.New()
	queue.PushBack(t.RootId)
	visited := make(map[int]bool)

	for queue.Len() > 0 {

		item := queue.Front()
		queue.Remove(item)
		nodeId := item.Value.(int)
		visited[nodeId] = true

		for _, edgeId := range t.GetNode(nodeId).Edges {

			edgeList = append(edgeList, test.Edge{
				From:  nodeId,
				To:    t.GetEdge(edgeId).To,
				Label: t.GetEdgeLabelStringByEdgeId(edgeId),
			})

			if !visited[t.GetEdge(edgeId).To] {
				queue.PushBack(t.GetEdge(edgeId).To)
			}
		}

		if includeSuffixLinks && t.GetNode(nodeId).Link > 0 {
			edgeList = append(edgeList, test.Edge{
				From:  nodeId,
				To:    t.GetNode(nodeId).Link,
				Label: "link",
			})
		}
	}

	return edgeList
}

func (t *SuffixTree) PrintEdgeList(includeSuffixLinks bool) {

	edgeList := t.toTestEdgeList(includeSuffixLinks)

	for _, edge := range edgeList {
		fmt.Println(edge.From, edge.To, edge.Label)
	}
}

func (t *SuffixTree) toTestPositionList() map[int][]test.Position {

	positionList := make(map[int][]test.Position)

	for _, node := range t.Nodes {

		if len(node.Positions) > 0 {
			positionList[node.Id] = make([]test.Position, 0)
			for _, position := range node.Positions {
				positionList[node.Id] = append(positionList[node.Id], test.Position{
					Index: position.Index,
					Start: position.Start,
				})
			}
		}
	}

	return positionList
}

func (t *SuffixTree) PrintNodes() {

	positionMap := t.toTestPositionList()

	for i, positions := range positionMap {

		fmt.Print(i, ": ")
		posStrings := make([]string, len(positions))
		for i, p := range positions {
			posStrings[i] = strconv.Itoa(p.Index) + "-" + strconv.Itoa(p.Start)
		}
		fmt.Println(strings.Join(posStrings, ","))
	}
}
