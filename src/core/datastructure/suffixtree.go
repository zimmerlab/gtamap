package datastructure

import (
	"container/list"
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core"
	"github.com/sirupsen/logrus"
	"strconv"
	"strings"
)

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

type Node struct {
	Id            int
	Edges         map[byte]int
	Link          int
	SequenceIndex int
	SequenceStart int
	Positions     []SuffixPosition
}

type Edge struct {
	Label  *SubString
	To     int
	isLast bool // edge connects the node with a leaf node and does not have a label
}

type SubString struct {
	SequenceIndex int
	Start         int
	End           int
}

func (s *SubString) copy() *SubString {
	return &SubString{
		SequenceIndex: s.SequenceIndex,
		Start:         s.Start,
		End:           s.End,
	}
}

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

func (t *SuffixTree) PropagatePositionsToInnerNodes(nodeId int) []SuffixPosition {

	if len(t.GetNode(nodeId).Edges) == 0 {
		return t.GetNode(nodeId).Positions
	}

	var positions []SuffixPosition

	for _, edgeId := range t.GetNode(nodeId).Edges {
		positions = append(positions, t.PropagatePositionsToInnerNodes(t.GetEdge(edgeId).To)...)
	}

	if nodeId != t.RootId {
		t.GetNode(nodeId).Positions = positions
	}

	return positions
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
		Id:            nodeId,
		Edges:         make(map[byte]int),
		Link:          0,
		SequenceIndex: -1,
		SequenceStart: -1,
		Positions:     make([]SuffixPosition, 0),
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
	logrus.Debug("--------------------")

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

		rest = &SubString{
			SequenceIndex: sequenceIndex,
			Start:         i,
			End:           len(sequence),
		}

		logrus.Debug()
		t.ToEdgeList(false)
		logrus.WithFields(logrus.Fields{
			"i":     i,
			"text":  t.GetSubstring(text),
			"char":  string(sequence[i]),
			"rest":  t.GetSubstring(rest),
			"queue": t.PositionQueue,
		}).Debug("adding char to suffix tree")

		t.PositionQueue = append(t.PositionQueue, SuffixPosition{
			Index: sequenceIndex,
			Start: i,
		})

		t.insertedNode = false

		sId, text = t.update(sId, text, sequence[i], rest)

		if sId != 1 && text.length() == 0 {
			fmt.Println("ended on node: ", sId)
			fmt.Println("add position: ", t.PositionQueue[0])

			position := t.PositionQueue[0].copy()
			t.GetNode(sId).Positions = append(t.GetNode(sId).Positions, position)
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

	fmt.Println("position queue")
	fmt.Println(t.PositionQueue)

	if !t.insertedNode {
		logrus.Debug("last suffix is already contained in the tree")
		logrus.Debug("adding positions to node and all of its links")

		activeNodeId := sId

		for activeNodeId != t.RootId {
			logrus.WithFields(logrus.Fields{
				"nodeId":        activeNodeId,
				"sequenceIndex": sequenceIndex,
				"sequenceStart": len(sequence) - t.remainder,
			}).Debug("adding position to node")

			t.GetNode(activeNodeId).Positions = append(t.GetNode(activeNodeId).Positions, SuffixPosition{
				Index: sequenceIndex,
				Start: len(sequence) - t.remainder,
			})
			t.remainder--

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
			// remove the position from the queue
			t.PositionQueue = t.PositionQueue[1:]

			logrus.WithFields(logrus.Fields{
				"removed": position,
				"queue":   t.PositionQueue,
			}).Debug("removing first item from position queue")

			// add the position to the new leaf node
			t.GetNode(leafId).Positions = append(t.GetNode(leafId).Positions, position)

			t.AddEdge(rest, rId, leafId)

			logrus.WithFields(logrus.Fields{
				"leaf":          leafId,
				"sequenceIndex": t.GetNode(leafId).SequenceIndex,
				"sequenceStart": t.GetNode(leafId).SequenceStart,
				"position":      position,
			}).Debug("created new leaf")
		}

		if t.ActiveLeafId != t.RootId {
			// update suffix link for new leaf
			t.GetNode(t.ActiveLeafId).Link = leafId

			logrus.WithFields(logrus.Fields{
				"oldActiveLeaf": t.ActiveLeafId,
				"link":          leafId,
			}).Debug("active leaf != root -> updating suffix link for new leaf")
		}

		t.ActiveLeafId = leafId

		logrus.WithFields(logrus.Fields{
			"activeLeaf": t.ActiveLeafId,
		}).Debug("updating active leaf")

		if oldRootId != t.RootId {
			t.GetNode(oldRootId).Link = rId

			logrus.WithFields(logrus.Fields{
				"oldRoot": oldRootId,
				"link":    rId,
			}).Debug("old root != root -> updating suffix link for old root")
		}
		oldRootId = rId
		logrus.WithFields(logrus.Fields{
			"oldRoot": oldRootId,
		}).Debug("updating old root")

		if t.GetNode(sId).Link == 0 {
			// special case
			k = k.shortenStart()

			logrus.WithFields(logrus.Fields{
				"k":         k,
				"kSequence": t.GetSubstring(k),
			}).Debug("s.Link == nil -> special case")

		} else {

			logrus.Debug("s.Link != nil")

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
		logrus.Debug("searchString: ", t.GetSubstring(searchString))
		logrus.Debug(string(t.GetFirstCharOfSubstring(searchString)))

		edgeId := t.GetNode(startNodeId).Edges[t.GetFirstCharOfSubstring(searchString)]
		edgeSequence := t.GetEdgeLabelStringByEdgeId(edgeId)
		searchSequence := t.GetSubstring(searchString)

		logrus.WithFields(logrus.Fields{
			"edgeSequence":   edgeSequence,
			"searchSequence": searchSequence,
		}).Debug("in testAndSplit !searchString.isEmpty loop")

		// TODO: maybe assert g != null

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

	activeNodeId := t.GetRoot().Id
	activeEdgeId := 0
	activeLength := 0

	for indexPattern := 0; indexPattern < len(*pattern); indexPattern++ {

		logrus.WithFields(logrus.Fields{
			"indexPattern": indexPattern,
			"activeNode":   activeNodeId,
			"activeEdge":   activeEdgeId,
			"activeLength": activeLength,
			"char":         string((*pattern)[indexPattern]),
		}).Debug("find pattern exact loop")

		if activeEdgeId == 0 {

			logrus.Debug("no active edge")

			activeEdgeId = t.GetNode(activeNodeId).Edges[(*pattern)[indexPattern]]

			if activeEdgeId == 0 {
				return result
			}

			logrus.WithFields(logrus.Fields{
				"activeEdge": activeEdgeId,
				"edgeLabel":  t.GetEdgeLabelStringByEdgeId(activeEdgeId),
				"char":       string((*pattern)[indexPattern]),
			}).Debug("found edge with first char of pattern")
		}

		logrus.WithFields(logrus.Fields{
			"activeEdge": activeEdgeId,
			"edgeLabel":  t.GetEdgeLabelStringByEdgeId(activeEdgeId),
		}).Debug("active edge")

		if t.GetEdgeLabelStringByEdgeId(activeEdgeId)[activeLength] != (*pattern)[indexPattern] {
			return result
		}

		activeLength++

		if t.GetEdge(activeEdgeId).Label.length() == activeLength {
			activeNodeId = t.GetEdge(activeEdgeId).To
			activeLength = 0
			activeEdgeId = 0

			logrus.WithFields(logrus.Fields{
				"activeNode":   activeNodeId,
				"activeEdge":   activeEdgeId,
				"activeLength": activeLength,
			}).Debug("end of active edge reached")

			continue
		}
	}

	// determine the node from which to start the search for leaf nodes
	startNodeId := activeNodeId
	if activeEdgeId != 0 {
		// the active edge is not fully traversed, the leaf search must start from its target node
		startNodeId = t.GetEdge(activeEdgeId).To
	}

	leafNodeIds := t.findLeafNodesRecursive(startNodeId)

	for _, leafNodeId := range leafNodeIds {

		leafNode := t.GetNode(leafNodeId)

		result.Matches = append(result.Matches, core.SequenceMatch{
			SequenceIndex: leafNode.SequenceIndex,
			FromTarget:    leafNode.SequenceStart,
			ToTarget:      leafNode.SequenceStart + len(*pattern),
		})
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

func (t *SuffixTree) PrintNodes() {
	for _, node := range t.Nodes {
		fmt.Print(node.Id, ": ")
		posStrings := make([]string, len(node.Positions))
		for i, p := range node.Positions {
			posStrings[i] = strconv.Itoa(p.Index) + "-" + strconv.Itoa(p.Start)
		}
		fmt.Println(strings.Join(posStrings, ","))
	}
}

func (t *SuffixTree) ToEdgeList(printLinks bool) {

	queue := list.New()
	queue.PushBack(t.RootId)
	visited := make(map[int]bool)

	for queue.Len() > 0 {

		item := queue.Front()
		queue.Remove(item)
		nodeId := item.Value.(int)
		visited[nodeId] = true

		for _, edgeId := range t.GetNode(nodeId).Edges {

			fmt.Println(nodeId, t.GetEdge(edgeId).To, t.GetEdgeLabelStringByEdgeId(edgeId))

			if !visited[t.GetEdge(edgeId).To] {
				queue.PushBack(t.GetEdge(edgeId).To)
			}
		}

		if printLinks && t.GetNode(nodeId).Link > 0 {
			fmt.Println(nodeId, t.GetNode(nodeId).Link, "link")
		}
	}
}
