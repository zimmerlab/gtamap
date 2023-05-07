package datastructure

import (
	"container/list"
	"fmt"
	"github.com/sirupsen/logrus"
	"strconv"
)

type SuffixTree struct {
	Root      *Node
	Sequences []string
}

type Node struct {
	// the id of the node
	Id int
	// maps start char of the edge to the edge
	Edges map[byte]*Edge
	// suffix link to other node, nil if no link present
	Link *Node
	// sequence index -> start positions of the substring in the sequence
	Locations map[int][]int
}

type Edge struct {
	// the index of the sequence in the suffix tree sequence array
	SequenceIndex int
	// start position of the sequence on this edge
	Start int
	// end position of the sequence on this edge
	End *int
	// the node the edge points to
	To *Node
}

// PrintEdgeList prints (STDOUT) the suffix tree as an edge list for visualization
func (tree *SuffixTree) PrintEdgeList() {

	queue := list.New()
	visited := make(map[*Node]bool)

	queue.PushBack(tree.Root)

	for queue.Len() > 0 {

		item := queue.Front()
		queue.Remove(item)
		node := item.Value.(*Node)

		visited[node] = true

		// prints the suffix link of the current node
		if node.Link != nil {
			fmt.Println(node.Id, node.Link.Id, "link")
		}

		for _, edge := range node.Edges {
			if visited[edge.To] {
				continue
			}
			queue.PushBack(edge.To)

			// the subsequence on this edge
			edgeSequence := tree.Sequences[edge.SequenceIndex][edge.Start : *edge.End-1]

			// prints the edge between the current node and the node the edge points to
			fmt.Println(node.Id, strconv.Itoa(edge.To.Id), edgeSequence)
		}
	}
}

// getEdgeSequence returns the sequence on the given edge
// if the edge is nil, "/" is returned
func (tree *SuffixTree) getEdgeSequence(edge *Edge) string {
	if edge == nil {
		return "/"
	}
	return tree.Sequences[edge.SequenceIndex][edge.Start : *edge.End-1]
}

// BuildSuffixTree builds a suffix tree from the given sequences using Ukkonen's linear time algorithm
func BuildSuffixTree(sequences []string) *SuffixTree {

	logrus.WithFields(logrus.Fields{
		"sequences": sequences,
	}).Debug("Building suffix tree...")

	// used to assign ids to nodes (not required for production)
	numNodes := 1

	rootNode := &Node{
		Id:        0,
		Edges:     make(map[byte]*Edge),
		Link:      nil,
		Locations: nil,
	}

	tree := &SuffixTree{
		Root:      rootNode,
		Sequences: sequences,
	}

	for sequenceId, sequence := range tree.Sequences {

		logrus.WithFields(logrus.Fields{
			"sequenceId": sequenceId,
			"sequence":   sequence,
		}).Debug("adding new sequence to tree")

		var activeNode *Node = tree.Root
		var activeEdge *Edge = nil
		var activeLength int = 0
		var remainder = 0

		// the main iteration over the sequence to build all suffixes
		for stop := 1; stop <= len(sequence); stop++ {

			remainder += 1

			logrus.WithFields(logrus.Fields{
				"step":         stop,
				"suffix":       sequence[0:stop],
				"activeNode":   activeNode.Id,
				"activeEdge":   tree.getEdgeSequence(activeEdge),
				"activeLength": activeLength,
				"remainder":    remainder,
			}).Debug("new step")

			// stores the last node which splits and edge that was added in this step
			var lastSplitNode *Node = nil

			// determines and adds actual suffixes to add which are not added implicitly
		addremainder:
			for remainder > 0 {

				// the current suffix to add
				suffix := sequence[stop-remainder : stop]
				charToAdd := suffix[len(suffix)-1]

				locationToAdd := stop
				if remainder > 1 {
					locationToAdd = stop - remainder + 1
				}

				logrus.WithFields(logrus.Fields{
					"remainder": remainder,
					"suffix":    suffix,
					"charToAdd": string(charToAdd),
					"location":  locationToAdd,
					"step":      stop,
				}).Debug("in substep")

				if activeLength == 0 {
					// no edge must be split
					// it must be checked if a new node+edge must be inserted or if edge with suffix exists

					logrus.Debug("check current node for path")

					// check if path with char exists
					// -> if yes, set active edge and active length
					if activeNode.Edges[charToAdd] != nil {

						logrus.Debug("path exists")

						// check if path ends after this char -> update active node, reset active edge
						if *activeNode.Edges[charToAdd].End-activeNode.Edges[charToAdd].Start == activeLength+1 {

							activeNode = activeNode.Edges[charToAdd].To
							activeEdge = nil
							activeLength = 0

							logrus.WithFields(logrus.Fields{
								"activeNode":   activeNode.Id,
								"activeEdge":   tree.getEdgeSequence(activeEdge),
								"activeLength": activeLength,
							}).Debug("path ends on node, update:")

							break addremainder
						}

						activeEdge = activeNode.Edges[charToAdd]
						activeLength = 1

						logrus.WithFields(logrus.Fields{
							"activeEdge":   tree.getEdgeSequence(activeEdge),
							"activeLength": activeLength,
						}).Debug("path does not end on node, update:")

						break addremainder
					}

					logrus.Debug("path does not exist -> add new node and edge")

					// edge with suffix does not exist -> add new node and edge
					newNode := &Node{
						Id:        numNodes,
						Edges:     make(map[byte]*Edge),
						Link:      nil,
						Locations: make(map[int][]int),
					}
					numNodes += 1
					newNode.Locations[sequenceId] = []int{locationToAdd}

					logrus.WithFields(logrus.Fields{
						"newNode":   newNode.Id,
						"locations": newNode.Locations,
					}).Debug("created new leaf node")

					newEdge := &Edge{
						SequenceIndex: sequenceId,
						Start:         stop - remainder,
						End:           &stop,
						To:            newNode,
					}
					activeNode.Edges[suffix[0]] = newEdge

					logrus.WithFields(logrus.Fields{
						"newEdge": tree.getEdgeSequence(newEdge),
						"from":    activeNode.Id,
						"to":      newNode.Id,
					}).Debug("created new edge leading to new leaf node")

				} else {
					// current position is on edge
					// check if edge must be split or if suffix is already on edge

					nextCharOnEdge := tree.Sequences[activeEdge.SequenceIndex][activeEdge.Start+activeLength]

					logrus.WithFields(logrus.Fields{
						"activeEdge":     tree.getEdgeSequence(activeEdge),
						"activeLength":   activeLength,
						"charToAdd":      string(charToAdd),
						"nextCharOnEdge": string(nextCharOnEdge),
					}).Debug("check if edge must be split or if char is already on edge")

					if nextCharOnEdge == charToAdd {

						logrus.Debug("next char is on edge")

						// check if path ends after this char -> update active node, reset active edge
						if *activeNode.Edges[charToAdd].End-activeNode.Edges[charToAdd].Start == activeLength+1 {

							activeNode = activeNode.Edges[charToAdd].To
							activeEdge = nil
							activeLength = 0

							logrus.WithFields(logrus.Fields{
								"activeNode":   activeNode.Id,
								"activeEdge":   tree.getEdgeSequence(activeEdge),
								"activeLength": activeLength,
							}).Debug("path ends on node, update:")

							break addremainder
						}

						// suffix is already on edge -> update active length and remainder
						activeLength += 1

						logrus.WithFields(logrus.Fields{
							"activeEdge":   tree.getEdgeSequence(activeEdge),
							"activeLength": activeLength,
						}).Debug("path does not end on node, update:")

						break addremainder
					}

					logrus.Debug("next char is not on edge -> split edge")

					// create new node that splits the activeEdge
					newSplitNode := &Node{
						Id:        numNodes,
						Edges:     make(map[byte]*Edge),
						Link:      nil,
						Locations: nil,
					}
					numNodes += 1

					logrus.WithFields(logrus.Fields{
						"newSplitNode": newSplitNode.Id,
					}).Debug("created new node that splits the active edge")

					// edge between the newly inserted split node and the old end node
					newSplitEdge := &Edge{
						SequenceIndex: sequenceId,
						Start:         activeEdge.Start + activeLength,
						End:           &stop,
						To:            activeEdge.To,
					}
					newSplitNode.Edges[nextCharOnEdge] = newSplitEdge

					logrus.WithFields(logrus.Fields{
						"newSplitEdge": tree.getEdgeSequence(newSplitEdge),
						"from":         newSplitNode.Id,
						"to":           newSplitEdge.To.Id,
					}).Debug("created new edge that connects the split node to the old end node")

					// updates the active edge to point to the new split node
					newEnd := activeEdge.Start + activeLength + 1
					activeEdge.End = &newEnd
					activeEdge.To = newSplitNode

					logrus.WithFields(logrus.Fields{
						"oldEdge": tree.getEdgeSequence(activeEdge),
						"from":    activeNode.Id,
						"to":      activeEdge.To.Id,
					}).Debug("updated the old edge to end at the new split node")

					// create new node and edge for the new suffix
					newNode := &Node{
						Id:        numNodes,
						Edges:     make(map[byte]*Edge),
						Link:      nil,
						Locations: make(map[int][]int),
					}
					numNodes += 1
					newNode.Locations[sequenceId] = []int{locationToAdd}

					logrus.WithFields(logrus.Fields{
						"newNode":   newNode.Id,
						"locations": newNode.Locations,
					}).Debug("created new leaf node")

					// new edge between the new split node and the new node
					newEdge := &Edge{
						SequenceIndex: sequenceId,
						Start:         stop - 1,
						End:           &stop,
						To:            newNode,
					}
					newSplitNode.Edges[charToAdd] = newEdge

					logrus.WithFields(logrus.Fields{
						"newEdge": tree.getEdgeSequence(newEdge),
						"from":    newSplitNode.Id,
						"to":      newNode.Id,
					}).Debug("created new edge between the new split node and the new node")

					// check if rule 1 applies:
					// active node is root: active node = root, active edge = first char of new suffix, active length -= 1
					if activeNode == tree.Root {

						activeEdge = activeNode.Edges[suffix[1]]
						activeLength -= 1

						logrus.WithFields(logrus.Fields{
							"activeEdge":   tree.getEdgeSequence(activeEdge),
							"activeLength": activeLength,
						}).Debug("rule 1 applies, update:")
					}

					// check if rule 2 applies:
					// inserted split node is not the first node added in this step, links previous node to it
					if lastSplitNode != nil {

						lastSplitNode.Link = newSplitNode

						logrus.WithFields(logrus.Fields{
							"from": lastSplitNode.Id,
							"to":   newSplitNode.Id,
						}).Debug("rule 2 applies, create suffix link:")
					}

					// check if rule 3 applies:
					// active node is not root:
					// if active node has no suffix link -> active node = root
					// if active node has suffix link -> set active node to suffix link node
					if activeNode != tree.Root {

						if activeNode.Link == nil {
							activeNode = tree.Root
						} else {
							activeNode = activeNode.Link
						}

						logrus.WithFields(logrus.Fields{
							"activeNode": activeNode.Id,
						}).Debug("rule 3 applies, update:")
					}

					// keeps track of the currently added split node
					lastSplitNode = newSplitNode
				}

				remainder -= 1
			}

		}

		break
	}

	return tree
}

// FindLeafNodesRecursive searches for all leaf nodes in the tree recursively, starting from the given node
func (tree *SuffixTree) FindLeafNodesRecursive(node *Node) []*Node {

	logrus.WithFields(logrus.Fields{
		"node": node.Id,
	}).Debug("search leaf nodes recursively from given node")

	if len(node.Edges) == 0 {
		// found leaf node -> add to list

		logrus.WithFields(logrus.Fields{
			"leaf": node.Id,
		}).Debug("found leaf node")

		return []*Node{node}
	}

	logrus.Debug("is inner node, continue recursively for all edges")

	var leafNodes []*Node

	for _, edge := range node.Edges {

		logrus.WithFields(logrus.Fields{
			"edgeSequence": tree.getEdgeSequence(edge),
			"from":         node.Id,
			"to":           edge.To.Id,
		}).Debug("traversing edge")

		leafNodes = append(leafNodes, tree.FindLeafNodesRecursive(edge.To)...)
	}

	return leafNodes
}

func (tree *SuffixTree) Search(pattern string) {

	logrus.Debug("searching for pattern: ", pattern)

	var activeNode *Node = tree.Root
	var activeEdge *Edge = nil
	var indexEdge int = 0

	for indexPattern := 0; indexPattern < len(pattern); indexPattern++ {

		currentChar := pattern[indexPattern]

		logrus.Debug("searching next char: ", string(currentChar))

		if activeEdge == nil {
			// find new active edge

			logrus.Debug("currently no active edge -> find new active edge")

			if activeNode.Edges[currentChar] == nil {
				// pattern not found
				logrus.Debug("no edge with char on current node -> pattern not found")
				return
			}

			activeEdge = activeNode.Edges[currentChar]
			indexEdge = 1
			logrus.Debug("found new active edge: ", tree.Sequences[activeEdge.SequenceIndex][activeEdge.Start:*activeEdge.End-1])

		} else {

			// check if char is on active edge
			nextCharOnEdge := tree.Sequences[activeEdge.SequenceIndex][activeEdge.Start+indexEdge]

			if nextCharOnEdge != currentChar {
				// pattern not found
				logrus.Debug("char not on active edge -> pattern not found")
				return
			}

			indexEdge += 1

			logrus.Debug("char on active edge -> continue search on active edge")
		}

		// check if end of active edge reached -> update active node and reset active edge
		if activeEdge.Start+indexEdge == *activeEdge.End-1 {
			logrus.Debug("end of active edge reached -> update active node and reset active edge")

			activeNode = activeEdge.To
			activeEdge = nil
		}
	}

	logrus.Debug("pattern found!")
	logrus.Debug("determining leaf nodes..")

	// find leaf nodes
	startNode := activeNode
	if activeEdge != nil {
		startNode = activeEdge.To
	}

	logrus.Debug("starting from node ", startNode.Id)

	leafs := tree.FindLeafNodesRecursive(startNode)

	for _, leaf := range leafs {
		for sequenceId, indices := range leaf.Locations {
			for _, index := range indices {
				logrus.WithFields(logrus.Fields{
					"sequenceId": sequenceId,
					"index":      index,
				}).Info("found pattern")
			}
		}
	}
}
