package datastructure

import (
	"container/list"
	"fmt"
	"github.com/KleinSamuel/gtamap/src/core"
	"github.com/sirupsen/logrus"
	"time"
)

var EndSymbol byte = '$'

type SuffixTree struct {
	// the root node of the suffix tree
	Root *Node
	// the sequences the suffix tree was built from
	Sequences []string
	// (DEBUG) used for printing the edge list
	//Nodes []*Node
	// (DEBUG) used for printing the edge list
	//Edges []*Edge
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
	// (DEBUG) used for printing the edge list
	//From *Node
}

// getEdgeSequence returns the sequence on the given edge
// if the edge is nil, "/" is returned
func (tree *SuffixTree) getEdgeSequence(edge *Edge) string {
	if edge == nil {
		return "/"
	}
	return tree.Sequences[edge.SequenceIndex][edge.Start:*edge.End]
}

// BuildSuffixTree builds a suffix tree from the given sequences using Ukkonen's linear time algorithm
// There are two different cases when adding a new node:
// 1. The new node can be added directly to the active node through an edge carrying only the current char
// -> happens when activeLength == 0 and char is not on any existing edge
// 2. The new node cannot be added directly to the active node, so an edge split is required
// -> happens when activeLength > 0 and char is not on any existing edge
// The 3 addition rules only apply when a new node was inserted using case 2 (edge split):
// Rule 1: If a new node was inserted using case 2 and the active node is not the root node
// -> activeNode stays root, activeLength -= 1, activeEdge = char at (stop - remainder + 1) (next char of current suffix)
// Rule 2: If a new node was inserted using case 2 and it is not the first node inserted during this substep:
// -> connect previously inserted node to the new node using a suffix link
// Rule 3: If a new node was inserted using case 2 and the active node is not the root node:
// -> set activeNode to target of suffix link of activeNode, to root if no suffix link is present
func BuildSuffixTree(sequences []string) *SuffixTree {

	logrus.WithFields(logrus.Fields{
		"sequences": sequences,
	}).Debug("Building suffix tree...")

	timerStart := time.Now()

	for i := 0; i < len(sequences); i++ {
		sequences[i] += string(EndSymbol)
	}

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
		//Nodes:     make([]*Node, 0),
		//Edges:     make([]*Edge, 0),
	}
	//tree.Nodes = append(tree.Nodes, rootNode)

	for sequenceId, sequence := range tree.Sequences {

		logrus.WithFields(logrus.Fields{
			"sequenceId": sequenceId,
			"length":     len(sequence),
		}).Debug("adding new sequence to tree")

		var activeNode *Node = tree.Root
		var activeEdge *Edge = nil
		var activeLength int = 0
		var remainder int = 0
		var stop int = 0

		// the main iteration over the sequence to build all suffixes
		for stop = 1; stop <= len(sequence); stop++ {

			remainder += 1

			logrus.WithFields(logrus.Fields{
				"step":         stop,
				"prefix":       sequence[0:stop],
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
				// the new character (last char of current suffix's prefix) to be added in this step
				charToAdd := suffix[len(suffix)-1]

				// the location of the current suffix within the sequence
				locationToAdd := stop
				if remainder > 1 {
					locationToAdd = stop - remainder + 1
				}

				logrus.WithFields(logrus.Fields{
					"remainder":    remainder,
					"subsequence":  suffix,
					"charToAdd":    string(charToAdd),
					"location":     locationToAdd,
					"step":         stop,
					"activeNode":   activeNode.Id,
					"activeEdge":   tree.getEdgeSequence(activeEdge),
					"activeLength": activeLength,
				}).Debug("in substep")

				if activeLength == 0 {
					// no edge must be split
					// it must be checked if a new node+edge must be inserted or if edge with suffix exists

					logrus.Debug("check current node for path")

					// check if path with char exists
					// -> if yes, set active edge and active length
					if activeNode.Edges[charToAdd] != nil {

						logrus.Debug("path exists")

						activeEdge = activeNode.Edges[charToAdd]

						// check if path ends after this char -> update active node, reset active edge
						if *activeEdge.End-activeEdge.Start == 1 {

							activeNode = activeEdge.To
							activeEdge = nil
							activeLength = 0

							logrus.WithFields(logrus.Fields{
								"activeNode":   activeNode.Id,
								"activeEdge":   tree.getEdgeSequence(activeEdge),
								"activeLength": activeLength,
							}).Debug("path ends on node, update:")

							break addremainder
						}

						activeLength = 1

						logrus.WithFields(logrus.Fields{
							"activeNode":   activeNode.Id,
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

					//tree.Nodes = append(tree.Nodes, newNode)

					logrus.WithFields(logrus.Fields{
						"newNode":   newNode.Id,
						"locations": newNode.Locations,
					}).Debug("created new leaf node")

					newEdge := &Edge{
						SequenceIndex: sequenceId,
						Start:         stop - 1,
						End:           &stop,
						To:            newNode,
						//From:          activeNode,
					}
					activeNode.Edges[sequence[newEdge.Start]] = newEdge

					//tree.Edges = append(tree.Edges, newEdge)

					logrus.WithFields(logrus.Fields{
						"newEdge": tree.getEdgeSequence(newEdge),
						"from":    activeNode.Id,
						"to":      newNode.Id,
						"start":   newEdge.Start,
						"end":     *newEdge.End,
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
						if *activeEdge.End-activeEdge.Start == activeLength+1 {

							activeNode = activeEdge.To
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
							"activeNode":   activeNode.Id,
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

					//tree.Nodes = append(tree.Nodes, newSplitNode)

					logrus.WithFields(logrus.Fields{
						"newSplitNode": newSplitNode.Id,
					}).Debug("created new node that splits the active edge")

					// edge between the newly inserted split node and the old end node
					// contains information from the old edge (careful with sequence index etc.)
					newSplitEdge := &Edge{
						SequenceIndex: activeEdge.SequenceIndex,
						Start:         activeEdge.Start + activeLength,
						End:           activeEdge.End,
						To:            activeEdge.To,
						//From:          newSplitNode,
					}
					// retrieve char from respective sequence of the "old" edge to form new edge
					var edgeStart byte = sequences[newSplitEdge.SequenceIndex][newSplitEdge.Start]
					newSplitNode.Edges[edgeStart] = newSplitEdge

					//tree.Edges = append(tree.Edges, newSplitEdge)

					logrus.WithFields(logrus.Fields{
						"newSplitEdge": tree.getEdgeSequence(newSplitEdge),
						"from":         newSplitNode.Id,
						"to":           newSplitEdge.To.Id,
						"start":        newSplitEdge.Start,
						"end":          *newSplitEdge.End,
					}).Debug("created new edge that connects the split node to the old end node")

					// updates the active edge to point to the new split node
					newEnd := activeEdge.Start + activeLength
					activeEdge.End = &newEnd
					activeEdge.To = newSplitNode

					logrus.WithFields(logrus.Fields{
						"oldEdge": tree.getEdgeSequence(activeEdge),
						"from":    activeNode.Id,
						"to":      activeEdge.To.Id,
						"start":   activeEdge.Start,
						"end":     *activeEdge.End,
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

					//tree.Nodes = append(tree.Nodes, newNode)

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
						//From:          newSplitNode,
					}
					// current sequence is correct because the edge contains the char of the current sequence
					newSplitNode.Edges[sequence[newEdge.Start]] = newEdge

					//tree.Edges = append(tree.Edges, newEdge)

					logrus.WithFields(logrus.Fields{
						"newEdge": tree.getEdgeSequence(newEdge),
						"from":    newSplitNode.Id,
						"to":      newNode.Id,
						"start":   newEdge.Start,
						"end":     *newEdge.End,
					}).Debug("created new edge between the new split node and the new node")

					// check if rule 2 applies:
					// inserted split node is not the first node added in this step, links previous node to it
					if lastSplitNode != nil {

						logrus.Debug("rule 2 applies, but suffix links are ignored for now")

						/*
							lastSplitNode.Link = newSplitNode

							logrus.WithFields(logrus.Fields{
								"from": lastSplitNode.Id,
								"to":   newSplitNode.Id,
							}).Debug("rule 2 applies, create suffix link:")
						*/
					}

					// keeps track of the currently added split node
					lastSplitNode = newSplitNode
				}

				// check if rule 1 applies:
				// active node is root: active node = root, active edge = first char of new suffix, active length -= 1
				if activeNode == tree.Root && activeLength > 0 {

					logrus.Debug("rule 1 applies, computing new activeNode, activeEdge, activeLength")

					activeLength -= 1

					if activeLength == 0 {
						activeEdge = nil
					} else {
						// the active edge is the first char of the new suffix
						// this edge must exist because the activeLength was > 1
						activeEdge = activeNode.Edges[sequence[(stop-remainder)+1]]
					}

					logrus.WithFields(logrus.Fields{
						"activeNode":   activeNode.Id,
						"activeEdge":   tree.getEdgeSequence(activeEdge),
						"activeLength": activeLength,
					}).Debug("rule 1 applied, update:")
				}

				followedSuffixLink := false

				// check if rule 3 applies:
				// a new node was inserted and active node is not root:
				// if active node has no suffix link -> active node = root
				// if active node has suffix link -> set active node to suffix link node
				// set activeEdge to edge starting with same character as current active edge
				// (such an edge must exist because active node was an inner node)
				if activeNode != tree.Root {
					logrus.Debug("rule 3 applies, computing new activeNode, activeEdge, activeLength")

					if activeNode.Link == nil {
						activeNode = tree.Root

						logrus.WithFields(logrus.Fields{
							"activeNode": activeNode.Id,
						}).Debug("active node has no suffix link, update activeNode:")
					} else {
						activeNode = activeNode.Link
						followedSuffixLink = true

						logrus.WithFields(logrus.Fields{
							"activeNode": activeNode.Id,
						}).Debug("active node has suffix link, update activeNode:")
					}
				}

				remainder -= 1

				logrus.WithFields(logrus.Fields{
					"remainder": remainder,
				}).Debug("decrementing remainder")

				rescanCounter := remainder
				if followedSuffixLink {
					rescanCounter = activeLength

					logrus.Debug("setting rescanCounter to activeLength because followed a suffix link")
				}

				if rescanCounter <= 1 {

					//logrus.Debug("setting activeEdge to nil because activeLength == 0")
					logrus.Debug("setting activeEdge to nil because rescanCounter <= 1")

					activeEdge = nil

				} else {

					// the active edge must be updated
					activeEdge = activeNode.Edges[sequence[stop-remainder]]

					// number of chars left to scan
					tmpRemainder := remainder - 1

					logrus.WithFields(logrus.Fields{
						"activeNode":   activeNode.Id,
						"activeEdge":   tree.getEdgeSequence(activeEdge),
						"tmpRemainder": tmpRemainder,
					}).Debug("sub-rescanning tree to find new activePoint")

					for tmpRemainder > 0 && *activeEdge.End-activeEdge.Start <= tmpRemainder {

						logrus.WithFields(logrus.Fields{
							"tmpRemainder": tmpRemainder,
							"edge length":  *activeEdge.End - activeEdge.Start,
							"suffix":       sequence[stop-tmpRemainder-1 : stop],
						}).Debug("in sub-rescan loop, because edge length <= remainder")

						tmpRemainder -= *activeEdge.End - activeEdge.Start
						activeNode = activeEdge.To

						if tmpRemainder == 0 {
							// path ends on a node -> activeEdge = nil
							activeEdge = nil
						} else {
							// path ends on an edge -> activeEdge = edge
							activeEdge = activeNode.Edges[sequence[stop-tmpRemainder-1]]
						}
					}

					activeLength = tmpRemainder

					/* used when suffix links are implemented

					// traverse edge -> activeLength -= edge length
					// update activeNode
					// choose new edge from activeNode that starts with char at this pos in suffix
					// repeat until activeLength < edge length
					for activeLength > 0 && *activeEdge.End-activeEdge.Start <= activeLength {

						activeLength -= *activeEdge.End - activeEdge.Start
						activeNode = activeEdge.To

						newChar := sequence[stop-activeLength-1]

						activeEdge = activeNode.Edges[newChar]

						fmt.Println("new activeLength", activeLength)
						fmt.Println("new activeNode", activeNode.Id)
						fmt.Println("index: ", stop-activeLength-1)
						fmt.Println("suffix: ", sequence[stop-activeLength-1:stop])
						fmt.Println("new char", string(newChar))
						fmt.Println("new activeEdge", tree.getEdgeSequence(activeEdge))

					}
					*/
				}

				logrus.WithFields(logrus.Fields{
					"activeNode":   activeNode.Id,
					"activeLength": activeLength,
					"activeEdge":   tree.getEdgeSequence(activeEdge),
				}).Debug("traversed tree to find new active node, active edge and active length")
			}

		}

		stop -= 1

		// TODO: remove when using multiple sequences
		// break
	}

	logrus.WithFields(logrus.Fields{
		"duration": time.Since(timerStart),
	}).Info("Added all sequences into the suffix tree")

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

// Search searches for the given pattern in the suffix tree and reports the positions
// of occurrences within the sequences of the suffix tree
// Iterates over each char of the given pattern and traverses the tree until the pattern is found or not
// Each leaf node of the tree contains the information about the start position of the subsequence
// (from root to leaf) within the set of sequences as well as the sequence index.
// The pattern is scanned byte by byte and compared to the edges of the tree, which is traversed until
// there are no more bytes (chars) within the pattern. If any byte is not found within the tree, the pattern
// is not contained within the tree. Otherwise, all leaf nodes are returned as occurrences. This is true
// even when the leaf node is not reached because then the pattern is a prefix of a suffix which is still
// contained within the sequences of the tree.
func (tree *SuffixTree) Search(pattern *string) *core.ExactMatchResult {

	logrus.Debug("searching for pattern: ", pattern)

	//result := core.PatternSearchResult{
	//	Pattern: pattern,
	//	Matches: make([]core.PatternMatch, 0),
	//}

	result := core.ExactMatchResult{
		Matches: make([]core.ExactMatch, 0),
	}

	// the last node that was traversed
	var activeNode *Node = tree.Root
	// the currently traversed edge
	var activeEdge *Edge = nil
	// the current position within the active edge
	var indexEdge int = 0

	// iterates over the entire pattern byte by byte
	for indexPattern := 0; indexPattern < len(*pattern); indexPattern++ {

		// the current char of the pattern to be found in the tree
		currentChar := (*pattern)[indexPattern]

		logrus.WithFields(logrus.Fields{
			"indexPattern": indexPattern,
			"currentChar":  string(currentChar),
			"activeNode":   activeNode.Id,
			"activeEdge":   tree.getEdgeSequence(activeEdge),
			"indexEdge":    indexEdge,
		}).Debug("compare next char of pattern")

		if activeEdge == nil {
			// find new active edge

			logrus.Debug("currently no active edge -> find new active edge")

			if activeNode.Edges[currentChar] == nil {
				// pattern not found
				logrus.Debug("no edge with char on current node -> pattern not found")
				return nil
			}

			activeEdge = activeNode.Edges[currentChar]
			indexEdge = 1

			logrus.WithFields(logrus.Fields{
				"activeEdge": tree.getEdgeSequence(activeEdge),
				"indexEdge":  indexEdge,
			}).Debug("update active edge")

		} else {

			// check if char is on active edge
			nextCharOnEdge := tree.Sequences[activeEdge.SequenceIndex][activeEdge.Start+indexEdge]

			logrus.WithFields(logrus.Fields{
				"nextCharOnEdge": string(nextCharOnEdge),
				"currenChar":     string(currentChar),
			}).Debug("compare next char of pattern to next char on edge")

			if nextCharOnEdge != currentChar {
				// pattern not found
				logrus.Debug("char not on active edge -> pattern not found")
				return nil
			}

			indexEdge += 1

			logrus.WithFields(logrus.Fields{
				"indexEdge": indexEdge,
			}).Debug("char found on active edge")
		}

		// check if end of active edge reached -> update active node and reset active edge
		if activeEdge.Start+indexEdge == *activeEdge.End {
			logrus.Debug("end of active edge reached -> update active node and reset active edge")

			activeNode = activeEdge.To
			activeEdge = nil
			indexEdge = 0

			logrus.WithFields(logrus.Fields{
				"activeNode": activeNode.Id,
				"activeEdge": tree.getEdgeSequence(activeEdge),
			}).Debug("end of active edge reached, update:")
		}
	}

	logrus.Debug("pattern found")

	// determines the node from which to start the leaf search
	startNode := activeNode
	if activeEdge != nil {
		// if the current position is not at the end of an edge, the leaf search must start from the next node
		startNode = activeEdge.To
	}

	logrus.WithFields(logrus.Fields{
		"startNode": startNode.Id,
	}).Debug("determining leaf nodes to extract positions")

	leafs := tree.FindLeafNodesRecursive(startNode)

	for _, leaf := range leafs {
		for sequenceId, indices := range leaf.Locations {
			for _, index := range indices {

				logrus.WithFields(logrus.Fields{
					"sequenceId": sequenceId,
					"index":      index,
				}).Debug("pattern match")

				result.Matches = append(result.Matches, core.ExactMatch{
					SequenceIndex: sequenceId,
					FromTarget:    index - 1,
					ToTarget:      index - 1 + len(*pattern),
				})
			}
		}
	}

	return &result
}

/*
// (DEBUG) used for printing the edge list for visualization
func (tree *SuffixTree) ToSif() {

	for _, edge := range tree.Edges {
		fmt.Println(edge.From.Id, tree.getEdgeSequence(edge), edge.To.Id)
	}

	for _, node := range tree.Nodes {
		if node.Link != nil {
			fmt.Println(node.Id, "link", node.Link.Id)
		}
	}
}
*/

func (tree *SuffixTree) ToEdgeList() {

	queue := list.New()
	queue.PushBack(tree.Root)
	visited := make(map[int]bool)

	for queue.Len() > 0 {

		item := queue.Front()
		queue.Remove(item)
		node := item.Value.(*Node)
		visited[node.Id] = true

		for _, edge := range node.Edges {

			fmt.Println(node.Id, edge.To.Id, tree.getEdgeSequence(edge))

			if !visited[edge.To.Id] {
				queue.PushBack(edge.To)
			}
		}
	}
}

/*
// (DEBUG) prints the suffix tree as an edge list for visualization
func (tree *SuffixTree) PrintEdgeList() {

	logrus.Info("Printing edge list")

	fmt.Println("num nodes: ", len(tree.Nodes))
	fmt.Println("num edges: ", len(tree.Edges))

	for _, edge := range tree.Edges {
		fmt.Println(edge.From.Id, edge.To.Id, tree.getEdgeSequence(edge))
	}

	for _, node := range tree.Nodes {
		if node.Link != nil {
			fmt.Println(node.Id, node.Link.Id, "link")
		}
	}
}
*/
