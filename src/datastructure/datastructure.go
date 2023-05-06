package datastructure

import (
	"container/list"
	"fmt"
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

func (tree *SuffixTree) PrintEdgeList() {

	queue := list.New()
	visited := make(map[*Node]bool)

	queue.PushBack(tree.Root)

	for queue.Len() > 0 {

		item := queue.Front()
		queue.Remove(item)
		node := item.Value.(*Node)

		visited[node] = true

		for _, edge := range node.Edges {
			if visited[edge.To] {
				continue
			}
			queue.PushBack(edge.To)
			fmt.Println(node.Id, strconv.Itoa(edge.To.Id), tree.Sequences[edge.SequenceIndex][edge.Start:*edge.End-1])
		}
	}
}

// TODO: check if rule 2 or 3 must be applied even when activeLength == 0
func BuildSuffixTree(sequences []string) *SuffixTree {

	fmt.Println("Building suffix tree...")
	fmt.Println("sequences: ", sequences)

	numNodes := 1

	rootNode := &Node{
		Id:    0,
		Edges: make(map[byte]*Edge),
		Link:  nil,
	}

	tree := &SuffixTree{
		Root:      rootNode,
		Sequences: sequences,
	}

	for sequenceId, sequence := range tree.Sequences {
		fmt.Println("adding sequence nr ", sequenceId, sequence)

		var activeNode *Node = tree.Root
		var activeEdge *Edge = nil
		var activeLength int = 0
		var remainder = 0

		for stop := 1; stop <= len(sequence); stop++ {

			remainder += 1

			fmt.Println("-------------------")
			fmt.Println("new step: ", stop)
			fmt.Println("suffix: ", sequence[0:stop])
			fmt.Println("--")
			fmt.Println("active node: ", activeNode)
			fmt.Println("active edge: ", activeEdge)
			fmt.Println("active length: ", activeLength)
			fmt.Println("remainder: ", remainder)
			fmt.Println("--")

			// stores the last node which splits and edge that was added in this step
			var lastSplitNode *Node = nil

			// determines and adds actual suffixes to add which are not added implicitly
		addremainder:
			//for i := 0; i < remainderTmp; i++ {
			for remainder > 0 {

				// the current suffix to add
				suffix := sequence[stop-remainder : stop]
				charToAdd := suffix[len(suffix)-1]

				fmt.Println("in remainder addition: ", remainder)
				fmt.Println("suffix: ", suffix)
				fmt.Println("char to add: ", string(charToAdd))

				if activeLength == 0 {
					// no edge must be split
					// it must be checked if a new node+edge must be inserted or if edge with suffix exists

					fmt.Println("active length == 0 -> check current node for path with ", string(charToAdd))

					// check if path with char exists
					// -> if yes, set active edge and active length
					if activeNode.Edges[charToAdd] != nil {

						fmt.Println("path exists -> update active edge, length and remainder")

						// check if path ends after this char -> update active node, reset active edge
						if *activeNode.Edges[charToAdd].End-activeNode.Edges[charToAdd].Start == activeLength+1 {

							fmt.Println("path ends after this char -> update active node, reset active edge")

							activeNode = activeNode.Edges[charToAdd].To
							activeEdge = nil
							activeLength = 0
							//remainder += 1
							break addremainder
						}

						fmt.Println("path does not end after this char -> update active edge, length and remainder")

						activeEdge = activeNode.Edges[charToAdd]
						activeLength = 1
						//remainder += 1
						break addremainder
					}

					fmt.Println("path does not exist -> add new node and edge")

					// edge with suffix does not exist -> add new node and edge
					newNode := &Node{
						Id:    numNodes,
						Edges: make(map[byte]*Edge),
						Link:  nil,
					}
					numNodes += 1
					newEdge := &Edge{
						SequenceIndex: sequenceId,
						Start:         stop - remainder,
						End:           &stop,
						To:            newNode,
					}
					activeNode.Edges[suffix[0]] = newEdge
				} else {
					// current position is on edge
					// check if edge must be split or if suffix is already on edge

					fmt.Println("active length > 0 -> check if edge must be split or if suffix is already on edge")

					nextCharOnEdge := tree.Sequences[activeEdge.SequenceIndex][activeEdge.Start+activeLength]

					fmt.Println("next char on edge: ", string(nextCharOnEdge))

					if nextCharOnEdge == charToAdd {

						fmt.Println("next char is on edge, check if path ends after this char..")

						// check if path ends after this char -> update active node, reset active edge
						if *activeNode.Edges[charToAdd].End-activeNode.Edges[charToAdd].Start == activeLength+1 {

							fmt.Println("path ends after this char -> update active node, reset active edge")

							activeNode = activeNode.Edges[charToAdd].To
							activeEdge = nil
							activeLength = 0
							//remainder += 1
							break addremainder
						}

						fmt.Println("path does not end after this char -> update active length and remainder")

						// suffix is already on edge -> update active length and remainder
						activeLength += 1
						//remainder += 1
						break addremainder
					}

					fmt.Println("next char on edge != charToAdd -> split edge")

					// create new node that splits the activeEdge
					newSplitNode := &Node{
						Id:    numNodes,
						Edges: make(map[byte]*Edge),
						Link:  nil,
					}
					numNodes += 1

					fmt.Println("added new node that splits the active edge id: ", newSplitNode.Id)

					// edge between the newly inserted split node and the old end node
					newSplitEdge := &Edge{
						SequenceIndex: sequenceId,
						Start:         activeEdge.Start + activeLength,
						End:           &stop,
						To:            activeEdge.To,
					}
					newSplitNode.Edges[nextCharOnEdge] = newSplitEdge

					fmt.Println("added new edge that connects the split node to the old end node id: ", newSplitEdge.To.Id)

					// updates the active edge to point to the new split node
					newEnd := activeEdge.Start + activeLength + 1
					activeEdge.End = &newEnd
					activeEdge.To = newSplitNode

					fmt.Println("updated the old edge to end at the new split node: ", activeEdge.Start, *activeEdge.End)
					fmt.Println(activeEdge.Start, *activeEdge.End)
					fmt.Println(sequences[activeEdge.SequenceIndex][activeEdge.Start:*activeEdge.End])

					// create new node and edge for the new suffix
					newNode := &Node{
						Id:    numNodes,
						Edges: make(map[byte]*Edge),
						Link:  nil,
					}
					numNodes += 1

					fmt.Println("added new node for the new suffix edge id: ", newNode.Id)

					// new edge between the new split node and the new node
					newEdge := &Edge{
						SequenceIndex: sequenceId,
						Start:         stop - 1,
						End:           &stop,
						To:            newNode,
					}
					newSplitNode.Edges[charToAdd] = newEdge

					// check if rule 1 applies:
					// active node is root: active node = root, active edge = first char of new suffix, active length -= 1
					if activeNode == tree.Root {
						fmt.Println("rule 1 applies")

						fmt.Println("-> active edge: ", string(suffix[1]))
						fmt.Println("-> active length -1 to ", activeLength-1)

						activeEdge = activeNode.Edges[suffix[1]]
						activeLength -= 1
					}

					// check if rule 2 applies:
					// inserted split node is not the first node added in this step, links previous node to it
					if lastSplitNode != nil {
						fmt.Println("rule 2 applies")

						lastSplitNode.Link = newSplitNode
						fmt.Println("added link from last split node to new split node")
					}

					// check if rule 3 applies:
					// active node is not root:
					// if active node has no suffix link -> active node = root
					// if active node has suffix link -> set active node to suffix link node
					if activeNode != tree.Root {

						fmt.Println("rule 3 applies")

						if activeNode.Link == nil {
							activeNode = tree.Root
							fmt.Println("active node has no suffix link -> set active node to root")
						} else {
							activeNode = activeNode.Link
							fmt.Println("active node has suffix link -> set active node to suffix link node")
						}
					}

					// keeps track of the currently added split node
					lastSplitNode = newSplitNode
				}

				remainder -= 1
			}

		}

		fmt.Println("###################")
		break
	}

	return tree
}
