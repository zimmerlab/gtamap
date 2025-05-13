package graph

import (
	"github.com/KleinSamuel/gtamap/src/core/index"
	"github.com/KleinSamuel/gtamap/src/core/mapper/mappedreadpair"
	"github.com/sirupsen/logrus"
	"strconv"
	"strings"
)

//
// IDEA
//
// 1. As input, we have one or more regions of ambiguity (window)
// 2. For each window create a mapping of pos to readID s.t. we
//    for each pos we know which reads are covering it map[int]int
// 3. Then we iterate from window start to window end. For each pos i
//    in the window, we extract the refPos and all readPositions of reads
//    spanning across pos i
// 4. We create a graph where nodes are readpairs. The interaction between two
//    readpairs is updated for each position using a scoring function until we reach
//    the end of the window
// 5. For the scoring, we give common mismatches a higher score and mismatches
//    inbetween clusters will decrease the weight between nodes

type Node struct {
	Id            int // unique identifier for the readpair
	ReadPairMatch *mappedreadpair.ReadPairMatchResult
	Edges         map[int]*Edge
}

type Edge struct {
	FromNodeId int
	ToNodeId   int
	Weight     int
}

type OverlapGraph struct {
	Nodes     map[int]*Node // adjacency table
	Index     *index.GenomeIndex
	NodeCount int
}

type NodeRegion struct {
	FromReadPair *mappedreadpair.ReadPairMatchResult
	NodeId       int
	FromFw       bool
	BaseAtI      byte
}

func (i *OverlapGraph) InsertNode(readpair *mappedreadpair.ReadPairMatchResult) {
	i.Nodes[i.NodeCount] = &Node{
		Id:            i.NodeCount,
		ReadPairMatch: readpair,
		Edges:         make(map[int]*Edge),
	}
	i.NodeCount += 1
}
func (i *OverlapGraph) InsertEdge(fromId, toId int, weight int) {
	fromNode, fromExists := i.Nodes[fromId]
	_, toExists := i.Nodes[toId]

	if !fromExists || !toExists {
		// Handle error: node does not exist
		logrus.Fatalf("Inserting edge of non existing nodes %s extist: %s, %s exists: %s", fromId, fromExists, toId, toExists)
		return
	}

	if _, exists := fromNode.Edges[toId]; !exists {
		edge := &Edge{
			FromNodeId: fromId,
			ToNodeId:   toId,
			Weight:     weight,
		}
		fromNode.Edges[toId] = edge
	}
}

func (i *OverlapGraph) ToString() string {
	builder := strings.Builder{}
	for id, node := range i.Nodes {

		builder.WriteString("Node " + strconv.Itoa(id) + " Edges: \n")
		for _, edge := range node.Edges {
			builder.WriteString("\t" + strconv.Itoa(edge.FromNodeId) + " -> " + "(" + strconv.Itoa(edge.Weight) + ") -> " + strconv.Itoa(edge.ToNodeId) + "\n")
		}
	}
	return builder.String()
}

func UpdateCluster(og *OverlapGraph, clusteredBaseNodes []int, delta int) {
	for i := 0; i < len(clusteredBaseNodes); i++ {
		for j := i + 1; j < len(clusteredBaseNodes); j++ {
			og.IncrementEdgeWeight(clusteredBaseNodes[i], clusteredBaseNodes[j], delta)
		}
	}
}

func (og *OverlapGraph) IncrementEdgeWeight(fromId, toId, delta int) {
	if fromId == -1 || toId == -1 {
		return
	}

	fromNode, fromExists := og.Nodes[fromId]
	if !fromExists {
		return
	}

	_, toExists := og.Nodes[toId]
	if !toExists {
		return
	}

	edge, exists := fromNode.Edges[toId]
	if exists {
		edge.Weight += delta
	} else {
		fromNode.Edges[toId] = &Edge{
			FromNodeId: fromId,
			ToNodeId:   toId,
			Weight:     delta,
		}
	}
}

type Community struct {
	Id      int
	Members map[int]bool
}

//func (og *OverlapGraph) LouvainClustering() map[int]int {
//	// Phase 1: Optimize modularity by moving individual nodes
//	numNodes := len(og.Nodes)
//	nodeCommunity := make(map[int]int, numNodes)
//	communityTotalWeight := make(map[int]int, numNodes) // Sum of weights of all links inside community
//	nodeDegree := make(map[int]int, numNodes)           // Sum of weights of links incident to node
//	communityDegree := make(map[int]int, numNodes)      // Sum of weights of all links incident to nodes in the community
//	totalEdgeWeight := 0
//
//	// Pre-allocate maps and compute degrees in a single pass
//	for nodeId, node := range og.Nodes {
//		nodeCommunity[nodeId] = nodeId // Each node in its own community
//		communityTotalWeight[nodeId] = 0
//		communityDegree[nodeId] = 0
//
//		selfLoopWeight := 0
//		for _, edge := range node.Edges {
//			nodeDegree[nodeId] += edge.Weight
//			// Count self-loops only once
//			if edge.ToNodeId == nodeId {
//				selfLoopWeight += edge.Weight
//			}
//		}
//
//		communityTotalWeight[nodeId] = selfLoopWeight
//		communityDegree[nodeId] = nodeDegree[nodeId]
//		totalEdgeWeight += nodeDegree[nodeId]
//	}
//
//	totalEdgeWeight /= 2 // Undirected graph: each edge counted twice
//
//	// Precompute a map of node neighbors for faster access
//	nodeNeighbors := make(map[int]map[int]int, numNodes)
//	for nodeId, node := range og.Nodes {
//		nodeNeighbors[nodeId] = make(map[int]int, len(node.Edges))
//		for _, edge := range node.Edges {
//			nodeNeighbors[nodeId][edge.ToNodeId] = edge.Weight
//		}
//	}
//
//	// Track node order and shuffle to avoid local minimum
//	nodeIds := make([]int, 0, numNodes)
//	for nodeId := range og.Nodes {
//		nodeIds = append(nodeIds, nodeId)
//	}
//
//	// Early termination if no improvement for several iterations
//	maxNoImprovementIterations := 3
//	noImprovementCount := 0
//
//	improved := true
//	iterationCount := 0
//	maxIterations := 50 // Set a hard limit on iterations
//
//	for improved && iterationCount < maxIterations {
//		improved = false
//		iterationCount++
//
//		// Shuffle node order for each iteration to avoid getting stuck in local optima
//		rand.Shuffle(len(nodeIds), func(i, j int) {
//			nodeIds[i], nodeIds[j] = nodeIds[j], nodeIds[i]
//		})
//
//		// Track nodes moved in this iteration for faster convergence check
//		nodesMoved := 0
//
//		// Iterate over all nodes in random order
//		for _, nodeId := range nodeIds {
//			currentCommunity := nodeCommunity[nodeId]
//
//			// Fast calculation of community weights for neighbors
//			neighborCommunityWeights := make(map[int]int)
//			for neighborId, weight := range nodeNeighbors[nodeId] {
//				neighborCommunity := nodeCommunity[neighborId]
//				neighborCommunityWeights[neighborCommunity] += weight
//			}
//
//			// Remove node from its current community
//			communityDegree[currentCommunity] -= nodeDegree[nodeId]
//			communityTotalWeight[currentCommunity] -= neighborCommunityWeights[currentCommunity]
//
//			// Find best community with local variables for speed
//			bestCommunity := currentCommunity
//			bestGain := 0.0
//			k_i := nodeDegree[nodeId]
//			m2 := float64(2 * totalEdgeWeight) // Precompute for faster math
//
//			// Only check communities that are actually connected to this node
//			// This is a major optimization for sparse graphs
//			for communityId, k_i_in := range neighborCommunityWeights {
//				sum_tot := communityDegree[communityId]
//
//				// Inline modularity gain calculation for speed
//				dQ := float64(k_i_in)/float64(totalEdgeWeight) - (float64(sum_tot)*float64(k_i))/(m2*float64(totalEdgeWeight))
//
//				if dQ > bestGain {
//					bestGain = dQ
//					bestCommunity = communityId
//				}
//			}
//
//			// Place node in best community
//			nodeCommunity[nodeId] = bestCommunity
//			communityDegree[bestCommunity] += nodeDegree[nodeId]
//			communityTotalWeight[bestCommunity] += neighborCommunityWeights[bestCommunity]
//
//			// Check if we made a change
//			if bestCommunity != currentCommunity {
//				improved = true
//				nodesMoved++
//			}
//		}
//
//		// Early termination check
//		if nodesMoved == 0 || float64(nodesMoved)/float64(numNodes) < 0.001 { // Less than 0.1% of nodes moved
//			noImprovementCount++
//			if noImprovementCount >= maxNoImprovementIterations {
//				break
//			}
//		} else {
//			noImprovementCount = 0
//		}
//	}
//
//	// Phase 2: Aggregate communities into super-nodes (optional for large graphs)
//	if len(og.Nodes) > 1000 { // Only do phase 2 for large graphs
//		// Build new graph where nodes are communities
//		newGraph := NewOverlapGraph()
//		communityToNewId := make(map[int]int)
//		newIdCounter := 0
//
//		// Create new nodes for each community
//		for _, nodeId := range nodeIds {
//			communityId := nodeCommunity[nodeId]
//			if _, exists := communityToNewId[communityId]; !exists {
//				communityToNewId[communityId] = newIdCounter
//				newGraph.AddNode(newIdCounter)
//				newIdCounter++
//			}
//		}
//
//		// Add edges between new community nodes
//		edgesAdded := make(map[string]bool)
//		for nodeId, neighbors := range nodeNeighbors {
//			srcCommunityId := nodeCommunity[nodeId]
//			srcNewId := communityToNewId[srcCommunityId]
//
//			for neighborId, weight := range neighbors {
//				dstCommunityId := nodeCommunity[neighborId]
//				dstNewId := communityToNewId[dstCommunityId]
//
//				// Create a unique edge identifier
//				edgeKey := ""
//				if srcNewId < dstNewId {
//					edgeKey = fmt.Sprintf("%d-%d", srcNewId, dstNewId)
//				} else {
//					edgeKey = fmt.Sprintf("%d-%d", dstNewId, srcNewId)
//				}
//
//				// Add weight to existing edge or create new edge
//				if _, exists := edgesAdded[edgeKey]; !exists {
//					newGraph.AddEdge(srcNewId, dstNewId, weight)
//					edgesAdded[edgeKey] = true
//				} else {
//					// Update existing edge weight
//					for i, edge := range newGraph.Nodes[srcNewId].Edges {
//						if edge.ToNodeId == dstNewId {
//							newGraph.Nodes[srcNewId].Edges[i].Weight += weight
//							break
//						}
//					}
//				}
//			}
//		}
//
//		// Run Louvain on the new graph
//		newCommunities := newGraph.LouvainClustering()
//
//		// Map back to original nodes
//		finalCommunities := make(map[int]int, len(nodeCommunity))
//		for nodeId, communityId := range nodeCommunity {
//			newCommunityId := communityToNewId[communityId]
//			finalCommunityId := newCommunities[newCommunityId]
//			finalCommunities[nodeId] = finalCommunityId
//		}
//
//		return finalCommunities
//	}
//
//	return nodeCommunity
//}
//
//// NewOverlapGraph creates a new overlap graph
//func NewOverlapGraph() *OverlapGraph {
//	return &OverlapGraph{
//		Nodes: make(map[int]*Node),
//	}
//}
//
//// AddNode adds a node to the graph
//func (og *OverlapGraph) AddNode(nodeId int) {
//	og.Nodes[nodeId] = &Node{
//		Id:    nodeId,
//		Edges: make(map[int]*Edge, 0),
//	}
//}
//
//// AddEdge adds an edge between two nodes
//func (og *OverlapGraph) AddEdge(srcId, dstId, weight int) {
//	// Ensure nodes exist
//	if _, exists := og.Nodes[srcId]; !exists {
//		og.AddNode(srcId)
//	}
//	if _, exists := og.Nodes[dstId]; !exists {
//		og.AddNode(dstId)
//	}
//
//	// Add edge
//	og.Nodes[srcId].Edges = append(og.Nodes[srcId].Edges[], &Edge{
//		ToNodeId: dstId,
//		Weight:   weight,
//	})
//
//	// Add reverse edge for undirected graph
//	if srcId != dstId {
//		og.Nodes[dstId].Edges = append(og.Nodes[dstId].Edges, Edge{
//			ToNodeId: srcId,
//			Weight:   weight,
//		})
//	}
//}
//
//// Calculate modularity gain for moving a node to a community
//// k_i_in: Sum of weights of links from node i to nodes in community C
//// sum_tot: Sum of weights of links incident to nodes in community C
//// k_i: Sum of weights of links incident to node i
//// m: Sum of weights of all links in the network
//func modularityGain(k_i_in int, sum_tot int, k_i int, m int) float64 {
//	m2 := float64(2 * m) // Denominator term used twice
//	return float64(k_i_in)/float64(m) - (float64(sum_tot)*float64(k_i))/(m2*float64(m))
//}
