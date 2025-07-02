<template>
  <div ref="treeContainer" class="upgma-tree"></div>
</template>

<script>

import { agnes } from 'ml-hclust'

import * as d3 from 'd3'
import { ref, onMounted, inject } from "vue";

export default {
  name: "MapperSimilarityUpgmaPlot",
  setup() {

    const ApiService = inject("http")

    const treeContainer = ref(null)

    async function draw() {

      try {

        const response = await ApiService.get("/api/mapperDistances")
        console.log(response.data)

        const distanceMatrix = response.data.distances;
        const labels = response.data.mappers;

        const data = agnes(distanceMatrix, {
          method: 'upgma',
          labels: labels
        });

        console.log(data)

        const container = treeContainer.value
        if (!container) {
          return
        }

        // Clear previous content
        d3.select(container).selectAll("*").remove()

        const margin = { top: 20, right: 120, bottom: 20, left: 40 }
        const width = 800 - margin.left - margin.right
        const height = 400 - margin.top - margin.bottom

        const svg = d3.select(container)
            .append('svg')
            .attr('width', width + margin.left + margin.right)
            .attr('height', height + margin.top + margin.bottom)

        const g = svg.append('g')
            .attr('transform', `translate(${margin.left},${margin.top})`)

        // Custom layout function that uses UPGMA heights
        function layoutTree(node, maxHeight, leafCount) {
          const leaves = []

          function collectLeaves(n) {
            if (n.isLeaf || n.children.length === 0) {
              leaves.push(n)
            } else {
              n.children.forEach(collectLeaves)
            }
          }
          collectLeaves(node)

          const leafSpacing = height / (leaves.length - 1 || 1)

          function assignPositions(n, leafIndex = { value: 0 }) {
            if (n.isLeaf || n.children.length === 0) {
              n.x = leafIndex.value * leafSpacing
              n.y = width // Leaves at the right edge
              leafIndex.value++
            } else {
              // Process children first
              n.children.forEach(child => assignPositions(child, leafIndex))

              // Position internal node
              n.x = d3.mean(n.children, d => d.x) // Average of children's x positions
              n.y = width * (1 - n.height / maxHeight) // Use height for horizontal positioning
            }
          }

          assignPositions(node)
        }

        // Find maximum height for scaling
        function findMaxHeight(node) {
          let maxHeight = node.height
          if (node.children && node.children.length > 0) {
            node.children.forEach(child => {
              maxHeight = Math.max(maxHeight, findMaxHeight(child))
            })
          }
          return maxHeight
        }

        const maxHeight = findMaxHeight(data)
        layoutTree(data, maxHeight)

        // Collect all nodes for D3 selection
        function collectNodes(node, nodes = []) {
          nodes.push(node)
          if (node.children && node.children.length > 0) {
            node.children.forEach(child => collectNodes(child, nodes))
          }
          return nodes
        }

        function collectLinks(node, links = []) {
          if (node.children && node.children.length > 0) {
            node.children.forEach(child => {
              links.push({ source: node, target: child })
              collectLinks(child, links)
            })
          }
          return links
        }

        const nodes = collectNodes(data)
        const links = collectLinks(data)

        // Create links (branches) - horizontal then vertical lines for dendrogram style
        const link = g.selectAll('.link')
            .data(links)
            .enter()
            .append('path')
            .attr('class', 'link')
            .attr('d', d => {
              // Create right-angle connections typical of dendrograms
              return `M${d.source.y},${d.source.x}
                  L${d.target.y},${d.source.x}
                  L${d.target.y},${d.target.x}`
            })
            .style('fill', 'none')
            .style('stroke', '#555')
            .style('stroke-width', '2px')

        // Create nodes
        const node = g.selectAll('.node')
            .data(nodes)
            .enter()
            .append('g')
            .attr('class', 'node')
            .attr('transform', d => `translate(${d.y},${d.x})`)

        // Add circles for nodes
        node.append('circle')
            .attr('r', 4)
            .style('fill', d => (d.children && d.children.length > 0) ? '#555' : '#999')
            .style('stroke', '#fff')
            .style('stroke-width', '2px')

        // Add labels for leaf nodes
        node.filter(d => d.isLeaf || (d.children && d.children.length === 0))
            .append('text')
            .attr('dx', 8)
            .attr('dy', 3)
            .style('text-anchor', 'start')
            .style('font-size', '12px')
            .style('font-family', 'Arial, sans-serif')
            .text(d => d.index >= 0 ? labels[d.index] : `Leaf ${d.index}`)

        // Add height labels on the vertical segments of branches
        const heightLabels = g.selectAll('.height-label')
            .data(links)
            .enter()
            .append('text')
            .attr('class', 'height-label')
            .attr('x', d => d.target.y - 5)
            .attr('y', d => d.source.x - 5)
            .style('text-anchor', 'end')
            .style('font-size', '10px')
            .style('font-family', 'Arial, sans-serif')
            .style('fill', '#666')
            .text(d => d.target.height.toFixed(3))

        // const response = await ApiService.get("/api/mapperDistances")
        //
        // const distanceMatrix = response.data.distances;
        // const labels = response.data.mappers;
        //
        // const data = agnes(distanceMatrix, {
        //   method: 'upgma', // or 'wpgma'
        //   labels: labels
        // });
        //
        // const container = treeContainer.value
        // if (!container) {
        //   return
        // }
        //
        // // Clear previous content
        // d3.select(container).selectAll("*").remove()
        //
        // const margin = { top: 20, right: 120, bottom: 20, left: 40 }
        // const width = 800 - margin.left - margin.right
        // const height = 400 - margin.top - margin.bottom
        //
        // const svg = d3.select(container)
        //     .append('svg')
        //     .attr('width', width + margin.left + margin.right)
        //     .attr('height', height + margin.top + margin.bottom)
        //
        // const g = svg.append('g')
        //     .attr('transform', `translate(${margin.left},${margin.top})`)
        //
        // // Create hierarchy from clustering data
        // const root = d3.hierarchy(data)
        //
        // // Create cluster layout (dendrogram)
        // const cluster = d3.cluster()
        //     .size([height, width])
        //
        // const nodes = cluster(root)
        //
        // // Create links (branches)
        // const link = g.selectAll('.link')
        //     .data(nodes.links())
        //     .enter()
        //     .append('path')
        //     .attr('class', 'link')
        //     .attr('d', d => {
        //       return `M${d.source.y},${d.source.x}
        //             C${(d.source.y + d.target.y) / 2},${d.source.x}
        //             ${(d.source.y + d.target.y) / 2},${d.target.x}
        //             ${d.target.y},${d.target.x}`
        //     })
        //     .style('fill', 'none')
        //     .style('stroke', '#555')
        //     .style('stroke-width', '2px')
        //
        // // Create nodes
        // const node = g.selectAll('.node')
        //     .data(nodes.descendants())
        //     .enter()
        //     .append('g')
        //     .attr('class', 'node')
        //     .attr('transform', d => `translate(${d.y},${d.x})`)
        //
        // // Add circles for nodes
        // node.append('circle')
        //     .attr('r', 4)
        //     .style('fill', d => d.children ? '#555' : '#999')
        //     .style('stroke', '#fff')
        //     .style('stroke-width', '2px')
        //
        // // Add labels for leaf nodes
        // node.filter(d => !d.children)
        //     .append('text')
        //     .attr('dx', 8)
        //     .attr('dy', 3)
        //     .style('text-anchor', 'start')
        //     .style('font-size', '12px')
        //     .style('font-family', 'Arial, sans-serif')
        //     .text(d => d.data.index >= 0 ? labels[d.data.index] : d.data.index)
        //
        // // Add distance labels on branches
        // const linkLabels = g.selectAll('.link-label')
        //     .data(nodes.links())
        //     .enter()
        //     .append('text')
        //     .attr('class', 'link-label')
        //     .attr('x', d => (d.source.y + d.target.y) / 2)
        //     .attr('y', d => (d.source.x + d.target.x) / 2 + 5)
        //     .style('text-anchor', 'middle')
        //     .style('font-size', '10px')
        //     .style('font-family', 'Arial, sans-serif')
        //     .style('fill', '#666')
        //     .text(d => d.target.data.height.toFixed(2))

      } catch (error) {
        console.error('Error loading data:', error)
      }
    }

    onMounted(() => {
      draw()
    })

    return {
      draw,
      treeContainer
    }
  }
}
</script>

<style scoped>
.upgma-tree {
  width: 100%;
  height: 100%;
  border: 1px solid #ddd;
  border-radius: 6px;
  overflow: hidden;
}
</style>