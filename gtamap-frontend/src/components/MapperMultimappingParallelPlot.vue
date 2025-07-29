<template>
  <w-button @click="toggleFilter" style="margin-bottom: 10px;">
    {{ filterEnabled ? 'Show concordant results' : 'Hide concordant results' }}
  </w-button>
  <div ref="mapperMultimappingParallelContainer"></div>
</template>

<script>
import { onMounted, ref, inject } from 'vue'
import * as d3 from 'd3'
import { sankey, sankeyLinkHorizontal } from 'd3-sankey'

export default {
  setup() {

    const ApiService = inject("http")

    const mapperMultimappingParallelContainer = ref(null)

    const filterEnabled = ref(false)

    function toggleFilter() {
      filterEnabled.value = !filterEnabled.value

      // Clear previous chart
      d3.select(mapperMultimappingParallelContainer.value).selectAll('svg').remove()

      // Redraw chart with new filter setting
      draw()
    }

    async function draw() {

      const margin = { top: 40, right: 20, bottom: 20, left: 20 }
      const width = mapperMultimappingParallelContainer.value.clientWidth
      const height = 500

      try {

        // const data = [
        //   { datapoint: 'A', method1: 1, method2: 0, method3: 1 },
        //   { datapoint: 'B', method1: 2, method2: 2, method3: 2 },
        //   { datapoint: 'C', method1: 0, method2: 1, method3: 1 },
        //   { datapoint: 'D', method1: 1, method2: 1, method3: 0 },
        // ]

        // const dimensions = ['method1', 'method2', 'method3']

        const response = await ApiService.get("/api/mapperMultimappingParallel")

        // transform the response data to the format needed for the sankey diagram
        let data = response.data.map(d => {

          let res = {
            datapoint: d.qname
          }

          for (const mapper in d.mappers) {
            res[mapper] = d.mappers[mapper]
          }

          return res
        })

        // Apply optional filtering
        if (filterEnabled.value) {
          data = data.filter(d => {
            for (const val in d) {
              if (val !== 'datapoint' && d[val] !== 2) {
                return true
              }
            }
            return false
          })
        }

        // extract the dimensions (mappers) from the data
        const dimensions = Object.keys(data[0]).filter(d => d !== 'datapoint')
        dimensions.sort((a, b) => {
          if (a === "gtamap") return -1
          if (b === "gtamap") return 1
          return a.localeCompare(b)
        })

        const width = mapperMultimappingParallelContainer.value.clientWidth
        const height = 500

        const svg = d3
            .select(mapperMultimappingParallelContainer.value)
            .append('svg')
            .attr('width', width)
            .attr('height', height)

        const sankeyData = []

        data.forEach(row => {
          for (let i = 0; i < dimensions.length - 1; i++) {
            const source = `${dimensions[i]}:${row[dimensions[i]]}`
            const target = `${dimensions[i + 1]}:${row[dimensions[i + 1]]}`
            const key = `${source}->${target}`
            const existing = sankeyData.find(d => d.key === key)
            if (existing) {
              existing.value += 1
            } else {
              sankeyData.push({
                source,
                target,
                key,
                value: 1,
              })
            }
          }
        })

        const nodesSet = new Set()
        sankeyData.forEach(d => {
          nodesSet.add(d.source)
          nodesSet.add(d.target)
        })
        const nodes = Array.from(nodesSet).map((id, index) => ({ id, index }))
        const nodeIndex = Object.fromEntries(nodes.map(n => [n.id, n.index]))
        const links = sankeyData.map(d => ({
          source: nodeIndex[d.source],
          target: nodeIndex[d.target],
          value: d.value,
        }))

        const sankeyLayout = sankey()
            .nodeWidth(15)
            .nodePadding(10)
            .extent([
              [0, 0],
              [width - margin.left - margin.right, height - margin.top - margin.bottom],
            ])

        const graph = sankeyLayout({
          nodes: nodes.map(d => Object.assign({}, d)),
          links: links.map(d => Object.assign({}, d)),
        })

        const g = svg.append('g')
            .attr('transform', `translate(${margin.left},${margin.top})`)

        // nodes
        g.append('g')
            .selectAll('rect')
            .data(graph.nodes)
            .join('rect')
            .attr('x', d => d.x0)
            .attr('y', d => d.y0)
            .attr('height', d => d.y1 - d.y0)
            .attr('width', d => d.x1 - d.x0)
            .attr('fill', '#69b3a2')
            .append('title')
            .text(d => `${d.id}\n${d.value}`)

        // labels for each node that show the value represented by the node
        g.append('g')
            .selectAll('text')
            .data(graph.nodes)
            .join('text')
            .attr('x', d => d.x0 < (width / 2) ? d.x1 + 6 : d.x0 - 6)
            .attr('y', d => (d.y0 + d.y1) / 2)
            .attr('dy', '0.35em')
            .attr('text-anchor', d => d.x0 < (width / 2) ? 'start' : 'end')
            .attr('font-size', 10)
            .text(d => d.id.split(':')[1])

        // links between nodes
        g.append('g')
            .attr('fill', 'none')
            .attr('stroke-opacity', 0.5)
            .selectAll('g')
            .data(graph.links)
            .join('path')
            .attr('d', sankeyLinkHorizontal())
            .attr('stroke', '#000')
            .attr('stroke-width', d => Math.max(1, d.width))
            .append('title')
            .text(d => `${d.source.id} → ${d.target.id}\n${d.value}`)

        // title
        // svg.append('text')
        //     .attr('x', width / 2)
        //     .attr('y', margin.top / 2)
        //     .attr('text-anchor', 'middle')
        //     .attr('font-size', '12px')
        //     .attr('font-weight', 'bold')
        //     .text('Number of Mappings per Read and Mapper')

        // subtitle
        svg.append("text")
            .attr('x', width / 2)
            .attr('y', margin.top / 2)
            .attr("text-anchor", "middle")
            .style("font-size", "11px")
            .style("fill", "grey")
            .style("max-width", 400)
            .text("Number of mapping locations per read and mapper");

        const methodPositions = {}

        dimensions.forEach(method => {
          // Find all nodes matching this method (they're named like 'method1:2')
          const matchingNodes = graph.nodes.filter(n => n.id.startsWith(`${method}:`))

          if (matchingNodes.length > 0) {
            const avgX = d3.mean(matchingNodes, d => (d.x0 + d.x1) / 2)
            methodPositions[method] = avgX
          }
        })

        // Add method name labels above each column
        svg.append('g')
            .selectAll('text')
            .data(dimensions)
            .join('text')
            .attr('x', d => margin.left + methodPositions[d] -100)
            .attr('y', margin.top -1)
            .attr('text-anchor', 'middle')
            .attr('font-size', 11)
            .attr('font-weight', 'bold')
            .attr('fill', 'black')
            // rotate 90 deg
            .attr('transform', d => `rotate(-90, ${margin.left + methodPositions[d]}, ${margin.top - 4})`)
            .text(d => d)

      } catch (error) {
        console.error('Error loading data:', error)
      }

    }

    onMounted(() => {
      draw()
    })

    return {
      mapperMultimappingParallelContainer,
      filterEnabled,
      toggleFilter
    }
  },
}
</script>
