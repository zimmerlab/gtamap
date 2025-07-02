<template>
  <div ref="chartContainer" class="pca-chart"></div>
</template>

<script>

import { ref, inject, onMounted } from 'vue'
import * as d3 from 'd3'
import { PCA } from 'ml-pca'

export default {
  name: 'MultidimensionalScalingPlot',
  setup() {

    const ApiService = inject("http")

    const chartContainer = ref(null)

    async function draw() {

      try {

        const response = await ApiService.get("/api/mapperEmbedding")

        console.log(response)

        // const methods = [
        //   {name: 'Method A', results: [[0.1, 0.2, 0.7], [0.15, 0.25, 0.6]]},
        //   {name: 'Method B', results: [[0.8, 0.1, 0.1], [0.75, 0.15, 0.1]]},
        //   {name: 'Method C', results: [[0.3, 0.3, 0.4], [0.35, 0.25, 0.4]]}
        // ]

        const methods = response.data.map(m => ({
          name: m.mapper,
          results: m.embedding
        }))

        methods.sort((a, b) => a.mapper - a.mapper)

        // Compute centroids and PCA
        const centroids = methods.map(m =>
            m.results.reduce((acc, r) => acc.map((v, i) => v + r[i]), [0, 0, 0])
                .map(v => v / m.results.length)
        )

        const pca = new PCA(centroids)
        const transformed = pca.predict(centroids)

        // D3 setup
        const margin = {top: 50, right: 30, bottom: 40, left: 40}
        const width = 600 - margin.left - margin.right
        const height = 400 - margin.top - margin.bottom

        const svg = d3.select(chartContainer.value)
            .append('svg')
            .attr('width', width + margin.left + margin.right)
            .attr('height', height + margin.top + margin.bottom)
            .append('g')
            .attr('transform', `translate(${margin.left},${margin.top})`)

        // Scales
        const xExtent = d3.extent(transformed.data, d => d[0])
        const yExtent = d3.extent(transformed.data, d => d[1])
        const xPadding = (xExtent[1] - xExtent[0]) * 0.1 || 1
        const yPadding = (yExtent[1] - yExtent[0]) * 0.1 || 1
        const xScale = d3.scaleLinear()
            .domain([xExtent[0] - xPadding, xExtent[1] + xPadding])
            .range([0, width])

        const yScale = d3.scaleLinear()
            .domain([yExtent[0] - yPadding, yExtent[1] + yPadding])
            .range([height, 0])

        const colorScale = d3.scaleOrdinal(d3.schemeCategory10)

        // Axes
        svg.append('g')
            .attr('transform', `translate(0,${height})`)
            .call(d3.axisBottom(xScale).tickFormat('').tickSize(0))

        svg.append('g')
            .call(d3.axisLeft(yScale).tickFormat('').tickSize(0))

        // Title
        svg.append('text')
            .attr('x', width / 2)
            .attr('y', -20) // place at the top of the plot area
            .attr('text-anchor', 'middle')
            .style('font-size', '12px')
            .style('font-weight', 'bold')
            .text('Similarity of Mapping Methods in PCA Space')

        // Points
        svg.selectAll('.point')
            .data(transformed.data)
            .enter()
            .append('circle')
            .attr('class', 'point')
            .attr('cx', d => xScale(d[0]))
            .attr('cy', d => yScale(d[1]))
            .attr('r', 8)
            .attr('fill', (d, i) => colorScale(i))
            .attr('stroke', 'white')
            .attr('stroke-width', 2)
            .on('mouseover', function (event, d, i) {
              d3.select(this).attr('r', 12)

              // Tooltip
              const tooltip = d3.select('body').append('div')
                  .attr('class', 'tooltip')
                  .style('opacity', 0)
                  .style('position', 'absolute')
                  .style('background', 'rgba(0,0,0,0.8)')
                  .style('color', 'white')
                  .style('padding', '8px')
                  .style('border-radius', '4px')
                  .style('pointer-events', 'none')

              tooltip.transition().duration(200).style('opacity', .9)
              tooltip.html(methods[i].name)
                  .style('left', (event.pageX + 10) + 'px')
                  .style('top', (event.pageY - 28) + 'px')
            })
            .on('mouseout', function () {
              d3.select(this).attr('r', 8)
              d3.selectAll('.tooltip').remove()
            })

        svg.selectAll('.label')
            .data(transformed.data)
            .enter()
            .append('text')
            .attr('class', 'label')
            .attr('x', d => xScale(d[0]))
            .attr('y', d => yScale(d[1]) - 15)
            .attr('text-anchor', 'middle')
            .style('font-size', '12px')
            .style('font-weight', 'bold')
            .text((d, i) => methods[i].name)

      } catch (error) {
        console.error('Error fetching or processing data:', error)
      }
    }

    onMounted(() => {
      draw()
    })

    return {
      chartContainer
    }
  }
}

</script>

<style scoped>
.pca-chart {
  width: 100%;
  height: 500px;
}
</style>