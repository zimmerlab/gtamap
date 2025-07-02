<template>
  <div class="w-full h-full flex flex-col items-center justify-center">
    <div id="umap-plot" />
  </div>
</template>

<script>

import { UMAP } from 'umap-js'
import seedrandom from 'seedrandom'

import { ref, inject, onMounted } from 'vue'

export default {
  name: 'UmapTestPlot',
  setup() {

    const ApiService = inject("http")

    const width = 600
    const height = 400
    const margin = 40
    const numPoints = 100
    const dims = 10

    const drawUmap = async () => {

      try {
        const response = await ApiService.get("/api/mapperEmbedding")

        const mapperNames = response.data.map(d => d.mapper)
        const data = response.data.map(d => d.embedding)

        const rng = seedrandom('gtamap');

        // const umap = new UMAP()
        const umap = new UMAP({
          nComponents: 2,
          nEpochs: 400,
          nNeighbors: 3,
          random: rng
        });

        const embedding = umap.fit(data)

        const xExtent = [
          Math.min(...embedding.map((d) => d[0])),
          Math.max(...embedding.map((d) => d[0])),
        ]
        const yExtent = [
          Math.min(...embedding.map((d) => d[1])),
          Math.max(...embedding.map((d) => d[1])),
        ]

        const scaleX = (x) =>
            margin * 2 + ((x - xExtent[0]) / (xExtent[1] - xExtent[0])) * (width - 4 * margin)
        const scaleY = (y) =>
            height - margin * 2 - ((y - yExtent[0]) / (yExtent[1] - yExtent[0])) * (height - 4 * margin)

        const svgNS = 'http://www.w3.org/2000/svg'
        const svg = document.createElementNS(svgNS, 'svg')
        svg.setAttribute('width', width)
        svg.setAttribute('height', height)
        svg.setAttribute('style', 'border:1px solid #ccc; background:#fafafa; font-family:monospace')

        // Title
        const title = document.createElementNS(svgNS, 'text')
        title.setAttribute('x', '50%')
        title.setAttribute('y', '40')
        title.setAttribute('text-anchor', 'middle')
        title.setAttribute('font-size', '16')
        title.setAttribute('font-weight', 'bold')
        title.setAttribute('fill', 'black')
        title.textContent = 'UMAP Projection of Mapper Similarity'
        svg.appendChild(title)

        // X Axis
        const xAxis = document.createElementNS(svgNS, 'line')
        xAxis.setAttribute('x1', margin)
        xAxis.setAttribute('y1', height - margin)
        xAxis.setAttribute('x2', width - margin)
        xAxis.setAttribute('y2', height - margin)
        xAxis.setAttribute('stroke', 'black')
        svg.appendChild(xAxis)

        // Y Axis
        const yAxis = document.createElementNS(svgNS, 'line')
        yAxis.setAttribute('x1', margin)
        yAxis.setAttribute('y1', margin)
        yAxis.setAttribute('x2', margin)
        yAxis.setAttribute('y2', height - margin)
        yAxis.setAttribute('stroke', 'black')
        svg.appendChild(yAxis)

        // Points and labels
        embedding.forEach((point, index) => {
          const [x, y] = [scaleX(point[0]), scaleY(point[1])]

          const circle = document.createElementNS(svgNS, 'circle')
          circle.setAttribute('cx', x)
          circle.setAttribute('cy', y)
          circle.setAttribute('r', 4)
          circle.setAttribute('fill', 'steelblue')
          circle.setAttribute('opacity', 0.7)
          svg.appendChild(circle)

          const label = document.createElementNS(svgNS, 'text')
          label.setAttribute('x', x + 6)
          label.setAttribute('y', y - 6)
          label.setAttribute('font-size', '10')
          label.setAttribute('fill', 'black')
          label.textContent = `${mapperNames[index]}`
          svg.appendChild(label)
        })

        // Mount to placeholder
        const container = document.getElementById('umap-plot')
        container.innerHTML = ''
        container.appendChild(svg)

      } catch (error) {
        console.error("Error fetching data:", error)
      }
    }

    onMounted(() => {
      drawUmap()
    })

    return {
      drawUmap
    }
  }
}
</script>

<style scoped>

</style>