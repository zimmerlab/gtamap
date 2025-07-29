<!--<template>-->
<!--    <div ref="chartContainer" id="my_dataviz"></div>-->
<!--</template>-->

<!--<script>-->

<!--import {inject, ref} from 'vue'-->

<!--import * as d3 from 'd3'-->

<!--export default {-->
<!--  name: 'MapperMultimappingRidgePlot',-->
<!--  setup() {-->

<!--    const ApiService = inject("http")-->

<!--    const chartContainer = ref(null)-->

<!--    const createChart = async () => {-->

<!--      // Clear any existing chart-->
<!--      d3.select(chartContainer.value).selectAll("*").remove()-->

<!--      console.log(chartContainer.value)-->

<!--      // Set dimensions and margins-->
<!--      const margin = { top: 60, right: 30, bottom: 20, left: 110 }-->
<!--      const width = 460 - margin.left - margin.right-->
<!--      const height = 400 - margin.top - margin.bottom-->

<!--      // Create SVG-->
<!--      const svg = d3.select(chartContainer.value)-->
<!--          .append("svg")-->
<!--          .attr("width", width + margin.left + margin.right)-->
<!--          .attr("height", height + margin.top + margin.bottom)-->
<!--          .append("g")-->
<!--          .attr("transform", `translate(${margin.left}, ${margin.top})`)-->

<!--      try {-->

<!--        const response = await ApiService.get("/api/readMultimappingDensityData")-->

<!--        console.log(response)-->

<!--        const data = response.data.data-->
<!--        // sort data by mapper-->
<!--        data.sort((a, b) => a.mapper.localeCompare(b.mapper))-->

<!--        const categories = data.map(e => e.mapper)-->

<!--        console.log(categories)-->

<!--        // const x = d3.scaleLog()-->
<!--        //     .domain([1, 140])-->
<!--        //     .range([0, width])-->

<!--        const x = d3.scaleLinear()-->
<!--            .domain([0.5, 10])-->
<!--            .range([0, width])-->

<!--        svg.append("g")-->
<!--            .attr("transform", `translate(0, ${height})`)-->
<!--            .call(d3.axisBottom(x))-->

<!--        // const y = d3.scaleLinear()-->
<!--        //     .domain([0, 0.4])-->
<!--        //     .range([height, 0])-->

<!--        // const y = d3.scaleLinear()-->
<!--        //     .domain([0, 40])-->
<!--        //     .range([height, 0])-->

<!--        const y = d3.scaleSymlog()-->
<!--            .domain([0, 2000])-->
<!--            .range([height, 0])-->

<!--        const yName = d3.scaleBand()-->
<!--            .domain(categories)-->
<!--            .range([0, height])-->
<!--            .paddingInner(1)-->

<!--        svg.append("g")-->
<!--            .call(d3.axisLeft(yName))-->

<!--        const kde = kernelDensityEstimator(kernelEpanechnikov(0.3), x.ticks(200))-->
<!--        // const kde = kernelDensityEstimator(kernelEpanechnikov(70), x.ticks(100))-->
<!--        const allDensity = []-->

<!--        for (let i = 0; i < data.length; i++) {-->
<!--          const key = categories[i]-->
<!--          const density = kde(data[i].readAssignments)-->
<!--          allDensity.push({ key: key, density: density })-->
<!--        }-->

<!--        // Add areas-->
<!--        svg.selectAll("areas")-->
<!--            .data(allDensity)-->
<!--            .join("path")-->
<!--            .attr("transform", d => `translate(0, ${(yName(d.key) - height)})`)-->
<!--            .datum(d => d.density)-->
<!--            .attr("fill", "#69b3a2")-->
<!--            .attr("stroke", "#000")-->
<!--            .attr("stroke-width", 1)-->
<!--            .attr("d", d3.line()-->
<!--                .curve(d3.curveBasis)-->
<!--                .x(d => x(d[0]))-->
<!--                .y(d => y(d[1]))-->
<!--            )-->

<!--      } catch (error) {-->
<!--        console.error('Error loading data:', error)-->
<!--      }-->
<!--    }-->

<!--    // Kernel density estimation functions-->
<!--    const kernelDensityEstimator = (kernel, X) => {-->
<!--      return (V) => {-->
<!--        return X.map(x => {-->
<!--          return [x, d3.mean(V, v => kernel(x - v))]-->
<!--        })-->
<!--      }-->
<!--    }-->

<!--    const kernelEpanechnikov = (k) => {-->
<!--      return (v) => {-->
<!--        return Math.abs(v /= k) <= 1 ? 0.75 * (1 - v * v) / k : 0-->
<!--      }-->
<!--    }-->

<!--    return {-->
<!--      chartContainer,-->
<!--      createChart-->
<!--    }-->
<!--  },-->
<!--  mounted() {-->
<!--    this.createChart()-->
<!--  }-->
<!--}-->

<!--</script>-->

<!--<style scoped>-->

<!--</style>-->


<template>
  <div ref="chartContainer" id="my_dataviz"></div>
</template>
<script>
import {inject, ref} from 'vue'
import * as d3 from 'd3'
export default {
  name: 'MapperMultimappingRidgePlot',
  setup() {
    const ApiService = inject("http")
    const chartContainer = ref(null)
    const createChart = async () => {
      // Clear any existing chart
      d3.select(chartContainer.value).selectAll("*").remove()
      console.log(chartContainer.value)
      // Set dimensions and margins
      const margin = { top: 60, right: 30, bottom: 20, left: 110 }
      const width = 460 - margin.left - margin.right
      const height = 400 - margin.top - margin.bottom
      // Create SVG
      const svg = d3.select(chartContainer.value)
          .append("svg")
          .attr("width", width + margin.left + margin.right)
          .attr("height", height + margin.top + margin.bottom)
          .append("g")
          .attr("transform", `translate(${margin.left}, ${margin.top})`)
      try {
        const response = await ApiService.get("/api/readMultimappingDensityData")
        console.log(response)
        const data = response.data.data
        // sort data by mapper
        data.sort((a, b) => a.mapper.localeCompare(b.mapper))
        const categories = data.map(e => e.mapper)
        console.log(categories)

        const x = d3.scaleLinear()
            .domain([0.5, 10])
            .range([0, width])

        svg.append("g")
            .attr("transform", `translate(0, ${height})`)
            .call(d3.axisBottom(x))

        const y = d3.scaleSymlog()
            .domain([0, 2000])
            .range([height, 0])

        const yName = d3.scaleBand()
            .domain(categories)
            .range([0, height])
            .paddingInner(1)

        svg.append("g")
            .call(d3.axisLeft(yName))

        // Create histogram for each category
        const histogram = d3.histogram()
            .domain(x.domain())
            .thresholds(x.ticks(40)) // Adjust number of bins

        const allHistograms = []
        for (let i = 0; i < data.length; i++) {
          const key = categories[i]
          const bins = histogram(data[i].readAssignments)
          // Normalize by total count to get density
          const total = data[i].readAssignments.length
          const normalizedBins = bins.map(bin => ({
            x0: bin.x0,
            x1: bin.x1,
            length: bin.length / total // Convert to density
          }))
          allHistograms.push({ key: key, bins: normalizedBins })
        }

        // Add histogram bars for each category
        svg.selectAll(".histogram-group")
            .data(allHistograms)
            .enter()
            .append("g")
            .attr("class", "histogram-group")
            .attr("transform", d => `translate(0, ${yName(d.key) - height})`)
            .selectAll("rect")
            .data(d => d.bins)
            .enter()
            .append("rect")
            .attr("x", d => x(d.x0))
            .attr("width", d => Math.max(0, x(d.x1) - x(d.x0) - 1))
            .attr("y", d => y(d.length))
            .attr("height", d => y(0) - y(d.length))
            .attr("fill", "#69b3a2")
            .attr("stroke", "#000")
            .attr("stroke-width", 0.5)

      } catch (error) {
        console.error('Error loading data:', error)
      }
    }

    return {
      chartContainer,
      createChart
    }
  },
  mounted() {
    this.createChart()
  }
}
</script>
<style scoped>
</style>