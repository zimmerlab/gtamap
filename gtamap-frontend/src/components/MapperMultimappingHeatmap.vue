<template>
  <div ref="chartContainer"></div>
</template>

<script>
import { ref, inject } from 'vue'
import * as d3 from 'd3'

export default {
  name: 'MapperMultimappingHeatmap',
  setup() {

    const ApiService = inject("http")

    const chartContainer = ref(null)

    async function createHeatmap() {

      // const margin = { top: 80, right: 25, bottom: 30, left: 40 },
      const margin = { top: 40, right: 25, bottom: 30, left: 120 },
          width = 500 - margin.left - margin.right,
          height = 250 - margin.top - margin.bottom;

      const svg = d3.select(chartContainer.value)
          .append("svg")
          .attr("width", width + margin.left + margin.right)
          .attr("height", height + margin.top + margin.bottom)
          .append("g")
          .attr("transform", `translate(${margin.left}, ${margin.top})`);

      try {

        const response = await ApiService.get("/api/readMultimappingDensityDataHeatmap")
        const data = response.data

        let mappers = Array.from(new Set(data.map(d => d.mapper)))
        mappers.sort((a, b) => {
          if (a === "gtamap") return 1
          if (b === "gtamap") return -1
          return b.localeCompare(a)
        })

        const numMappings = Array.from(new Set(data.map(d => d.numMappings)))

        const x = d3.scaleBand()
            .range([0, width])
            .domain(numMappings)
            .padding(0.05)

        svg.append("g")
            .style("font-size", 11)
            .attr("transform", `translate(0, ${height})`)
            .call(d3.axisBottom(x).tickSize(0))
            .select(".domain").remove()

        const y = d3.scaleBand()
            .range([height, 0])
            .domain(mappers)
            .padding(0.05)

        svg.append("g")
            .style("font-size", 11)
            .call(d3.axisLeft(y).tickSize(0))
            .select(".domain").remove()

        const myColor = d3.scaleSequential()
            .interpolator(d3.interpolateInferno)
            .domain([1, 50]);

        function getTextColor(backgroundColor) {
          const rgb = d3.color(backgroundColor).rgb();
          const luminance = 0.299 * rgb.r + 0.587 * rgb.g + 0.114 * rgb.b;
          return luminance > 140 ? "black" : "white";
        }

        // const tooltip = d3.select(chartContainer.value)
        //     .append("div")
        //     .style("opacity", 0)
        //     .attr("class", "tooltip")
        //     .style("background-color", "white")
        //     .style("border", "solid")
        //     .style("border-width", "2px")
        //     .style("border-radius", "5px")
        //     .style("padding", "5px")
        //     .style("position", "absolute")

        const mouseover = function (event, d) {
          // tooltip.style("opacity", 1)
          d3.select(this).style("stroke", "orange").style("opacity", 1)
        }

        // const mousemove = function (event, d) {
        //   tooltip
        //       .html(`The exact value of<br>this cell is: ${d.count}`)
        //       .style("left", event.pageX + 10 + "px")
        //       .style("top", event.pageY - 28 + "px")
        // }

        const mouseleave = function (event, d) {
          // tooltip.style("opacity", 0)
          d3.select(this).style("stroke", "none").style("opacity", 0.8)
        }

        function handleCellClick(event, d) {
          console.log("Cell clicked:", d)
          // you can add more logic here (e.g., emit an event, navigate, etc.)
        }

        // Vertical grid lines
        svg.selectAll("line.vert-grid")
            .data(numMappings)
            .join("line")
            .attr("class", "vert-grid")
            .attr("x1", d => x(d) + x.bandwidth() / 2)
            .attr("x2", d => x(d) + x.bandwidth() / 2)
            .attr("y1", 0)
            .attr("y2", height)
            .style("stroke", "#e0e0e0")
            .style("stroke-width", 1);
        // Horizontal grid lines
        svg.selectAll("line.horiz-grid")
            .data(mappers)
            .join("line")
            .attr("class", "horiz-grid")
            .attr("x1", 0)
            .attr("x2", width)
            .attr("y1", d => y(d) + y.bandwidth() / 2)
            .attr("y2", d => y(d) + y.bandwidth() / 2)
            .style("stroke", "#e0e0e0")
            .style("stroke-width", 1);

        svg.selectAll("g.cell")
            .data(data, d => d.numMappings + ':' + d.mapper)
            .join("g")
            .attr("class", "cell")
            .on("mouseover", mouseover)
            // .on("mousemove", mousemove)
            .on("mouseleave", mouseleave)
            .on("click", handleCellClick)
            .each(function (d) {

              const cellGroup = d3.select(this)

              if (d.count === 0) {
                return;
              }

              const fillColor = myColor(d.count);
              cellGroup.append("rect")
                  .attr("x", x(d.numMappings))
                  .attr("y", y(d.mapper))
                  .attr("rx", 4)
                  .attr("ry", 4)
                  .attr("width", x.bandwidth())
                  .attr("height", y.bandwidth())
                  .style("fill", fillColor)
                  .style("stroke-width", 1)
                  .style("stroke", "#ccc")
                  .style("opacity", 1)
                  .style("cursor", "pointer")

              cellGroup.append("text")
                  .attr("x", x(d.numMappings) + x.bandwidth() / 2)
                  .attr("y", y(d.mapper) + y.bandwidth() / 2 + 4)
                  .attr("text-anchor", "middle")
                  .style("fill", getTextColor(fillColor))
                  .style("font-size", "9px")
                  .text(d.count)
                  .style("cursor", "pointer")
            })

        // // Add title
        // svg.append("text")
        //     .attr("x", 0)
        //     .attr("y", -50)
        //     .attr("text-anchor", "left")
        //     .style("font-size", "22px")
        //     .text("A d3.js heatmap");

        // Add subtitle
        svg.append("text")
            .attr("x", 0)
            .attr("y", -20)
            .attr("text-anchor", "left")
            .style("font-size", "11px")
            .style("fill", "grey")
            .style("max-width", 400)
            .text("Number of mapping locations per read and mapper");

      } catch (error) {
        console.error("Error fetching data:", error);
      }
    }

    return {
      chartContainer,
      createHeatmap
    }
  },
  mounted() {
    this.createHeatmap()
  },
}
</script>

<style scoped>
.tooltip {
  pointer-events: none;
  position: absolute;
  z-index: 10;
  font-size: 10px;
}
</style>
