<template>

  <div class="tw:flex tw:flex-col">

    <div class="tw:flex tw:flex-row tw:gap-6 tw:justify-center tw:p-6">
      <div class="tw:w-fit my-card">
        <div class="tw:text-sm tw:font-semibold tw:text-gray-600 tw:mb-2">FASTQ Info</div>
        <div class="tw:text-sm tw:mb-1">Number of initial reads: --</div>
        <div class="tw:text-sm tw:mb-1">Average read length: --</div>
        <div class="tw:text-sm tw:mb-3">Paired reads?: --</div>
        <w-button
          class="ma1"
          xs
          >Open FASTQC Report</w-button
        >
      </div>

      <div class="tw:w-fit my-card">
        <div class="tw:text-sm tw:font-semibold tw:text-gray-600 tw:mb-2">Mapping Info</div>
        <div class="tw:text-sm tw:mb-1">Number of mapped reads: --</div>
        <div class="tw:text-sm tw:mb-1">Number of uniquely mapped: --</div>
        <div class="tw:text-sm tw:mb-1">Number of multimapped: --</div>
        <div class="tw:text-sm tw:mb-1">Number of junctions detected: --</div>
        <div class="tw:text-sm tw:mb-1">Coverage: --</div>
        <div class="tw:text-sm">Variants: --</div>
      </div>

      <div class="tw:w-fit my-card">
        <div class="tw:text-sm tw:font-semibold tw:text-gray-600 tw:mb-2">Mapper Consensus</div>
        <div class="tw:text-sm">Number of mapped reads: --</div>
      </div>
    </div>

    <!--    <div class="tw:flex tw:flex-row">-->
    <!--      <div class="tw:flex-1"></div>-->
    <!--      <MapperMultimappingHeatmap id="mapper-multimap-heatmap" class="tw:flex-1"></MapperMultimappingHeatmap>-->
    <!--      <div class="tw:flex-1"></div>-->
    <!--    </div>-->

    <!--    <MapperMultimappingParallelPlot id="mapper-multimap-parallel"></MapperMultimappingParallelPlot>-->

    <div class="my-card tw:mx-6 tw:my-6">
      <h3 class="tw:text-lg tw:font-bold tw:mb-4 tw:text-center">Mapper Agreement Analysis</h3>
      <div class="tw:flex tw:justify-center">
        <UpsetPlotReadPos
          title=""
          url="/api/upsetData"
        ></UpsetPlotReadPos>
      </div>
    </div>

    <div class="my-card tw:mx-6 tw:my-6">
      <h3 class="tw:text-lg tw:font-bold tw:mb-4 tw:text-center tw:text-gray-500">
        Read Selection
      </h3>
      <ReadSummaryTable
        ref="readSummaryTableRef"
        class="tw:flex-1"
        @open-read-details="openReadDetailsPage"
        @accepted-table-update="tableUpdated"
        @summary-table-update="summaryTableUpdate"
      >
      </ReadSummaryTable>
      <Igv
        ref="igvSummary"
        url="/api/summary/igvConfig"
        :callback="setIgvBrowserSummary"
        @read-click="handleReadClick"
      ></Igv>
    </div>

    <div class="my-card tw:mx-6 tw:my-6">
      <h3 class="tw:text-lg tw:font-bold tw:mb-4 tw:text-center">Accepted Reads</h3>
      <Igv
        ref="igvAccepted"
        url="/api/accepted/igvConfig"
        :callback="setIgvBrowserAccepted"
        @read-click="handleReadClickAccepted"
      ></Igv>
    </div>

    <div></div>
  </div>
</template>

<script>
  import { ref, inject } from 'vue'
  import { useRouter } from 'vue-router'

  import { useDataStore } from '@/store/data.store'

  import * as UpSetJS from '@upsetjs/bundle'

  import * as d3 from 'd3'

  import MapperMultimappingHeatmap from '../components/MapperMultimappingHeatmap.vue'
  import MapperMultimappingParallelPlot from '../components/MapperMultimappingParallelPlot.vue'
  import ReadSummaryTable from '../components/ReadSummaryTable.vue'
  import UpsetPlotMapperAgreement from '../components/UpsetPlotMapperAgreement.vue'
  import Igv from '../components/Igv.vue'

  export default {
    name: 'OverviewPage',
    components: {
      UpsetPlotReadPos: UpsetPlotMapperAgreement,
      ReadSummaryTable,
      MapperMultimappingHeatmap,
      MapperMultimappingParallelPlot,
      Igv,
    },
    setup() {
      const ApiService = inject('http')
      const router = useRouter()

      const dataStore = useDataStore()

      let igvSummary = ref(null)
      let igvAccepted = ref(null)

      const readSummaryTableRef = ref(null)

      /**
       * Handle the click event on a read in the igv by sending the read information to the read summary table
       * which brings this read into view and expands it to show all mapping locations.
       * @param readInfo
       */
      const handleReadClick = function (readInfo) {
        readSummaryTableRef.value.selectAndScrollToRead(readInfo)
      }

      const handleReadClickAccepted = function (readInfo) {
        console.log('handle read accepted')
        console.log(readInfo)
      }

      // RIDGE
      // const chartContainer = ref(null)
      //
      // const createChart = async () => {
      //
      //   // Clear any existing chart
      //   d3.select(chartContainer.value).selectAll("*").remove()
      //
      //   console.log(chartContainer.value)
      //
      //   // Set dimensions and margins
      //   const margin = { top: 60, right: 30, bottom: 20, left: 110 }
      //   const width = 460 - margin.left - margin.right
      //   const height = 400 - margin.top - margin.bottom
      //
      //   // Create SVG
      //   const svg = d3.select(chartContainer.value)
      //       .append("svg")
      //       .attr("width", width + margin.left + margin.right)
      //       .attr("height", height + margin.top + margin.bottom)
      //       .append("g")
      //       .attr("transform", `translate(${margin.left}, ${margin.top})`)
      //
      //   try {
      //     // Load data
      //     const data = await d3.csv("https://raw.githubusercontent.com/zonination/perceptions/master/probly.csv")
      //
      //     console.log(data)
      //
      //     // Get categories
      //     const categories = data.columns
      //     const n = categories.length
      //
      //     // Add X axis
      //     const x = d3.scaleLinear()
      //         .domain([-10, 140])
      //         .range([0, width])
      //
      //     svg.append("g")
      //         .attr("transform", `translate(0, ${height})`)
      //         .call(d3.axisBottom(x))
      //
      //     // Create Y scale for densities
      //     const y = d3.scaleLinear()
      //         .domain([0, 0.4])
      //         .range([height, 0])
      //
      //     // Create Y axis for names
      //     const yName = d3.scaleBand()
      //         .domain(categories)
      //         .range([0, height])
      //         .paddingInner(1)
      //
      //     svg.append("g")
      //         .call(d3.axisLeft(yName))
      //
      //     // Compute kernel density estimation
      //     const kde = kernelDensityEstimator(kernelEpanechnikov(7), x.ticks(40))
      //     const allDensity = []
      //
      //     for (let i = 0; i < n; i++) {
      //       const key = categories[i]
      //       const density = kde(data.map(d => d[key]))
      //       allDensity.push({ key: key, density: density })
      //     }
      //
      //     // Add areas
      //     svg.selectAll("areas")
      //         .data(allDensity)
      //         .join("path")
      //         .attr("transform", d => `translate(0, ${(yName(d.key) - height)})`)
      //         .datum(d => d.density)
      //         .attr("fill", "#69b3a2")
      //         .attr("stroke", "#000")
      //         .attr("stroke-width", 1)
      //         .attr("d", d3.line()
      //             .curve(d3.curveBasis)
      //             .x(d => x(d[0]))
      //             .y(d => y(d[1]))
      //         )
      //
      //   } catch (error) {
      //     console.error('Error loading data:', error)
      //   }
      // }
      //
      // // Kernel density estimation functions
      // const kernelDensityEstimator = (kernel, X) => {
      //   return (V) => {
      //     return X.map(x => {
      //       return [x, d3.mean(V, v => kernel(x - v))]
      //     })
      //   }
      // }
      //
      // const kernelEpanechnikov = (k) => {
      //   return (v) => {
      //     return Math.abs(v /= k) <= 1 ? 0.75 * (1 - v * v) / k : 0
      //   }
      // }

      const openReadDetailsPage = function (readItem) {
        router.push({
          name: 'readDetailsPage',
          query: {
            readId: readItem.qname,
          },
        })
      }

      const updateIgvSummary = function () {
        igvSummary.value.update()
      }

      const updateIgvAccepted = function () {
        igvAccepted.value.update()
      }

      const tableUpdated = function () {
        igvAccepted.value.update()
      }

      const summaryTableUpdate = function () {
        if (!igvSummary.value) {
          console.warn('IGV summary not initialized')
          return
        }
        updateIgvSummary()
      }

      const setIgvBrowserSummary = function(igvDiv, igvBrowser, config) {
        dataStore.setIgvSummary(igvDiv, igvBrowser, config)
      }

      const setIgvBrowserAccepted = function(igvDiv, igvBrowser, config) {
        dataStore.setIgvAccepted(igvDiv, igvBrowser, config)
      }

      return {
        igvSummary,
        igvAccepted,
        updateIgvSummary,
        updateIgvAccepted,
        summaryTableUpdate,
        tableUpdated,
        // createChart,
        // chartContainer,
        readSummaryTableRef,
        openReadDetailsPage,
        handleReadClick,
        handleReadClickAccepted,
        setIgvBrowserSummary,
        setIgvBrowserAccepted,
      }
    },
    mounted() {},
  }
</script>

<style scoped>
  #mapper-multimap-heatmap {
    transform: scale(1);
    transform-origin: top left;
  }

  #mapper-multimap-parallel {
    transform: scale(1);
    transform-origin: top left;
  }
</style>
