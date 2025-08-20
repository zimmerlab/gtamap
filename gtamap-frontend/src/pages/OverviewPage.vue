<template>
  <div class="tw:flex tw:flex-col">

    <div class="tw:flex tw:flex-row">

      <div class="tw:flex-1 tw:flex tw:flex-col tw:items-left tw:justify-center tw:px-12
                  tw:text-sm tw:border-r tw:border-[#eee]">
        <div class="tw:text-md tw:font-bold">FASTQ Info:</div>
        <div>Number of initial reads: --</div>
        <div>Average read length: --</div>
        <div>Paired reads?: --</div>
        <div>
            <w-button class="ma1" xs>Open FASTQC Report</w-button>
        </div>
      </div>

      <div class="tw:flex-1 tw:flex tw:flex-col tw:items-left tw:justify-center tw:px-12
                  tw:text-sm tw:border-r tw:border-[#eee]">
        <div class="tw:text-md tw:font-bold">Mapping Info:</div>
        <div>Number of mapped reads: --</div>
        <div>Number of uniquely mapped: --</div>
        <div>Number of multimapped: --</div>
        <div>Number of junctions detected: --</div>
        <div>Coverage: --</div>
        <div>Variants: --</div>
      </div>

      <div class="tw:flex-1 tw:flex tw:flex-col tw:items-left tw:justify-center tw:px-12
                  tw:text-sm tw:border-r tw:border-[#eee]">
        <div class="tw:text-md tw:font-bold">Mapper Consensus:</div>
        <div>Number of mapped reads: --</div>
      </div>
    </div>

<!--    <div class="tw:flex tw:flex-row">-->
<!--      <div class="tw:flex-1"></div>-->
<!--      <MapperMultimappingHeatmap id="mapper-multimap-heatmap" class="tw:flex-1"></MapperMultimappingHeatmap>-->
<!--      <div class="tw:flex-1"></div>-->
<!--    </div>-->

<!--    <MapperMultimappingParallelPlot id="mapper-multimap-parallel"></MapperMultimappingParallelPlot>-->

    <div class="tw:px-12 tw:mb-6">
      <div class="tw:border tw:border-gray-300 tw:rounded tw:p-4">
        <h3 class="tw:text-lg tw:font-bold tw:mb-4 tw:text-center">Mapper Agreement Analysis</h3>
        <div class="tw:flex tw:justify-center">
          <UpsetPlotReadPos title="" url="/api/upsetData"></UpsetPlotReadPos>
        </div>
      </div>
    </div>

<!--    <div class="tw:py-2 tw:flex tw:flex-row">-->
<!--      <div class="tw:flex-1"></div>-->
<!--      <div id="upset-recordposcigar-div" class="tw:flex-1"></div>-->
<!--      <div class="tw:flex-1"></div>-->
<!--    </div>-->

    <!-- <div class="tw:flex tw:flex-row tw:px-12 tw:mb-10"> -->
    <!--   <ReadSummaryTable ref="readSummaryTableRef" class="tw:flex-1" -->
    <!--     url="/api/readSummaryTable" -->
    <!--     @open-read-details="openReadDetailsPage"></ReadSummaryTable> -->
    <!-- </div> -->

    <!-- <div class="tw:my-5 tw:mb-32"> -->
    <!--   <div id="igv-div" class="tw:h-[1000px]"></div> -->
    <!-- </div> -->

    <div class="tw:px-12 tw:mb-6">
      <div class="tw:border tw:border-gray-300 tw:rounded tw:px-4 tw:pb-4 tw:pt-2">
        <h3 class="tw:text-lg tw:font-bold tw:mb-4 tw:text-center tw:text-gray-500">Read Selection</h3>
        <ReadSummaryTable
          ref="readSummaryTableRef"
          class="tw:flex-1"
          url="/api/readSummaryTable"
          @open-read-details="openReadDetailsPage"
          @content-changed="tableUpdated">
        </ReadSummaryTable>
        <Igv 
          ref="igvSummary"
          url="/api/igvConfigTarget"
          @read-click="handleReadClick"
        ></Igv>
      </div>
    </div>

    <div class="tw:px-12 tw:mb-6">
      <div class="tw:border tw:border-gray-300 tw:rounded tw:p-4">
        <h3 class="tw:text-lg tw:font-bold tw:mb-4 tw:text-center">Accepted Reads</h3>
        <Igv
          ref="igvAccepted"
          url="/api/accepted/igvConfig"
          @read-click="handleReadClickAccepted"
        ></Igv>
      </div>
    </div>

    <div>
    </div>

  </div>
</template>

<script>
import {ref, inject} from "vue";
import { useRouter } from "vue-router"

import igv from "../js/igv/dist/igv.esm.js"

import * as UpSetJS from '@upsetjs/bundle'

import * as d3 from 'd3'

import MapperMultimappingHeatmap from "../components/MapperMultimappingHeatmap.vue"
import MapperMultimappingParallelPlot from "../components/MapperMultimappingParallelPlot.vue"
import ReadSummaryTable from "../components/ReadSummaryTable.vue"
import UpsetPlotMapperAgreement from "../components/UpsetPlotMapperAgreement.vue"
import Igv from "../components/Igv.vue"

export default {
  name: "OverviewPage",
  components: {
    UpsetPlotReadPos: UpsetPlotMapperAgreement,
    ReadSummaryTable,
    MapperMultimappingHeatmap,
    MapperMultimappingParallelPlot,
    Igv
  },
  setup() {

    const ApiService = inject("http")
    const router = useRouter()

    let igvSummary = ref(null)
    let igvAccepted = ref(null)

    const readSummaryTableRef = ref(null);

    /**
     * Handle the click event on a read in the igv by sending the read information to the read summary table
     * which brings this read into view and expands it to show all mapping locations.
     * @param readInfo
     */
    const handleReadClick = function(readInfo) {
      console.log("handle read")
      console.log(readInfo)
      readSummaryTableRef.value.selectAndScrollToRead(readInfo)
    }

    const handleReadClickAccepted = function(readInfo) {
      console.log("handle read accepted")
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

    const openReadDetailsPage = function(readItem) {
      const url = router.resolve({
        name: 'readDetailsPage',
        query: {
          readId: readItem.qname
        }
      }).href

      window.open(url, '_blank')
    
      // INFO: this is used to open in the same tab
      // router.push({
      //   name: 'readDetailsPage',
      //   query: {
      //     readId: readItem.qname
      //   }
      // })
    }

    const updateIgvSummary = function() {
      igvSummary.value.update()
    }

    const updateIgvAccepted = function() {
      igvAccepted.value.update()
    }

    const tableUpdated = function() {
      igvAccepted.value.update()
    }

    return {
      igvSummary,
      igvAccepted,
      updateIgvSummary,
      updateIgvAccepted,
      tableUpdated,
      // getIgvConfig,
      // initializeIgv,
      // createChart,
      // chartContainer,
      readSummaryTableRef,
      openReadDetailsPage,
      handleReadClick,
      handleReadClickAccepted
    }
  },
  mounted() {
    // this.getIgvConfig()
    // this.getReadSummaryTableData()
    // this.getUpsetDataRead()
    // this.getUpsetDataRecordPos()
    // this.getUpsetDataRecordPosCigar()
    // this.createChart()
  },
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
