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

    <div class="tw:py-2 tw:flex tw:flex-row">
      <div class="tw:flex-1"></div>
      <UpsetPlotReadPos title="" url="/api/upsetData"></UpsetPlotReadPos>
      <div class="tw:flex-1"></div>
    </div>

<!--    <div class="tw:py-2 tw:flex tw:flex-row">-->
<!--      <div class="tw:flex-1"></div>-->
<!--      <div id="upset-recordposcigar-div" class="tw:flex-1"></div>-->
<!--      <div class="tw:flex-1"></div>-->
<!--    </div>-->

    <div class="tw:flex tw:flex-row tw:px-20 tw:mb-10">
      <ReadSummaryTable ref="readSummaryTableRef" class="tw:flex-1"
        @open-read-details="openReadDetailsPage"></ReadSummaryTable>
    </div>

    <div class="tw:my-5 tw:mb-32">
      <div id="igv-div" class="tw:h-[1000px]"></div>
    </div>

  </div>
</template>

<script>
import {ref, inject} from "vue";
import { useRouter } from "vue-router"

import igv from "../js/igv/dist/igv.esm.js"

import * as UpSetJS from '@upsetjs/bundle'

import * as d3 from 'd3'

import MapperMultimappingHeatmap from "../components/MapperMultimappingHeatmap.vue";
import MapperMultimappingParallelPlot from "../components/MapperMultimappingParallelPlot.vue";
import ReadSummaryTable from "../components/ReadSummaryTable.vue";
import UpsetPlotMapperAgreement from "../components/UpsetPlotMapperAgreement.vue";

export default {
  name: "OverviewPage",
  components: {
    UpsetPlotReadPos: UpsetPlotMapperAgreement,
    ReadSummaryTable,
    MapperMultimappingHeatmap,
    MapperMultimappingParallelPlot
  },
  setup() {

    const ApiService = inject("http")
    const router = useRouter()

    let igvInfo = ref({})
    let igvBrowser = ref(undefined)

    const readSummaryTableRef = ref(null);

    const getIgvConfig = function() {

      ApiService.get("/api/igvConfigTarget")
          .then((response) => {
            if (response.status === 200) {

              igvInfo.value ={
                genome: response.data.genomeConfig,
                tracks: response.data.tracks,
                location: response.data.location,
              }

              initializeIgv()
            }

          }).catch((err) => {
            console.log(err)
          })
    }

    const initializeIgv = function () {

      const igvDiv = ref(document.getElementById("igv-div"))

      // the genome information is loaded from API
      const options = {
        genome: igvInfo.value.genome,
        locus: igvInfo.value.location
      }

      igv.createBrowser(igvDiv.value, options)
        .then(function (browser) {

          igvBrowser = ref(browser)

          // initialize tracks by track config loaded from API
          for (const trackConfig of igvInfo.value.tracks) {
            browser.loadTrack(trackConfig)
          }

          browser.on('trackclick', function (track, popoverData) {

            if (track.type === "alignment" && popoverData && popoverData.length > 0) {

              let readInfo = {
                // the read name (qname) of the read
                qname: "",
                // the name of the mapper which produced this mapping
                mapperName: "",
                // the index of the mapping in the respective sam file
                index: 0,
              }
             
              popoverData.forEach(item => {
                switch (item.name) {
                  case "Read Name":
                    readInfo.qname = item.value;
                    break;
                  case "XM":
                    readInfo.mapperName = item.value;
                    break;
                  case "XI":
                    readInfo.index = item.value;
                    break;
                }
              })
              
              track.alignmentTrack.setHighlightedReads([readInfo.qname], "#0000ff")
              track.updateViews()

              handleReadClick(readInfo)
            }
          });
        })
    }

    /**
     * Handle the click event on a read in the igv by sending the read information to the read summary table
     * which brings this read into view and expands it to show all mapping locations.
     * @param readInfo
     */
    const handleReadClick = function(readInfo) {
      readSummaryTableRef.value.selectAndScrollToRead(readInfo)
    }

    // RIDGE
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
        // Load data
        const data = await d3.csv("https://raw.githubusercontent.com/zonination/perceptions/master/probly.csv")

        console.log(data)

        // Get categories
        const categories = data.columns
        const n = categories.length

        // Add X axis
        const x = d3.scaleLinear()
            .domain([-10, 140])
            .range([0, width])

        svg.append("g")
            .attr("transform", `translate(0, ${height})`)
            .call(d3.axisBottom(x))

        // Create Y scale for densities
        const y = d3.scaleLinear()
            .domain([0, 0.4])
            .range([height, 0])

        // Create Y axis for names
        const yName = d3.scaleBand()
            .domain(categories)
            .range([0, height])
            .paddingInner(1)

        svg.append("g")
            .call(d3.axisLeft(yName))

        // Compute kernel density estimation
        const kde = kernelDensityEstimator(kernelEpanechnikov(7), x.ticks(40))
        const allDensity = []

        for (let i = 0; i < n; i++) {
          const key = categories[i]
          const density = kde(data.map(d => d[key]))
          allDensity.push({ key: key, density: density })
        }

        // Add areas
        svg.selectAll("areas")
            .data(allDensity)
            .join("path")
            .attr("transform", d => `translate(0, ${(yName(d.key) - height)})`)
            .datum(d => d.density)
            .attr("fill", "#69b3a2")
            .attr("stroke", "#000")
            .attr("stroke-width", 1)
            .attr("d", d3.line()
                .curve(d3.curveBasis)
                .x(d => x(d[0]))
                .y(d => y(d[1]))
            )

      } catch (error) {
        console.error('Error loading data:', error)
      }
    }

    // Kernel density estimation functions
    const kernelDensityEstimator = (kernel, X) => {
      return (V) => {
        return X.map(x => {
          return [x, d3.mean(V, v => kernel(x - v))]
        })
      }
    }

    const kernelEpanechnikov = (k) => {
      return (v) => {
        return Math.abs(v /= k) <= 1 ? 0.75 * (1 - v * v) / k : 0
      }
    }

    const openReadDetailsPage = function(readItem) {
      router.push({
        name: 'readDetailsPage',
        query: {
          readId: readItem.qname
        }
      })
    }

    return {
      getIgvConfig,
      initializeIgv,
      createChart,
      chartContainer,
      readSummaryTableRef,
      openReadDetailsPage,
    }
  },
  mounted() {
    this.getIgvConfig()
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
