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

    <div class="tw:flex tw:flex-row">
      <div class="tw:flex-1"></div>
      <MapperMultimappingHeatmap id="mapper-multimap-heatmap"></MapperMultimappingHeatmap>
      <div class="tw:flex-1"></div>
    </div>

    <div class="tw:py-2 tw:flex tw:flex-row">
<!--      <div class="tw:flex-1"></div>-->
<!--      <div id="upset-read-div" class="tw:flex-1"></div>-->
      <div class="tw:flex-1"></div>
      <div id="upset-recordpos-div" class="tw:flex-1"></div>
      <div class="tw:flex-1"></div>
    </div>

    <div class="tw:flex tw:flex-row">

      <div class="tw:flex-1 tw:p-5">
        <div>GTAMap</div>
        <w-table
            :loading="tableData.loading"
            :headers="tableData.headers"
            :items="readSummaryTableData.items"
            fixed-headers
            :pagination="tableData.pagination"
            :selectable-rows="1"
            @row-select="selectedRead = $event"
            style="height: 250px"
            class="tw:text-xs">
          >
        </w-table>
      </div>

      <div class="tw:flex-1 tw:p-5">
        <div>Other Mappers</div>
        <w-table
            :loading="tableData.loading"
            :headers="tableData.headers"
            :items="readSummaryTableData.items"
            fixed-headers
            :pagination="tableData.pagination"
            :selectable-rows="1"
            @row-select="selectedRead = $event"
            style="height: 250px"
            class="tw:text-xs">
          >
        </w-table>
      </div>
    </div>


    <div class="tw:my-5 tw:mb-32">
      <div id="igv-div" class="tw:h-[1000px]"></div>
    </div>

  </div>
</template>

<script>
import {ref, inject} from "vue";

import igv from "../js/igv/dist/igv.esm.js"

import * as UpSetJS from '@upsetjs/bundle'

import * as d3 from 'd3'

import MapperMultimappingHeatmap from "../components/MapperMultimappingHeatmap.vue";

export default {
  name: "OverviewPage",
  components: {
    MapperMultimappingHeatmap
  },
  setup() {

    const ApiService = inject("http")

    let igvInfo = ref({})
    let igvBrowser = ref(undefined)

    const selectedRead = ref({});

    const tableData = ref({
      loading: true,
      headers: [
        {label: 'Read Name', key: 'qname'},
        {label: 'Read Length', key: 'readLength'},
        {label: 'Is Paired?', key: 'isPaired'},
        {label: 'Number Mismatches', key: 'numMismatches'},
        {label: 'Num Gaps', key: 'numGaps'},
        {label: 'Num Insertions', key: 'numInsertions'},
        {label: 'Num Matches GTAMap', key: 'numMatchesGtamap'},
        {label: 'Num Other Mappers Used', key: 'numOtherMappersUsed'}
      ],
      pagination: {
        itemsPerPage: 10,
        total: 0,
        itemsPerPageOptions: []
      }
    })
    const readSummaryTableData = ref({items: []})

    let getReadSummaryTableData = function() {

      tableData.value.loading = true

      ApiService.get("/api/readSummaryTable")
          .then(response => {
            if (response.status === 200) {

              readSummaryTableData.value.items = response.data

              tableData.value.pagination.total = response.data.length

            } else {
              console.error("Failed to fetch read summary table data")
            }

            tableData.value.loading = false

          }).catch(err => {
            console.log(err);
          })
    }

    const getIgvConfig = function() {

      ApiService.get("/api/genomeConfig")
          .then((response) => {
            if (response.status === 200) {

              console.log(response.data)

              igvInfo = ref({
                genome: response.data.genomeConfig,
                tracks: response.data.tracks,
              })

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
        locus: "3:45884425-45903174"
      }

      igv.createBrowser(igvDiv.value, options)
          .then(function (browser) {

            igvBrowser = ref(browser)

            // initialize tracks by track config loaded from API
            for (const trackConfig of igvInfo.value.tracks) {
              browser.loadTrack(trackConfig)
            }
          })
    }

    const upsetDataRead = ref({
      sets: [],
      combinations: []
    });

    const getUpsetDataRead = function() {

      ApiService.get("/api/upsetDataRead")
          .then(response => {
            if (response.status === 200) {

              const elems = response.data;
              const { sets, combinations } = UpSetJS.extractCombinations(elems, { type: 'distinctIntersection' });

              upsetDataRead.value.sets = sets;
              upsetDataRead.value.combinations = combinations;

              renderUpsetRead();

            } else {
              console.error("Failed to fetch upset data");
            }
          }).catch(err => {
            console.log(err);
          })
    }

    let selectionUpsetRead = ref(null);

    function onHoverUpsetRead(set) {
      selectionUpsetRead.value = set;
      renderUpsetRead();
    }

    function renderUpsetRead() {
      let sets = upsetDataRead.value.sets
      let combinations = upsetDataRead.value.combinations

      const dimensions = getChartDimensions(upsetReadContainer);

      console.log(dimensions);

      const props = {
        sets,
        combinations,
        width: 400,
        height: 300,
        selection: selectionUpsetRead.value,
        onHover: onHoverUpsetRead,
        color: "darkorchid",
        selectionColor: "#65cc32",
        //hoverHintColor: "#cc9932",
        //hasSelectionColor: "#cc9932",
        //alternatingBackgroundColor: true,
        //hasSelectionOpacity: 0.5,
        title: "Reads Mapped",
        // barPadding: 2,
        fontSizes: {
          barLabel: "8pt",
          chartLabel: "8pt",
          axisTick: "8pt",
          setLabel: "8pt",
          title: "11pt",
        }
      }
      UpSetJS.render(document.getElementById("upset-read-div"), props);
    }

    const upsetDataRecordPos = ref({
      sets: [],
      combinations: []
    });

    const getUpsetDataRecordPos = function() {

      ApiService.get("/api/upsetDataRecordPos")
          .then(response => {
            if (response.status === 200) {

              const elems = response.data;

              const sets = UpSetJS.extractSets(elems, elem => elem.sets);
              const combinations = UpSetJS.generateCombinations(sets, {
                type: 'distinctIntersection'
              });

              // sort the combinations by cardinality (set size) and name
              const nameToCombination = {};
              for (const combination of combinations) {
                nameToCombination[combination.name] = combination;
              }
              const sortSets = (a, b) => {
                const sizeA = nameToCombination[a.name] ? nameToCombination[a.name].cardinality : 0;
                const sizeB = nameToCombination[b.name] ? nameToCombination[b.name].cardinality : 0;
                const sizeDiff = sizeB - sizeA;
                return sizeDiff !== 0 ? sizeDiff : a.name.localeCompare(b.name);
              };
              const sortCombinations = (a, b) => {
                const sizeDiff = b.cardinality - a.cardinality;
                return sizeDiff !== 0 ? sizeDiff : a.name.localeCompare(b.name);
              };
              sets.sort(sortSets);
              combinations.sort(sortCombinations);
              // end sort by cardinality and name

              upsetDataRecordPos.value.sets = sets;
              upsetDataRecordPos.value.combinations = combinations;

              renderUpsetRecordPos();

            } else {
              console.error("Failed to fetch upset data");
            }
          }).catch(err => {
        console.log(err);
      })
    }

    let selectionUpsetRecordPos = ref(null);

    function onHoverUpsetRecordPos(set) {
      selectionUpsetRecordPos.value = set;
      renderUpsetRecordPos();
    }

    function renderUpsetRecordPos() {
      let sets = upsetDataRecordPos.value.sets
      let combinations = upsetDataRecordPos.value.combinations
      const props = {
        sets,
        combinations,
        width: 800,
        height: 300,
        selection: selectionUpsetRecordPos.value,
        onHover: onHoverUpsetRecordPos,
        color: "darkorchid",
        selectionColor: "#65cc32",
        //hoverHintColor: "#cc9932",
        //hasSelectionColor: "#cc9932",
        //alternatingBackgroundColor: true,
        //hasSelectionOpacity: 0.5,
        title: "Record Position Assignment",
        fontSizes: {
          barLabel: "8pt",
          chartLabel: "8pt",
          axisTick: "8pt",
          setLabel: "8pt",
          title: "11pt",
        }
      }
      UpSetJS.render(document.getElementById("upset-recordpos-div"), props);
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

    return {
      getIgvConfig,
      initializeIgv,
      getUpsetDataRead,
      renderUpsetRead,
      getUpsetDataRecordPos,
      renderUpsetRecordPos,
      selectionUpsetRecordPos,
      getReadSummaryTableData,
      readSummaryTableData,
      tableData,
      selectedRead,
      createChart,
      chartContainer
    }
  },
  mounted() {
    this.getIgvConfig()
    this.getReadSummaryTableData()
    // this.getUpsetDataRead()
    this.getUpsetDataRecordPos()
    this.createChart()
  },
}
</script>

<style scoped>
#mapper-multimap-heatmap {
  transform: scale(1);
  transform-origin: top left;
}
</style>