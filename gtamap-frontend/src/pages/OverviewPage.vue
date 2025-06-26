<template>
  <div class="tw:flex tw:flex-col tw:h-screen">

    <div class="tw:flex tw:flex-row tw:h-1/5">

      <div class="tw:flex-1 tw:flex tw:flex-col tw:items-left tw:justify-center tw:px-12
                  tw:text-sm tw:border-r tw:border-[#eee]">
        <div class="tw:text-md tw:font-bold">FASTQ Info:</div>
        <div>Number of initial reads: --</div>
        <div>Average read length: --</div>
        <div>Paired reads?: --</div>
        <div>RNA or DNA?: --</div>
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

    <div class="tw:h-1/3 tw:py-2 tw:flex tw:flex-row">
      <div class="tw:flex-1"></div>
      <div id="upset-div" class="tw:flex-1"></div>
      <div class="tw:flex-1"></div>
    </div>

    <div class="tw:h-1/3 tw:overflow-auto tw:py-2">
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

    <div class="overview-bottom-section">
      <div id="igv-div"></div>
    </div>

  </div>
</template>

<script>
import {ref, inject} from "vue";

import igv from "../js/igv/dist/igv.esm.js"

import * as UpSetJS from '@upsetjs/bundle'

export default {
  name: "OverviewPage",
  components: {
    // UpSetJS
  },
  data() {

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

    const upsetData = ref({
      sets: [],
      combinations: []
    });

    const getUpsetData = function() {

      ApiService.get("/api/upsetData")
          .then(response => {
            if (response.status === 200) {

              const elems = response.data;
              const { sets, combinations } = UpSetJS.extractCombinations(elems, { type: 'distinctIntersection' });

              console.log("Sets:");
              console.log(sets);

              console.log("Combinations:");
              console.log(combinations);

              upsetData.value.sets = sets;
              upsetData.value.combinations = combinations;

              render();

            } else {
              console.error("Failed to fetch upset data");
            }
          }).catch(err => {
            console.log(err);
          })
    }

    let selection = null;

    function onHover(set) {
      selection = set;
      render();
    }

    function render() {
      let sets = upsetData.value.sets
      let combinations = upsetData.value.combinations
      const props = {
        sets,
        combinations,
        width: 800,
        height: 300,
        selection,
        onHover,
        color: "darkorchid",
        selectionColor: "darkorchid",
        // barPadding: 2,
        fontSizes: {
          barLabel: "8pt",
          chartLabel: "8pt",
          axisTick: "8pt",
          setLabel: "8pt",
        }
      }
      UpSetJS.render(document.getElementById("upset-div"), props);
    }

    return {
      // elems,
      getIgvConfig,
      initializeIgv,
      getUpsetData,
      getReadSummaryTableData,
      readSummaryTableData,
      tableData,
      selectedRead,
      render
    }
  },
  computed: {
    sets() {
      return "Hello World"
      // return extractSets(this.elems)
    }
  },
  mounted() {
    //this.getIgvConfig()
    this.getReadSummaryTableData()
    // this.render()
    this.getUpsetData()
  },
}
</script>

<style>

</style>