<template>
  <div v-if="readId==='state:read-id-loading'">
    Loading
  </div>
  <div v-else-if="readId==='state:no-read-id-found' || !readId || readId === ''">
    No readId found in query parameters
  </div>
  <div v-else>

    <!-- Top row with 3 columns -->
    <div class="tw:flex tw:flex-row tw:gap-6 tw:p-6 tw:border-b tw:border-gray-200 tw:mb-6">
      
      <!-- Read ID Column -->
      <div class="tw:flex-1 tw:flex tw:flex-col">
        <div class="tw:text-sm tw:font-semibold tw:text-gray-600 tw:mb-2">Read ID</div>
        <div class="tw:text-lg tw:font-bold">{{ readId }}</div>
      </div>

      <!-- Confidence Score Column -->
      <div class="tw:flex-1 tw:flex tw:flex-col">
        <div class="tw:text-sm tw:font-semibold tw:text-gray-600 tw:mb-2">Confidence Score</div>
        <div class="tw:text-lg tw:font-bold">0.95</div>
      </div>

      <!-- Tags per Mapper Column -->
      <div class="tw:flex-1 tw:flex tw:flex-col">
        <div class="tw:text-sm tw:font-semibold tw:text-gray-600 tw:mb-2">Tags per Mapper</div>
        <div class="tw:flex tw:flex-wrap tw:gap-2">
          <w-tag color="blue" class="tw:text-xs">Mapper1: Tag1</w-tag>
          <w-tag color="green" class="tw:text-xs">Mapper2: Tag2</w-tag>
          <w-tag color="orange" class="tw:text-xs">Mapper3: Tag3</w-tag>
        </div>
      </div>
      
    </div>

    <!-- Table Container -->
    <div class="tw:px-6 tw:mb-6">
      <!-- Table placeholder -->
      <div class="tw:border tw:border-gray-300 tw:rounded-lg tw:p-4">
        <ReadSummaryTable
          ref="readSummaryTableRef"
          class="tw:flex-1"
          :url="'/api/readDetails/table?qname='+readId"
          @open-read-details="openReadDetailsPage">
        </ReadSummaryTable>
      </div>
    </div>

    <!-- IGV Viewer Container -->
    <div class="tw:px-6 tw:mb-6">
      <div class="tw:border tw:border-gray-300 tw:rounded-lg">
        <div id="read-details-igv-div" class="tw:h-[600px] tw:bg-white tw:rounded-lg">
          <!-- IGV viewer will be initialized here -->
        </div>
      </div>
    </div>

  </div>
</template>

<script>

import { ref, onMounted } from 'vue'
import { useRoute } from 'vue-router'

import ReadSummaryTable from '../components/ReadSummaryTable.vue'

import igv from "../js/igv/dist/igv.esm.js"

export default {
  name: 'ReadDetailsPage',
  data() {

    const route = useRoute()
    
    let readId = ref("state:read-id-loading")

    const getReadIdFromQuery = function() {
      if (route.query.readId) {
        readId.value = route.query.readId
      } else {
        console.error("No readId found in query parameters")
        readId.value = "state:no-read-id-found"
      }
    }


    // IGV
    let igvBrowser = ref(undefined)

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

          // browser.on('trackclick', function (track, popoverData) {
          //
          //   if (track.type === "alignment" && popoverData && popoverData.length > 0) {
          //
          //     let readInfo = {
          //       // the read name (qname) of the read
          //       qname: "",
          //       // the name of the mapper which produced this mapping
          //       mapperName: "",
          //       // the index of the mapping in the respective sam file
          //       index: 0,
          //     }
          //    
          //     popoverData.forEach(item => {
          //       switch (item.name) {
          //         case "Read Name":
          //           readInfo.qname = item.value;
          //           break;
          //         case "XM":
          //           readInfo.mapperName = item.value;
          //           break;
          //         case "XI":
          //           readInfo.index = item.value;
          //           break;
          //       }
          //     })
          //     
          //     track.alignmentTrack.setHighlightedReads([readInfo.qname], "#0000ff")
          //     track.updateViews()
          //
          //     handleReadClick(readInfo)
          //   }
          // });
        })
    }



    onMounted(() => {
      getReadIdFromQuery()
    })

    return {
      readId
    }
  },
  components: {
    ReadSummaryTable
  },
}

</script>

<style scoped>
</style>
