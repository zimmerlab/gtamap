<template>
  <div v-if="!props.url" class="error-message tw:m-10">
    Error: URL prop is required for IGV component
  </div>
  <div v-else ref="igvDiv"></div>
</template>

<script setup>

import { ref, inject, onMounted, defineExpose, defineEmits, nextTick } from 'vue'
import igv from "../js/igv/dist/igv.esm.js"

const ApiService = inject("http")

const props = defineProps({
  url: {
    type: String,
    required: true
  }
})

const emit = defineEmits(['read-click'])

let igvInfo = ref({})
const igvDiv = ref(null)
let igvBrowser = ref(null)

const getIgvConfig = function() {

  ApiService.get(props.url)
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

  // the genome information is loaded from API
  const options = {
    genome: igvInfo.value.genome,
    locus: igvInfo.value.location
  }

  igv.createBrowser(igvDiv.value, options)
    .then(function (browser) {

      // igvBrowser = ref(browser)
      igvBrowser.value = browser

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
              case "XG":
                readInfo.index = item.value;
                break;
            }
          })
          
          track.alignmentTrack.setHighlightedReads([readInfo.qname], "#0000ff")
          track.updateViews()

          handleReadClick(readInfo)
        }
      })
    })
}

const updateBrowserTracks = function() {
  if (!igvBrowser.value) {
    console.warn("IGV browser not initialized")
    return
  }

  if (!igvDiv.value) {
    console.warn("IGV div not available")
    return
  }

  try {
    igv.removeBrowser(igvBrowser.value)
    igvBrowser.value = null
    
    // Wait for DOM to be ready before reinitializing
    nextTick(() => {
      initializeIgv()
    })
  } catch (error) {
    console.error("Error updating IGV browser tracks:", error)
  }
}

const handleReadClick = function (readInfo) {
  // Emit an event or handle the read click as needed
  console.log("Read clicked:", readInfo)
  emit('read-click', readInfo)
}

const update = function() {
  updateBrowserTracks()
}

defineExpose({
  update
})

onMounted(() => {
  getIgvConfig()
})


</script>

<style scoped>
.error-message {
  color: red;
  font-weight: bold;
  padding: 1rem;
  border: 1px solid red;
  border-radius: 4px;
  background-color: #fee;
}
</style>
