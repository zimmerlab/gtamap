<template>
  <div
    v-if="!isValid"
    class="error-message tw:m-10"
  >
    Error: URL prop is required for IGV component
  </div>
  <div
    ref="igvDiv"
    class="tw:z-1"
  ></div>
</template>

<script setup>
  import { ref, inject, onMounted, defineExpose, defineEmits, nextTick } from 'vue'
  import igv from '../js/igv/dist/igv.esm.js'

  const ApiService = inject('http')

  const props = defineProps({
    url: {
      type: String,
      default: undefined,
    },
    config: {
      type: Object,
      default: () => undefined,
    },
    callback: {
      type: Function,
      default: undefined,
    }
  })

  const emit = defineEmits(['read-click'])

  let isValid = ref(false)

  let igvInfo = ref({})
  const igvDiv = ref(null)
  let igvBrowser = ref(null)

  const getIgvConfig = function () {
    if (!props.url && !props.config) {
      console.error('IGV prop url or config is required')
      return
    }
    if (props.url && props.config) {
      console.error('IGV prop url and config cannot be used together')
      return
    }

    isValid.value = true

    if (props.url) {
      ApiService.get(props.url)
        .then((response) => {
          if (response.status === 200) {
            igvInfo.value = {
              genome: response.data.genomeConfig,
              tracks: response.data.tracks,
              location: response.data.location,
            }
            initializeIgv()
          }
        })
        .catch((err) => {
          console.log(err)
        })
    } else if (props.config) {
      // Use the provided config directly
      igvInfo.value = {
        genome: props.config.genome,
        tracks: props.config.tracks,
        location: props.config.location,
      }
      initializeIgv()
    }
  }

  const initializeIgv = function () {
    // the genome information is loaded from API
    const options = {
      genome: igvInfo.value.genome,
      locus: igvInfo.value.location,
    }

    igv.createBrowser(igvDiv.value, options).then(function (browser) {
      // igvBrowser = ref(browser)
      igvBrowser.value = browser

      if (props.callback && typeof props.callback === 'function') {
        props.callback(igvDiv, browser, igvInfo.value)
      }

      // initialize tracks by track config loaded from API
      for (const trackConfig of igvInfo.value.tracks) {
        browser.loadTrack(trackConfig)
      }

      browser.on('trackclick', function (track, popoverData) {
        if (track.type === 'alignment' && popoverData && popoverData.length > 0) {
          let readInfo = {
            // the read name (qname) of the read
            qname: '',
            // the names of the mapper which produced this mapping
            mapperNames: [],
            // the indices of the mappings in the respective sam file
            samIndices: [],
            // list if indices of global records
            recordIndices: [],
          }

          popoverData.forEach((item) => {
            switch (item.name) {
              case 'Read Name':
                readInfo.qname = item.value
                break
              case 'XM':
                readInfo.mapperNames.push(item.value)
                break
              case 'XG':
                readInfo.recordIndices.push(parseInt(item.value))
                break
              case 'XI':
                readInfo.samIndices.push(parseInt(item.value))
                break
            }
          })

          // track.alignmentTrack.setHighlightedReads([readInfo.qname], '#0000ff')
          // track.updateViews()

          handleReadClick(readInfo)
        }
      })
    })
  }

  const updateBrowserTracks = function () {
    if (!igvBrowser.value) {
      console.warn('IGV browser not initialized')
      return
    }

    if (!igvDiv.value) {
      console.warn('IGV div not available')
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
      console.error('Error updating IGV browser tracks:', error)
    }
  }

  const handleReadClick = function (readInfo) {
    // Emit an event or handle the read click as needed
    console.log('Read clicked:', readInfo)
    emit('read-click', readInfo)
  }

  const init = function (config) {

    isValid.value = true

    igvInfo.value = {
      genome: config.genomeConfig,
      tracks: config.tracks,
      location: config.location,
    }

    initializeIgv()
  }

  const update = function () {
    updateBrowserTracks()
  }

  const highlightRead = function(tags, color) {

    for (const track of igvBrowser.value.findTracks()) {
      if (track.type !== 'alignment') {
        continue
      }

      track.setHighlightedTags(tags, color)
      track.updateViews()
    }
  }

  const getIgvBrowser = function() {
    return igvBrowser.value
  }

  defineExpose({
    init,
    update,
    highlightRead,
    getIgvBrowser,
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
