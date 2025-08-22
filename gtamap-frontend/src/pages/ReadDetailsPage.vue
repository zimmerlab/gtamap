<template>
  <div v-if="readId === 'state:read-id-loading'">Loading</div>
  <div v-else-if="readId === 'state:no-read-id-found' || !readId || readId === ''">
    No readId found in query parameters
  </div>
  <div v-else>
    <div class="tw:px-6 tw:py-4">
      <w-button @click="goBack"
        class="tw:bg-gray-100 hover:tw:bg-gray-200 tw:text-gray-700 tw:px-4 tw:py-2 tw:rounded tw:border tw:border-gray-300 tw:text-sm tw:font-medium">
        <w-icon class="mr1">mdi mdi-arrow-left</w-icon>
        Go back
      </w-button>
    </div>

    <div class="tw:flex tw:flex-row tw:gap-6 tw:p-6 tw:border-b tw:border-gray-200 tw:mb-6">
      <div class="tw:flex-1 tw:flex tw:flex-col">
        <div class="tw:text-sm tw:font-semibold tw:text-gray-600 tw:mb-2">Read ID</div>
        <div class="tw:text-lg tw:font-bold">{{ readId }}</div>
      </div>

      <div class="tw:flex-1 tw:flex tw:flex-col">
        <div class="tw:text-sm tw:font-semibold tw:text-gray-600 tw:mb-2">Confidence Score</div>
        <div class="tw:text-lg tw:font-bold">0.95</div>
      </div>

      <div class="tw:flex-1 tw:flex tw:flex-col">
        <div class="tw:text-sm tw:font-semibold tw:text-gray-600 tw:mb-2">Tags per Mapper</div>
        <div class="tw:flex tw:flex-wrap tw:gap-2">
          <w-tag color="blue" class="tw:text-xs">Mapper1: Tag1</w-tag>
          <w-tag color="green" class="tw:text-xs">Mapper2: Tag2</w-tag>
          <w-tag color="orange" class="tw:text-xs">Mapper3: Tag3</w-tag>
        </div>
      </div>
    </div>

    <div class="tw:px-6 tw:mb-6">
      <div class="tw:border tw:border-gray-300 tw:rounded tw:px-4 tw:pb-4 tw:pt-2">
        <h3 class="tw:text-lg tw:font-bold tw:mb-4 tw:text-center tw:text-gray-500">
          Read Quality
        </h3>
      </div>
    </div>

    <div class="tw:px-6 tw:mb-6">
      <w-table :loading="tableData.loading" :headers="tableData.headers" :items="tableItems" class="tw:text-xs"
        fixed-headers>
        <template #item-cell.isAccepted="{ item }">
          <w-checkbox :model-value="item.isAccepted"
            @update:model-value="(value) => updateItemAcceptance(item, value)" />
        </template>

        <template #item-cell.targetRegionOverlap="{ item }">
          <w-tag v-if="item.targetRegionOverlap === 0" xs bg-color="red">no overlap</w-tag>
          <w-tag v-else-if="item.targetRegionOverlap === 1" xs bg-color="green">contained</w-tag>
          <w-tag v-else xs bg-color="yellow">partially</w-tag>
        </template>

        <template #item-cell.pairType="{ item }">
          <span v-if="item.pairType === 'first'">R1</span>
          <span v-else-if="item.pairType === 'second'">R2</span>
          <span v-else>NA</span>
        </template>

        <template #item-cell.isForwardStrand="{ item }">
          <span v-if="item.isForwardStrand">+</span>
          <span v-else>-</span>
        </template>

        <template #item-cell.numMismatches="{ item }">
          {{ item.numMismatches !== undefined ? item.numMismatches : 'N/A' }}
        </template>

        <template #item-cell.numGaps="{ item }">
          {{ item.numGaps !== undefined ? item.numGaps : 'N/A' }}
        </template>

        <template #item-cell.mappedBy="{ item }">
          <span v-if="Array.isArray(item.mappedBy)">
            <w-tag v-for="(val, idx) in Array.from(new Set(item.mappedBy))" :key="idx" class="mr4" color="primary">{{
              val }}</w-tag>
          </span>
          <span v-else>{{ item.mappedBy }}</span>
        </template>
      </w-table>
    </div>

    <div class="tw:px-6 tw:mb-6">

      <div v-for="(config, index) in viewerData" :key="config" class="tw:border tw:border-gray-300 tw:rounded-lg tw:mb-4 tw:py-2">

        <h3 class="tw:text-lg tw:font-bold tw:mb-4 tw:text-center tw:text-gray-500">
          Read Group {{ index + 1 }}
        </h3>

        <div  class="tw:px-4">
          <Igv :ref="(el) => setIgvRef(el, index)" </Igv>
        </div>

      </div>
    </div>

  </div>
</template>

<script>
import { ref, onMounted, inject, nextTick } from 'vue'
import { useRoute } from 'vue-router'

import igv from '../js/igv/dist/igv.esm.js'

import Igv from '../components/Igv.vue'

export default {
  name: 'ReadDetailsPage',
  data() {
    const route = useRoute()
    const ApiService = inject('http')

    let readId = ref('state:read-id-loading')

    const getReadIdFromQuery = function () {
      if (route.query.readId) {
        readId.value = route.query.readId
      } else {
        console.error('No readId found in query parameters')
        readId.value = 'state:no-read-id-found'
      }
    }

    // NAVIGATION
    const goBack = function () {
      // Navigate back to the previous page by closing this tab
      window.close()
    }

    // READ TABLE
    let tableItems = ref([])

    const tableData = ref({
      loading: true,
      sortKey: '+qname',
      headers: [
        { label: 'Accepted', key: 'isAccepted', sortable: false },
        {
          label: 'In Target Region?',
          key: 'targetRegionOverlap',
          type: 'number',
          sortable: false,
        },
        { label: 'Pair Type', key: 'pairType' },
        { label: 'Contig', key: 'contigName' },
        { label: 'Strand', key: 'isForwardStrand' },
        { label: 'Position', key: 'position' },
        { label: 'Cigar String', key: 'cigar' },
        { label: 'Mismatches', key: 'numMismatches', type: 'number' },
        { label: 'Gaps', key: 'numGaps', type: 'number' },
        { label: 'Num Mapped By', key: 'numMappedBy' },
        { label: 'Mapped By', key: 'mappedBy', sortable: false },
      ],
    })

    const getReadDetails = function () {
      ApiService.get('/api/details/table', {
        params: {
          qname: readId.value,
        },
      })
        .then((response) => {
          if (response.status === 200) {
            tableItems.value = response.data.locations
            tableData.value.loading = false
          }
        })
        .catch((err) => {
          console.error(err)
          tableData.value.loading = false
        })
    }

    // IGV VIEWER DATA

    let viewerData = ref([])

    const igvRefs = ref([])

    const setIgvRef = (el, index) => {
      if (el) {
        igvRefs.value[index] = el
      }
    }

    const getViewerData = function () {
      ApiService.get('/api/details/viewer/data', {
        params: {
          qname: readId.value,
        }
      })
        .then(async (response) => {
          if (response.status === 200) {
            console.log(response)
            viewerData.value = response.data.viewerConfigs

            await nextTick()

            console.log("after config load")
            console.log(igvRefs.value)

            for (let i = 0; i < igvRefs.value.length; i++) {
              igvRefs.value[i].init(viewerData.value[i])
            }
          }
        })
        .catch((err) => {
          console.error(err)
        })
    }

    // IGV
    let igvBrowser = ref(undefined)

    const getIgvConfig = function () {
      ApiService.get('/api/igvConfigTarget')
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
    }

    const initializeIgv = function () {
      const igvDiv = ref(document.getElementById('igv-div'))

      // the genome information is loaded from API
      const options = {
        genome: igvInfo.value.genome,
        locus: igvInfo.value.location,
      }

      igv.createBrowser(igvDiv.value, options).then(function (browser) {
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
      getReadDetails()
      getViewerData()
    })

    return {
      readId,
      goBack,
      tableItems,
      tableData,
      getReadDetails,
      viewerData,
      getViewerData,
      igvRefs,
      setIgvRef,
    }
  },
  components: {
    Igv
  },
}
</script>

<style scoped></style>
