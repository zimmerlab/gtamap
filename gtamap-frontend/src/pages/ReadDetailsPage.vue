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

    <div class="tw:flex tw:flex-row tw:gap-6 tw:justify-center tw:p-6 tw:mb-6">

      <div class="tw:w-fit tw:border tw:border-gray-300 tw:rounded tw:px-6 tw:py-4 tw:text-center">
        <div class="tw:text-sm tw:font-semibold tw:text-gray-600">Read ID</div>
        <div class="tw:text-lg tw:font-bold">{{ readId }}</div>
      </div>

      <div class="tw:w-fit tw:border tw:border-gray-300 tw:rounded tw:px-6 tw:py-4 tw:text-center">
        <div class="tw:text-sm tw:font-semibold tw:text-gray-600 tw:mb-2">Confidence Score</div>
        <CustomTag :level="detailsData.confidence" size="md"></CustomTag>
      </div>

      <div class="tw:w-fit tw:border tw:border-gray-300 tw:rounded tw:px-6 tw:py-4 tw:text-center">
        <div class="tw:text-sm tw:font-semibold tw:text-gray-600 tw:mb-2">Mapped By</div>

        <!-- R1 Row -->
        <div class="tw:mb-2">
          <div class="tw:flex tw:flex-wrap tw:gap-2 tw:justify-center">
            <w-tag xs bg-color="blue" class="tw:w-8 tw:text-center">R1</w-tag>

            <w-tag v-for="mapper in detailsData.r1Mapped" :key="mapper.mapperName"
              :color="mapper.mapped ? 'green' : 'red'" class="tw:text-xs">
              {{ mapper.mapperName }}
            </w-tag>
          </div>
        </div>

        <!-- R2 Row -->
        <div>
          <div class="tw:flex tw:flex-wrap tw:gap-2 tw:justify-center">
            <w-tag xs bg-color="green" class="tw:w-8 tw:text-center">R2</w-tag>

            <w-tag v-for="mapper in detailsData.r2Mapped" :key="mapper.mapperName"
              :color="mapper.mapped ? 'green' : 'red'" class="tw:text-xs">
              {{ mapper.mapperName }}
            </w-tag>
          </div>
        </div>
      </div>
    </div>

    <div class="tw:px-6 tw:mb-6">
      <div class="tw:border tw:border-gray-300 tw:rounded tw:px-4 tw:pb-4 tw:pt-2">
        <h3 class="tw:text-lg tw:font-bold tw:mb-4 tw:text-center tw:text-gray-500">
          Read Quality
        </h3>
        <div class="tw:flex tw:justify-center">
          <w-tag>coming soon</w-tag>
        </div>
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

    <div class="tw:p-6 tw:mb-6">

      <div v-for="(group, index) in readGroups" :key="group"
        class="tw:border tw:border-gray-300 tw:rounded-lg tw:my-4 tw:py-2">

        <h3 class="tw:text-lg tw:font-bold tw:mb-4 tw:text-center tw:text-gray-500">
          Read Group {{ index + 1 }}
        </h3>

        <div class="tw:mx-4 tw:border tw:rounded-sm tw:border-gray-200 tw:relative">
          <!-- Table Header -->
          <div
            class="tw:bg-gray-50 tw:grid tw:grid-cols-9 tw:gap-2 tw:px-6 tw:py-2 tw:text-xs tw:font-semibold tw:text-gray-600 tw:border-b tw:border-gray-200">
            <div></div>
            <div>Pair Type</div>
            <div>In Target Region?</div>
            <div>Contig</div>
            <div>Strand</div>
            <div>Position</div>
            <div>Cigar String</div>
            <div>Mismatches</div>
            <div>Gaps</div>
          </div>

          <div v-for="byMapper in group.byMapperList" :key="byMapper"
            class="tw:border-b tw:border-gray-200 last:tw:border-b-0">

            <!-- MappedBy Group Header -->
            <div class="tw:bg-blue-50 tw:grid tw:grid-cols-9 tw:gap-2 tw:px-6 tw:border-b tw:border-gray-200">
              <h4 class="tw:text-sm tw:font-semibold tw:text-blue-800">
                {{ byMapper.mapperName }}
              </h4>
              <div></div>
              <div></div>
              <div></div>
              <div></div>
              <div></div>
              <div></div>
              <div></div>
              <div></div>
            </div>

            <!-- R1 Rows -->
            <div v-for="item in byMapper.r1" :key="'r1-' + item.position + item.contigName"
              class="tw:grid tw:grid-cols-9 tw:gap-2 tw:px-6 tw:py-2 tw:text-xs tw:border-b tw:border-gray-100 hover:tw:bg-gray-50 tw:select-none tw:cursor-pointer tw:hover:bg-gray-50"
              :data-contig-cigar="`${item.contigName}-${item.cigar}`" :data-group-id="`group-${index}`">

              <!-- Empty first column -->
              <div></div>

              <!-- Pair Type -->
              <div class="tw:flex tw:items-center">
                <w-tag xs bg-color="blue">R1</w-tag>
              </div>

              <!-- Target Region Overlap -->
              <div class="tw:flex tw:items-center">
                <w-tag v-if="item.targetRegionOverlap === 0" xs bg-color="red">no overlap</w-tag>
                <w-tag v-else-if="item.targetRegionOverlap === 1" xs bg-color="green">contained</w-tag>
                <w-tag v-else xs bg-color="yellow">partially</w-tag>
              </div>

              <!-- Contig -->
              <div class="tw:flex tw:items-center">{{ item.contigName }}</div>

              <!-- Strand -->
              <div class="tw:flex tw:items-center">
                <span v-if="item.isForwardStrand">+</span>
                <span v-else>-</span>
              </div>

              <!-- Position -->
              <div class="tw:flex tw:items-center">{{ item.position }}</div>

              <!-- Cigar String -->
              <div class="tw:flex tw:items-center tw:truncate" :title="item.cigar">{{ item.cigar }}</div>

              <!-- Mismatches -->
              <div class="tw:flex tw:items-center">
                {{ item.numMismatches !== undefined ? item.numMismatches : 'N/A' }}
              </div>

              <!-- Gaps -->
              <div class="tw:flex tw:items-center">
                {{ item.numGaps !== undefined ? item.numGaps : 'N/A' }}
              </div>
            </div>

            <!-- Space between R1 and R2 -->
            <div v-if="byMapper.r1.length > 0 && byMapper.r2.length > 0" class="tw:h-2 tw:bg-gray-100"></div>

            <!-- R2 Rows -->
            <div v-for="item in byMapper.r2" :key="'r2-' + item.position + item.contigName"
              class="tw:grid tw:grid-cols-9 tw:gap-2 tw:px-6 tw:py-2 tw:text-xs tw:border-b tw:border-gray-100 hover:tw:bg-gray-50"
              :data-contig-cigar="`${item.contigName}-${item.cigar}`" :data-group-id="`group-${index}`">

              <!-- Empty first column -->
              <div></div>

              <!-- Pair Type -->
              <div class="tw:flex tw:items-center">
                <w-tag xs bg-color="green">R2</w-tag>
              </div>

              <!-- Target Region Overlap -->
              <div class="tw:flex tw:items-center">
                <w-tag v-if="item.targetRegionOverlap === 0" xs bg-color="red">no overlap</w-tag>
                <w-tag v-else-if="item.targetRegionOverlap === 1" xs bg-color="green">contained</w-tag>
                <w-tag v-else xs bg-color="yellow">partially</w-tag>
              </div>

              <!-- Contig -->
              <div class="tw:flex tw:items-center">{{ item.contigName }}</div>

              <!-- Strand -->
              <div class="tw:flex tw:items-center">
                <span v-if="item.isForwardStrand">+</span>
                <span v-else>-</span>
              </div>

              <!-- Position -->
              <div class="tw:flex tw:items-center">{{ item.position }}</div>

              <!-- Cigar String -->
              <div class="tw:flex tw:items-center tw:truncate" :title="item.cigar">{{ item.cigar }}</div>

              <!-- Mismatches -->
              <div class="tw:flex tw:items-center">
                {{ item.numMismatches !== undefined ? item.numMismatches : 'N/A' }}
              </div>

              <!-- Gaps -->
              <div class="tw:flex tw:items-center">
                {{ item.numGaps !== undefined ? item.numGaps : 'N/A' }}
              </div>
            </div>
          </div>

          <!-- SVG Overlay for Dendrogram Lines -->
          <svg class="tw:absolute tw:top-0 tw:right-0 tw:w-16 tw:h-full tw:pointer-events-none"
            :id="`dendrogram-svg-${index}`">
            <!-- Lines will be drawn here by JavaScript -->
          </svg>
        </div>

        <div class="tw:px-4 tw:mt-6">
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
import CustomTag from '../components/CustomTag.vue'

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

    // READ GROUPS DATA

    let readGroups = ref([])

    let readGroupTableData = ref({
      loading: true,
      sortKey: '+qname',
      headers: [
        // { label: 'Accepted', key: 'isAccepted', sortable: false },
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
        // { label: 'Num Mapped By', key: 'numMappedBy' },
        { label: 'Mapped By', key: 'mappedBy', sortable: false },
      ],
    })

    // IGV VIEWER DATA

    let detailsData = ref({})
    let viewerData = ref([])

    const igvRefs = ref([])

    const setIgvRef = (el, index) => {
      if (el) {
        igvRefs.value[index] = el
      }
    }

    const getViewerData = function () {
      ApiService.get('/api/details/data', {
        params: {
          qname: readId.value,
        }
      })
        .then(async (response) => {
          if (response.status === 200) {

            console.log(response.data)

            detailsData.value.confidence = response.data.confidence
            detailsData.value.r1Mapped = response.data.r1Mapped
            detailsData.value.r2Mapped = response.data.r2Mapped

            for (let group of response.data.readGroups) {

              group.tableData = {
                loading: false
              }

              let byMapperList = []

              for (let location of group.locations) {

                let byMapperIndex = -1

                for (let i = 0; i < byMapperList.length; i++) {
                  if (byMapperList[i].mapperName === location.mappedBy[0]) {
                    byMapperIndex = i
                    break
                  }
                }

                if (byMapperIndex === -1) {
                  byMapperIndex = byMapperList.length
                  byMapperList.push({
                    mapperName: location.mappedBy[0],
                    r1: [],
                    r2: []
                  })
                }

                if (location.pairType === 'first') {
                  byMapperList[byMapperIndex].r1.push(location)
                } else if (location.pairType === 'second') {
                  byMapperList[byMapperIndex].r2.push(location)
                }
              }

              group.byMapperList = byMapperList

              readGroups.value.push(group)
            }

            console.log(readGroups.value)

            await nextTick()

            for (let i = 0; i < igvRefs.value.length; i++) {
              igvRefs.value[i].init(readGroups.value[i].viewerConfig)
            }

            // Draw dendrogram lines after data is loaded and DOM is updated
            await nextTick()
            drawDendrogramLines()

          }
        })
        .catch((err) => {
          console.error(err)
        })
    }

    const drawDendrogramLines = function () {
      readGroups.value.forEach((group, groupIndex) => {
        const svg = document.getElementById(`dendrogram-svg-${groupIndex}`)
        if (!svg) return

        // Clear existing lines
        svg.innerHTML = ''

        // Get all rows with data attributes within this group
        const rows = document.querySelectorAll(`[data-group-id="group-${groupIndex}"]`)

        // Group rows by contig-cigar combination
        const groups = {}
        const rowPositions = []

        rows.forEach((row, index) => {
          const contigCigar = row.getAttribute('data-contig-cigar')
          if (!groups[contigCigar]) {
            groups[contigCigar] = []
          }

          const rect = row.getBoundingClientRect()
          const containerRect = svg.parentElement.getBoundingClientRect()
          const relativeTop = rect.top - containerRect.top + rect.height / 2

          groups[contigCigar].push({
            element: row,
            y: relativeTop,
            index: index
          })

          rowPositions.push({ contigCigar, y: relativeTop, index })
        })

        // Draw connecting lines for groups with multiple rows
        const colors = ['#3b82f6', '#ef4444', '#10b981', '#f59e0b', '#8b5cf6', '#ec4899']
        let colorIndex = 0

        Object.entries(groups).forEach(([contigCigar, groupRows]) => {
          if (groupRows.length > 1) {
            const color = colors[colorIndex % colors.length]
            colorIndex++

            // Sort by Y position
            groupRows.sort((a, b) => a.y - b.y)

            // Draw connecting lines between consecutive rows in the group
            for (let i = 0; i < groupRows.length - 1; i++) {
              const y1 = groupRows[i].y
              const y2 = groupRows[i + 1].y

              // Create a curved path
              const path = document.createElementNS('http://www.w3.org/2000/svg', 'path')
              const d = `M 10 ${y1} Q 30 ${(y1 + y2) / 2} 10 ${y2}`

              path.setAttribute('d', d)
              path.setAttribute('stroke', color)
              path.setAttribute('stroke-width', '2')
              path.setAttribute('fill', 'none')
              path.setAttribute('opacity', '0.7')

              svg.appendChild(path)
            }

            // Add small circles at connection points
            groupRows.forEach(row => {
              const circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle')
              circle.setAttribute('cx', '10')
              circle.setAttribute('cy', row.y)
              circle.setAttribute('r', '3')
              circle.setAttribute('fill', color)
              circle.setAttribute('opacity', '0.8')

              svg.appendChild(circle)
            })
          }
        })
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
      detailsData,
      tableItems,
      tableData,
      readGroups,
      readGroupTableData,
      getReadDetails,
      viewerData,
      getViewerData,
      igvRefs,
      setIgvRef,
      drawDendrogramLines,
    }
  },
  components: {
    Igv,
    CustomTag,
  },
}
</script>

<style scoped></style>
