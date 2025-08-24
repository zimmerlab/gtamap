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

      <div class="tw:w-fit my-card">
        <div class="tw:text-sm tw:font-semibold tw:text-gray-600">Read ID</div>
        <div class="tw:text-lg tw:font-bold">{{ readId }}</div>
        <w-button bg-color="red" sm>discard</w-button>
      </div>

      <div class="tw:w-fit my-card">
        <div class="tw:text-sm tw:font-semibold tw:text-gray-600 tw:mb-2">Confidence Score</div>
        <CustomTag :level="detailsData.confidence" size="md"></CustomTag>
      </div>

      <div class="tw:w-fit my-card">
        <div class="tw:text-sm tw:font-semibold tw:text-gray-600 tw:mb-2">Genomic Locations</div>
        <div v-for="group in readGroups" :key="group" class="tw:mb-1">
          <w-tag>chr{{ group.interval.contig }} {{ group.interval.start }} -
            {{ group.interval.end }}</w-tag>
        </div>
      </div>


      <div class="tw:w-fit my-card">
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

    <div class="my-card tw:mx-6 tw:mb-6">
      <h3 class="tw:text-lg tw:font-bold tw:mb-4 tw:text-center tw:text-gray-500">
        Read Quality
      </h3>
      <div class="tw:flex tw:justify-center">
        <w-tag>coming soon</w-tag>
      </div>
    </div>

    <div class="my-card tw:mx-6 tw:mb-6">

      <h3 class="tw:text-lg tw:font-bold tw:mb-4 tw:text-center tw:text-gray-500">
        Read Overview
      </h3>

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

    <div v-for="(group, index) in readGroups" :key="group" class="my-card tw:mx-6 tw:mb-6">

      <h3 class="tw:mb-4 tw:text-center">
        <div class="tw:text-lg tw:font-bold tw:text-gray-500">
          Genomic Location
        </div>
        <w-tag>
          chr{{ group.interval.contig }} {{ group.interval.start }} -
          {{ group.interval.end }}
        </w-tag>
      </h3>

      <div class="tw:border tw:rounded-sm tw:border-gray-200 tw:relative">
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
          <div v-for="item in byMapper.r1" :key="'r1-' + item.position + item.cigar"
            class="tw:grid tw:grid-cols-9 tw:gap-2 tw:px-6 tw:py-2 tw:text-xs tw:border-b tw:border-gray-100 tw:select-none tw:cursor-pointer"
            :data-contig-cigar="`${item.contigName}-${item.position}-${item.cigar}`" :data-group-id="`group-${index}`"
            :data-item-index="`${item.index}`"
            @mouseenter="highlightEquivalentRows(`${item.contigName}-${item.position}-${item.cigar}`, item, group, index)"
            @mouseleave="clearHighlight">

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
          <div v-for="item in byMapper.r2" :key="'r2-' + item.position + item.cigar"
            class="tw:grid tw:grid-cols-9 tw:gap-2 tw:px-6 tw:py-2 tw:text-xs tw:border-b tw:border-gray-100 tw:select-none tw:cursor-pointer"
            :data-contig-cigar="`${item.contigName}-${item.position}-${item.cigar}`" :data-group-id="`group-${index}`"
            :data-item-index="`${item.index}`"
            @mouseenter="highlightEquivalentRows(`${item.contigName}-${item.position}-${item.cigar}`, item, group, index)"
            @mouseleave="clearHighlight">

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

      <div class="tw:mt-6">
        <Igv :ref="(el) => setIgvRef(el, index)" </Igv>
      </div>

    </div>

  </div>
</template>

<script>
import { ref, onMounted, inject, nextTick } from 'vue'
import { useRoute, useRouter } from 'vue-router'

import igv from '../js/igv/dist/igv.esm.js'

import Igv from '../components/Igv.vue'
import CustomTag from '../components/CustomTag.vue'

export default {
  name: 'ReadDetailsPage',
  data() {
    const route = useRoute()
    const router = useRouter()
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
      router.push({ name: 'index' })
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

            detailsData.value.confidence = response.data.confidence
            detailsData.value.r1Mapped = response.data.r1Mapped
            detailsData.value.r2Mapped = response.data.r2Mapped

            for (let group of response.data.readGroups) {

              group.tableData = {
                loading: false
              }

              let byMapperList = []

              for (let location of group.locations) {

                location.index = location.readIndices[0]

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

              // sort byMapperList by mapperName (gtamap always first, then alphabetically)
              byMapperList.sort((a, b) => {
                if (a.mapperName === 'gtamap') return -1
                if (b.mapperName === 'gtamap') return 1
                return a.mapperName.localeCompare(b.mapperName)
              })

              group.byMapperList = byMapperList

              readGroups.value.push(group)
            }

            // need to wait for the next tick to be able to access the igvRefs
            await nextTick()

            // update each igv with the respective read group data
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

    const colors = ['#3b82f6', '#ef4444', '#10b981', '#f59e0b', '#8b5cf6', '#ec4899']

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

    onMounted(() => {
      getReadIdFromQuery()
      getReadDetails()
      getViewerData()
    })

    const highlightEquivalentRows = function (contigCigar, item, group, groupIndex) {

      // find all rows with matching data-contig-cigar attribute
      const allRows = document.querySelectorAll(`[data-contig-cigar="${contigCigar}"]`)


      // all unique contig-cigar combinations in the current group
      const groupElement = allRows[0]?.closest('[data-group-id]')

      if (groupElement) {

        const groupId = groupElement.getAttribute('data-group-id')
        const groupRows = document.querySelectorAll(`[data-group-id="${groupId}"]`)

        // group rows by their data-contig-cigar attribute and count them
        let groups = []
        for (let row of groupRows) {
          const key = row.getAttribute('data-contig-cigar')
          let index = -1
          for (let i = 0; i < groups.length; i++) {
            if (groups[i].key === key) {
              groups[i].count += 1
              index = i
              break
            }
          }
          if (index == -1) {
            groups.push({ key: key, count: 1 })
          }
        }

        // find the color index based on rows with multiple entries
        let index = -1
        // check if the current row has multiple entries and deserves a color
        let mult = false
        for (const g of groups) {
          if (g.count > 1) {
            index += 1
          }
          if (g.key === contigCigar) {
            mult = g.count > 1
            break
          }
        }

        // apply default color (gray) if no equivalent rows found
        let color = "#d3d3d3"
        let readColor = "#FD7E14"
        // use the same color as in the SVG
        if (mult && index < colors.length) {
          color = colors[index]
          readColor = colors[index]
        }

        // apply the same color as background with opacity
        allRows.forEach(row => {
          row.style.backgroundColor = color + '33'
          row.classList.add('highlighted-row')
        })

        const readIndices = Array.from(allRows).map(row => "XG:" + row.getAttribute('data-item-index'))
        igvRefs.value[groupIndex].highlightRead(readIndices, readColor)
      }
    }

    const clearHighlight = function () {

      // clear highlighted igv reads
      for (const igvRef of igvRefs.value) {
        igvRef.highlightRead([], '#000000')
      }

      // remove highlight class from all previously highlighted rows
      const highlightedRows = document.querySelectorAll('.highlighted-row')
      highlightedRows.forEach(row => {
        row.classList.remove('highlighted-row')
        row.style.backgroundColor = ''
      })
    }

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
      highlightEquivalentRows,
      clearHighlight,
    }
  },
  components: {
    Igv,
    CustomTag,
  },
}
</script>

<style scoped>
.highlighted-row {
  /* @apply tw:bg-blue-100 tw:border-blue-300; */
  background-color: blue;
}
</style>
