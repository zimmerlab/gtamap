<template>
  <div>
    <div class="filter-section tw:mb-3 tw:px-3 tw:py-2 tw:border tw:border-gray-200 tw:rounded tw:bg-gray-25">
      <h4 class="tw:text-sm tw:font-medium tw:mb-2 tw:text-gray-600">Table Filters</h4>
      <div class="tw:flex tw:items-center tw:gap-3 tw:text-sm">
        <w-checkbox v-model="filterHideAcceptedRecords" label="Hide accepted records" xs />
        <w-checkbox v-model="filterHideDiscardedRecords" label="Hide discarded records" xs />
        <w-checkbox v-model="filterHideNonInProgressRecords" label="Show only records in progress" xs />
      </div>
    </div>

    <div v-if="dataStore.getSummaryTableSpecificFilterAsString !== undefined"
      class="filter-section tw:mb-3 tw:px-3 tw:py-2 tw:border tw:border-gray-200 tw:rounded tw:bg-gray-25">
      <h4 class="tw:text-sm tw:font-medium tw:mb-2 tw:text-gray-600">
        Table Filters from Upset Plot
      </h4>
      <div class="tw:flex tw:items-center tw:gap-3 tw:text-sm">
        Showing only Read Mappings for
        <span class="tw:font-bold">{{ dataStore.getSummaryTableSpecificFilterSetName }}</span> with
        filters
        <span class="tw:font-bold">{{ dataStore.getSummaryTableSpecificFilterAsString }}</span>
      </div>
      <w-button outline sm @click="resetSummaryTableSpecificFilters">Reset</w-button>
    </div>

    <div class="filter-section tw:mb-3 tw:px-3 tw:py-2 tw:border tw:border-gray-200 tw:rounded tw:bg-gray-25">
      <h4 class="tw:text-sm tw:font-medium tw:mb-2 tw:text-gray-600">Quick Actions</h4>
      <div class="tw:flex tw:items-center tw:gap-3">
        <w-button @click="acceptAllMaxConfidenceReads" xs><span class="tw:text-xs">Accept All Max Confidence
            Reads</span></w-button>
      </div>
    </div>

    <div class="filter-section tw:mb-3 tw:px-3 tw:py-2 tw:border tw:border-gray-200 tw:rounded tw:bg-gray-25">
      <h4 class="tw:text-sm tw:font-medium tw:mb-2 tw:text-gray-600">Selection</h4>
      <div class="tw:flex tw:items-center tw:gap-2 tw:flex-wrap">
        <w-button @click="deselectAll" xs><span class="tw:text-xs">Unselect All</span></w-button>
        <w-button xs><span class="tw:text-xs">Reset Acceptance for Selected</span></w-button>
        <w-button :disabled="!tableData.actions.accept" @click="acceptSelected" xs><span class="tw:text-xs">Accept
            Selected</span>
        </w-button>
        <w-button xs><span class="tw:text-xs">Discard Selected</span></w-button>
        <w-button @click="resetSelected" xs><span class="tw:text-xs">Reset Selected</span></w-button>
      </div>
    </div>

    <w-table loading="false" :headers="tableData.headers" :items="filteredReads" :pagination="tableData.pagination"
      :selectable-rows="1" :selected-rows="tableData.selectedRows" :expanded-rows="tableData.expandedRows"
      :sort-function="customSortFunction" @row-expand="expandRowClickHandler" class="tw:text-xs" fixed-headers
      expandable-rows>
      <template #item-cell.isSelected="{ item }">
        <w-checkbox :model-value="item.isSelected" @update:model-value="(value) => toggleItemSelection(item, value)" />
      </template>

      <template #item-cell.actions="{ item }">
        <w-button @click.stop="openReadDetails(item)" outline xs>
          <w-icon class="mr1">mdi mdi-open-in-new</w-icon>
          view details
        </w-button>
        <w-button v-if="!item.isDiscarded" fg-color="red" class="tw:ml-2" @click.stop="discardReads([item])" outline
          xs>discard read</w-button>
        <w-button v-else class="tw:ml-2" @click.stop="resetReads([item])" outline xs>reset</w-button>
      </template>

      <template #item-cell.state="{ item }">
        <w-tag v-if="item.isDiscarded" xs bg-color="red">discarded</w-tag>
        <w-tag v-else-if="item.isAcceptedR1 && item.isAcceptedR2" xs bg-color="green">accepted</w-tag>
        <w-tag v-else-if="item.isAcceptedR1 || item.isAcceptedR2" xs bg-color="yellow">in progress</w-tag>
      </template>

      <template #item-cell.confidence="{ item }">
        <CustomTag :level="item.confidenceLevel" size="xs"></CustomTag>
      </template>

      <template #item-cell.mappedBy="{ item }">
        <span v-if="Array.isArray(item.mappedBy)">
          <w-tag v-for="(val, idx) in item.mappedBy" :key="idx" class="mr4" color="primary">{{ val }}</w-tag>
        </span>
        <span v-else>{{ item.mappedBy }}</span>
      </template>

      <template #row-expansion="{ item }">
        <w-table :headers="tableData.expandedHeaders" :selectable-rows="1"
          :selected-rows="tableData.selectedRowsInExpanded" :items="readMappingTableData.items">
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

          <template #item-cell.select="{ item }">
            <w-checkbox :model-value="isItemSelected(item)"
              @update:model-value="(value) => updateItemAcceptance(item, value)" />
          </template>

          <template #item-cell.mappedBy="{ item }">
            <span v-if="Array.isArray(item.mappedBy)">
              <w-tag v-for="(val, idx) in Array.from(new Set(item.mappedBy))" :key="idx" class="mr4" color="primary">{{
                val }}</w-tag>
            </span>
            <span v-else>{{ item.mappedBy }}</span>
          </template>
        </w-table>
      </template>

      <template #pagination="{ range, total }">
        <div class="tw:flex tw:flex-1 tw:flex-row">
          <div class="tw:flex-1"></div>

          <div>
            <w-button :disabled="tableData.pagination.page <= 1" outline sm color="grey" class="tw:pl-0"
              @click="goToPageCustom(tableData.pagination.page - 1)"><w-icon class="tw:mx-0">mdi
                mdi-chevron-left</w-icon>Back</w-button>

            <w-button outline sm color="grey" class="tw:px-2 tw:ml-1" :class="paginationColor(1)"
              @click="goToPageCustom(1)">1</w-button>

            <span v-if="tableData.pagination.showLeft" class="tw:ml-1">...</span>

            <w-button v-for="p in tableData.pagination.middle" :key="p" outline sm color="grey" class="tw:px-2 tw:ml-1"
              :class="paginationColor(p)" @click="goToPageCustom(p)">{{ p }}</w-button>

            <span v-if="tableData.pagination.showRight" class="tw:ml-1">...</span>

            <w-button v-if="tableData.pagination.last > 1" outline sm color="grey" class="tw:px-1 tw:ml-1"
              :class="paginationColor(tableData.pagination.last)" @click="goToPageCustom(tableData.pagination.last)">{{
                tableData.pagination.last }}</w-button>

            <w-button :disabled="tableData.pagination.page >= tableData.pagination.last" outline sm color="grey"
              class="tw:pl-2 tw:pr-0 tw:ml-1" @click="goToPageCustom(tableData.pagination.page + 1)">Next<w-icon
                class="tw:mx-0">mdi mdi-chevron-right</w-icon></w-button>
          </div>

          <div class="tw:flex tw:flex-1 tw:justify-end">
            <div class="tw:ml-1">
              showing {{ range }} of {{ total }} reads (total {{ allReadsCount }})
            </div>
          </div>
        </div>
      </template>
    </w-table>
  </div>
</template>

<script setup>
import CustomTag from './CustomTag.vue'

import { useDataStore } from '../store/data.store'

import { defineExpose, inject, nextTick, onMounted, ref, defineEmits, computed, watch } from 'vue'

const emit = defineEmits([
  'open-read-details',
  'content-changed',
  'summary-table-update',
  'accepted-table-update',
])

const dataStore = useDataStore()

const ApiService = inject('http')
const DataService = inject('data')

let filteredReads = computed(() => dataStore.getSummaryTableReads)
let allReadsCount = computed(() => dataStore.getReads.length)

const tableData = ref({
  loading: true,
  // sortKey: '+qname',
  headers: [
    { label: 'Select', key: 'isSelected', sortable: false, type: 'boolean' },
    { label: '', key: 'actions', sortable: false },
    { label: 'State', key: 'state', sortable: false },
    { label: 'Confidence', key: 'confidence', sortable: true, type: 'number' },
    { label: 'Read Name', key: 'qname', sortable: true, type: 'string' },
    { label: 'Locations', key: 'numLocations', type: 'number' },
    { label: 'Num Mapped By', key: 'numMappedBy', type: 'number' },
    { label: 'Mapped By', key: 'mappedBy', sortable: false },
  ],
  expandedHeaders: [
    { label: 'Accepted', key: 'isAccepted', sortable: false },
    { label: 'In Target Region?', key: 'targetRegionOverlap', type: 'number', sortable: false },
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
  selectedRows: [],
  expandedRows: [],
  selectedRowsInExpanded: [],
  customSelectedItems: [],
  actions: {
    deselect: false,
    accept: true,
    discard: true,
  },
  pagination: {
    total: 0,
    itemsPerPage: 15,
    page: 1,
    pagesCount: 1,
    showLeft: false,
    showRight: false,
    first: 1,
    middle: [],
    last: 1,
  },
})

const updatePagination = function () {
  tableData.value.pagination.total = filteredReads.value.length

  if (tableData.value.pagination.page === 0) {
    tableData.value.pagination.page = 1
  }

  tableData.value.pagination.pagesCount = Math.ceil(
    tableData.value.pagination.total / tableData.value.pagination.itemsPerPage
  )

  tableData.value.pagination.last = tableData.value.pagination.pagesCount

  if (tableData.value.pagination.page > tableData.value.pagination.pagesCount) {
    tableData.value.pagination.page = tableData.value.pagination.pagesCount
  }

  if (tableData.value.pagination.pagesCount <= 6) {
    tableData.value.pagination.middle = []
    for (let i = 2; i < tableData.value.pagination.pagesCount; i++) {
      tableData.value.pagination.middle.push(i)
    }
  } else {
    tableData.value.pagination.middle = []
    let start = Math.max(2, tableData.value.pagination.page - 2)
    let end = Math.min(
      tableData.value.pagination.pagesCount - 1,
      tableData.value.pagination.page + 2
    )
    if (start === 2) {
      end = 5
    } else if (end === tableData.value.pagination.pagesCount - 1) {
      start = tableData.value.pagination.pagesCount - 4
    }
    for (let i = start; i <= end; i++) {
      tableData.value.pagination.middle.push(i)
    }
  }

  tableData.value.pagination.showLeft =
    tableData.value.pagination.middle.length > 0 && tableData.value.pagination.middle[0] > 2
  tableData.value.pagination.showRight =
    tableData.value.pagination.middle.length > 0 &&
    tableData.value.pagination.middle[tableData.value.pagination.middle.length - 1] <
    tableData.value.pagination.pagesCount - 1
}

const customSortFunction = function (sortKeys) {
  resetHighlight()

  let sortKey = sortKeys[0] ? sortKeys[0] : ''

  // change the sortKey to match the data property name for confidence level
  if (sortKey && sortKey.slice(1) === 'confidence') {
    sortKey = sortKey[0] + 'confidenceLevel'
  }

  dataStore.setSummaryTableSortKey(sortKey)
}

// content of the inner table (mapping location records per qname)
const readMappingTableData = ref({ items: [] })

// use filters from dataStore
const filterHideAcceptedRecords = computed({
  get() {
    return dataStore.isSummaryTableHideAcceptedRecords
  },
  set() {
    DataService.toggleSummaryTableFiltersHideAcceptedRecords()
  },
})

const filterHideDiscardedRecords = computed({
  get() {
    return dataStore.isSummaryTableHideDiscardedRecords
  },
  set() {
    DataService.toggleSummaryTableFiltersHideDiscardedRecords()
  },
})

const filterHideNonInProgressRecords = computed({
  get() {
    return dataStore.isSummaryTableHideNonInProgressRecords
  },
  set() {
    DataService.toggleSummaryTableFiltersHideNonInProgressRecords()
  },
})

const deselectRows = function () {
  tableData.value.selectedRows = []
  tableData.value.selectedRowsInExpanded = []
  tableData.value.expandedRows = []
  readMappingTableData.value.items = []
}

const expandRowClickHandler = function (rowInfo) {
  tableData.value.selectedRows = []
  tableData.value.selectedRowsInExpanded = []

  if (!rowInfo.expanded) {
    tableData.value.expandedRows = []
    readMappingTableData.value.items = []
    return
  }

  expandRow(rowInfo.item)
}

const expandRow = function (rowItem) {
  tableData.value.expandedRows = [rowItem._uid]
  readMappingTableData.value.items = rowItem.locations.slice() || []
}

let openReadDetails = function (readItem) {
  emit('open-read-details', readItem)
}

/**
 * Scroll to a record in the summary table, expand its parent qname row and highlight both.
 *
 * @param readInfo - Object containing the read name (qname) and record indices to identify the specific mapping record.
 */
const selectAndScrollToRead = async function (readInfo) {
  resetHighlight()

  const targetReadIndex = filteredReads.value.findIndex((item) => item.qname === readInfo.qname)

  if (!targetReadIndex || targetReadIndex === -1) {
    console.warn('Read not found in the table data')
    return
  }

  const targetPage = Math.floor(targetReadIndex / tableData.value.pagination.itemsPerPage) + 1

  goToPageCustom(targetPage)

  const targetRead = filteredReads.value[targetReadIndex]

  tableData.value.selectedRows = [targetRead._uid]
  expandRow(targetRead)

  await nextTick()

  const targetRecord = targetRead.locations.find((record) =>
    record.readIndices.includes(readInfo.recordIndices[0])
  )

  if (!targetRecord) {
    console.warn('Mapping record not found in the expanded table data')
    return
  }

  tableData.value.selectedRowsInExpanded = [targetRecord._uid]

  console.log('scrolled to record:', targetRecord)
}

const resetHighlight = function () {
  tableData.value.selectedRows = []
  tableData.value.selectedRowsInExpanded = []
  tableData.value.expandedRows = []
  readMappingTableData.value.items = []
}

const isItemSelected = function (item) {
  return tableData.value.customSelectedItems.includes(item._uid)
}

const discardReads = function (reads) {
  DataService.discardReads(reads)
}

const resetReads = function (reads) {
  DataService.resetReads(reads)
}

const acceptReads = function (reads) {
  const records = reads.map((item) => item.locations).flat()
  acceptRecords(records)
}

const acceptRecord = function (record) {
  acceptRecords([record])
}

const acceptRecords = function (records) {
  DataService.acceptRecords(records)
}

const unacceptRecord = function (record) {
  unacceptRecords([record])
}

const unacceptRecords = function (records) {
  for (const r of records) {
    r.isAccepted = false
  }

  emit('accepted-table-update')

  const recordIds = records.map((item) => item.readIndices).flat()

  ApiService.post('/api/summary/unacceptRecords', {
    recordIds: recordIds,
  })
    .then((response) => {
      console.log(response)
    })
    .catch((err) => {
      console.error('Error unaccepting record:', err)
    })
}

const updateItemAcceptance = function (record, isAccepted) {
  DataService.updateRecordAcceptance(record, isAccepted)
}

// QUICK ACTIONS

const acceptAllMaxConfidenceReads = function () {
  DataService.acceptAllMaxConfidenceReads()
}

// FILTERS

// SELECTION ACTIONS
/**
 * Toggle the selection state of an item in the summary table.
 * This function updates the item's selection state and checks if the
 * summary table can accept the current selection.
 * The selection can not be accepted if there are reads where it can not be
 * automatically decided which read to accept (i.e. multiple R1 or R2 mappings).
 * @param {Object} item - The item to toggle selection for.
 * @param {boolean} isSelected - The new selection state.
 */
const toggleItemSelection = function (item, isSelected) {
  // toggle the item selections
  item.isSelected = isSelected

  // check if the selection can be accepted
  let canAccept = true
  let numSelected = 0

  for (const read of summaryTableData.value.filteredItems) {
    if (!read.isSelected) {
      continue
    }

    numSelected++

    let numFirst = read.locations.filter((l) => l.pairType === 'first').length
    let numSecond = read.locations.filter((l) => l.pairType === 'second').length

    if (numFirst === 1 && numSecond === 1) {
      continue
    }

    canAccept = false
    break
  }

  tableData.value.actions.deselect = numSelected > 0
  tableData.value.actions.accept = canAccept
}

const deselectAll = function () {
  for (const item of summaryTableData.value.items) {
    item.isSelected = false
  }
  tableData.value.actions.deselect = false
}

const acceptSelected = function () {
  const selected = summaryTableData.value.filteredItems.filter((item) => item.isSelected)

  acceptReads(selected)
}

const resetSelected = function () {
  const selected = summaryTableData.value.filteredItems.filter((item) => item.isSelected)

  resetReads(selected)
}

const goToPageCustom = function (page) {
  resetHighlight()

  tableData.value.pagination.page = page
  updatePagination()
}

const paginationColor = function (page) {
  if (page === tableData.value.pagination.page) {
    return 'tw:bg-[#234781] tw:text-white tw:border-[#234781]'
  }
}

const resetSummaryTableSpecificFilters = function () {
  DataService.resetSummaryTableSpecificFilters()
}

watch(
  filteredReads,
  (newReads, oldReads) => {
    updatePagination()

    // if the number of reads changed, reset any selection/highlight
    if (newReads && oldReads && newReads.length !== oldReads.length) {
      // emit('summary-table-update', newReads.length)
      resetHighlight()
    }
  },
  { immediate: true }
)

defineExpose({ selectAndScrollToRead })
</script>

<style scoped></style>
