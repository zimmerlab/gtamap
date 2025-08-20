<template>
  <div>
    <div class="filter-section tw:mb-4 tw:p-4 tw:border tw:border-gray-300 tw:rounded-lg tw:bg-gray-50">
    <h3 class="tw:text-lg tw:font-medium tw:mb-3">Table Filters</h3>
    <div class="tw:flex tw:items-center tw:gap-4">
      <w-checkbox
        v-model="filterOptions.hideAcceptedRecords"
        label="Hide accepted records"
      />
      <w-checkbox
        v-model="filterOptions.hideDiscardedRecords"
        label="Hide discarded records"
      />
      <w-checkbox
        v-model="filterOptions.hideNonInProgressRecords"
        label="Show only records in progress"
      />
    </div>
    <div>
      <w-button @click="applyFilters" class="tw:mt-2">
        <span class="tw-text-xs">Apply Filters</span>
      </w-button>
    </div>
  </div>

  <div class="filter-section tw:mb-4 tw:px-4 tw:py-2 tw:border tw:border-gray-300 tw:rounded-lg tw:bg-gray-50">
    <h3 class="tw:text-md tw:font-medium tw:mb-3">Quick Actions</h3>
    <div class="tw:flex tw:items-center tw:gap-4">
      <w-button @click="acceptAllMaxConfidenceReads"><span class="tw-text-xs">Accept All Max Confidence Reads</span></w-button>
    </div>
  </div>

  <div class="filter-section tw:mb-4 tw:px-4 tw:py-2 tw:border tw:border-gray-300 tw:rounded-lg tw:bg-gray-50">
    <h3 class="tw:text-md tw:font-medium tw:mb-3">Selection</h3>
    <div class="tw:flex tw:items-center tw:gap-4">
      <w-button><span class="tw-text-xs">Deselect All</span></w-button>
      <w-button><span class="tw-text-xs">Reset Acceptance for Selected</span></w-button>
      <w-button><span class="tw-text-xs">Accept Selected</span></w-button>
      <w-button><span class="tw-text-xs">Discard Selected</span></w-button>
    </div>
  </div>

  <w-table
      :loading="tableData.loading"
      :headers="tableData.headers"
      :items="summaryTableData.filteredItems"
      :pagination="pagination"
      :selectable-rows="1"
      :selected-rows="tableData.selectedRows"
      :expanded-rows="tableData.expandedRows"
      @row-expand="expandRowClickHandler"
      :sort="tableData.sortKey"
      :sort-function="customSort"
      @update:sort="handleSort"
      class="tw:text-xs"
      fixed-headers
      expandable-rows>

    <template #item-cell.isSelected="{ item }">
      <w-checkbox
        :model-value="item.isSelected"
      />
    </template>

    <template #item-cell.actions="{ item }">
      <w-button @click.stop="openReadDetails(item)" outline xs>
        <w-icon class="mr1">mdi mdi-open-in-new</w-icon>
        view details
      </w-button>
      <w-button v-if="!item.isDiscarded" fg-color="red" class="tw:ml-2" @click.stop="discardRead(item)" outline xs>discard read</w-button>
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
          <w-tag
              v-for="(val, idx) in item.mappedBy"
              :key="idx"
              class="mr4"
              color="primary"
          >{{ val }}</w-tag>
        </span>
      <span v-else>{{ item.mappedBy }}</span>
    </template>

    <template #row-expansion="{ item }">
      <w-table
          :headers="tableData.expandedHeaders"
          :selectable-rows="1"
          :selected-rows="tableData.selectedRowsInExpanded"
          :items="readMappingTableData.items"
      >

        <template #item-cell.isAccepted="{ item }">
          <w-checkbox
            :model-value="item.isAccepted"
            @update:model-value="(value) => updateItemAcceptance(item, value)"
          />
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
          <w-checkbox
            :model-value="isItemSelected(item)"
            @update:model-value="(value) => updateItemAcceptance(item, value)"
          />
        </template>

        <template #item-cell.mappedBy="{ item }">
          <span v-if="Array.isArray(item.mappedBy)">
            <w-tag
              v-for="(val, idx) in Array.from(new Set(item.mappedBy))"
              :key="idx"
              class="mr4"
              color="primary"
            >{{ val }}</w-tag>
          </span>
        <span v-else>{{ item.mappedBy }}</span>
        </template>
      </w-table>
    </template>
  </w-table>
  </div>
</template>

<script setup>

import CustomTag from "./CustomTag.vue"

import {defineExpose, inject, nextTick, onMounted, ref, defineEmits, computed} from 'vue'

const emit = defineEmits(['open-read-details', 'content-changed', 'summary-table-update'])

const ApiService = inject("http")

const props = defineProps({
  url: {
    type: String,
    required: true
  }
})

const pagination = ref({
  total: 0,
  itemsPerPage: 30,
  itemsPerPageOptions: [{value: 15}, {value: 30}],
})

const tableData = ref({
  loading: true,
  sortKey: "+qname",
  headers: [
    {label: 'Select', key: 'isSelected', sortable: false, type: 'boolean'},
    {label: '', key: 'actions', sortable: false},
    {label: 'State', key: 'state', sortable: false},
    {label: 'Confidence', key: 'confidence', sortable: true, type: 'number'},
    {label: 'Read Name', key: 'qname', sortable: true, type: 'string'},
    {label: 'Locations', key: 'numLocations', type: 'number'},
    {label: 'Num Mapped By', key: 'numMappedBy', type: 'number'},
    {label: 'Mapped By', key: 'mappedBy', sortable: false},
  ],
  expandedHeaders: [
    {label: 'Accepted', key: 'isAccepted', sortable: false},
    {label: 'Pair Type', key: 'pairType'},
    {label: 'Contig', key: 'contigName'},
    {label: 'Strand', key: 'isForwardStrand'},
    {label: 'Position', key: 'position'},
    {label: 'Cigar String', key: 'cigar'},
    {label: 'Mismatches', key: 'numMismatches', type: 'number'},
    {label: 'Gaps', key: 'numGaps', type: 'number'},
    {label: 'Num Mapped By', key: 'numMappedBy'},
    {label: 'Mapped By', key: 'mappedBy', sortable: false},
  ],
  selectedRows: [],
  expandedRows: [],
  selectedRowsInExpanded: [],
  customSelectedItems: []
})

const handleSort = function(sortInfo) {
  tableData.value.sortKey = sortInfo[0] ? sortInfo[0] : ""
  customSort()
}

const customSort = function() {

  const items = summaryTableData.value.items

  if (!tableData.value.sortKey || tableData.value.sortKey === "") {
    return items
  }

  const desc = tableData.value.sortKey[0] === '-'
  const key = tableData.value.sortKey.slice(1)

  if (!key) {
    return items
  }

  const sortLikeNumber = function(a, b) {
    return a[key] < b[key] ? -1 : a[key] > b[key] ? 1 : 0
  }

  const sortLikeString = function(a, b) {
    return a[key].localeCompare(b[key])
  }

  const isNumber = tableData.value.headers.some((header) => {
    return header.key === key && header.type === 'number'
  })

  items.sort((a, b) => {
    if (isNumber) {
      return sortLikeNumber(a, b) * (desc ? -1 : 1)
    } else {
      return sortLikeString(a, b) * (desc ? -1 : 1)
    }
  })
  
  return items
}

// content of the outer table (qnames)
const summaryTableData = ref({
  items: [],
  filteredItems: []
})
// content of the inner table (mapping location records per qname)
const readMappingTableData = ref({items: []})

const filterOptions = ref({
  hideAcceptedRecords: false,
  hideDiscardedRecords: false,
  hideNonInProgressRecords: false,
})

const applyFilters = function() {

  let items = summaryTableData.value.items

  if (filterOptions.value.hideAcceptedRecords) {
    items = items.filter(item => !(item.isAcceptedR1 && item.isAcceptedR2))
  }
  if (filterOptions.value.hideDiscardedRecords) {
    items = items.filter(item => !item.isDiscarded)
  }
  if (filterOptions.value.hideNonInProgressRecords) {
    items = items.filter(item => (item.isAcceptedR1 || item.isAcceptedR2) && !(item.isAcceptedR1 && item.isAcceptedR2))
  }

  // emit the table update only once the server was notified about the new filters
  updateSummaryFilters().then(() => {
    emit("summary-table-update")
  })

  deselectRows()

  summaryTableData.value.filteredItems = items
}

const deselectRows = function() {
  tableData.value.selectedRows = []
  tableData.value.selectedRowsInExpanded = []
  tableData.value.expandedRows = []
  readMappingTableData.value.items = []
}

const expandRowClickHandler = function(rowInfo) {

  tableData.value.selectedRows = []
  tableData.value.selectedRowsInExpanded = []

  if (!rowInfo.expanded) {
    tableData.value.expandedRows = []
    readMappingTableData.value.items = []
    return
  }

  expandRow(rowInfo.item)
}

const expandRow = function(rowItem) {
  tableData.value.expandedRows = [rowItem._uid]
  readMappingTableData.value.items = rowItem.locations.slice() || []
}

let getSummaryTableData = function() {

  tableData.value.loading = true

  ApiService.get(props.url)
      .then(response => {

        if (response.status === 200) {

          summaryTableData.value.items = response.data
  
          applyFilters()

          pagination.value.total = summaryTableData.value.filteredItems.length
          pagination.value.itemsPerPage = 15

          // add parent information to each location
          for (const item of summaryTableData.value.items) {
            item.isSelected = false
            for (const l of item.locations) {
              l.parent = item
            }
          }

          customSort()

        } else {
          console.error("Failed to fetch read summary table data")
        }

        tableData.value.loading = false

      }).catch(err => {
    console.log(err);
  })
}

const updateSummaryFilters = function() {
  return ApiService.post("/api/summary/filterUpdate", {
    hideAccepted: filterOptions.value.hideAcceptedRecords,
    hideDiscarded: filterOptions.value.hideDiscardedRecords,
    hideNonInProgress: filterOptions.value.hideNonInProgressRecords,
  })
}

let openReadDetails = function(readItem) {
  emit("open-read-details", readItem)
}

const selectAndScrollToRead = async function(readInfo) {

  const targetRead = summaryTableData.value.items.find(item => item.qname === readInfo.qname)
  if (!targetRead) {
    console.warn("Read not found in the table data")
    return
  }

  pagination.value.page = Math.floor(targetRead._uid / pagination.value.itemsPerPage) + 1

  tableData.value.selectedRows = [targetRead._uid]

  expandRow(targetRead)
  await nextTick()

  let selectedMappingRow = undefined

  for (const mapping of readMappingTableData.value.items) {

    for (let i = 0; i < mapping.mappedBy.length; i++) {
      for (let j = 0; j < readInfo.mapperNames.length; j++) {
        if (mapping.mappedBy[i] === readInfo.mapperNames[j] && mapping.readIndices[i] === readInfo.samIndices[j]) {
          if (selectedMappingRow !== undefined) {
            console.warn("Multiple mappings found for the selected read, selecting the first one")
            break
          }
          selectedMappingRow = mapping
          break
        }
      }
    }
  }

  if (!selectedMappingRow) {
    console.warn("No matching mapping found for the selected read")
    return
  }

  tableData.value.selectedRowsInExpanded = [selectedMappingRow._uid]
}

const isItemSelected = function(item) {
  return tableData.value.customSelectedItems.includes(item._uid)
}

const discardRead = function(item) {
  discardReads([item])
}

const discardReads = function(items) {

  for (const item of items) {
    item.isDiscarded = true
    item.isAcceptedR1 = false
    item.isAcceptedR2 = false
    for (const l of item.locations) {
      l.isAccepted = false
    }
  }

  const qnames = items.map(item => item.qname).flat()
  ApiService.post("/api/summary/discardReads", {
    qnames: qnames
  }).then(response => {
      console.log(response)
  }).catch(err => {
    console.error("Error discarding records:", err)
  })
}

const acceptRecord = function(item) {
  acceptRecords([item])
}

const acceptRecords = function(items) {

  for (const item of items) {
    item.parent.isDiscarded = false
  }

  const recordIds = items.map(item => item.readIndices).flat()
  
  ApiService.post("/api/summary/acceptRecords", {
    recordIds: recordIds
  }).then(response => {
      console.log(response)
  }).catch(err => {
    console.error("Error accepting records:", err)
  })
}

const unacceptRecord = function(item) {
  ApiService.post("/api/summary/unacceptRecords", {
    recordIds: item.readIndices
  }).then(response => {
      console.log(response)
  }).catch(err => {
    console.error("Error unaccepting record:", err)
  })
}

const updateItemAcceptance = function(item, isAccepted) {

  for (const o of item.parent.locations) {
    if (o.pairType !== item.pairType) {
      continue
    }
    if (!o.isAccepted) {
      continue
    }
    o.isAccepted = false
    unacceptRecord(o)
  }
  
  if (isAccepted) {
    item.isAccepted = isAccepted
    acceptRecord(item)
  }

  // check if the parent record should be accepted if both r1 and r2 records are accepted
  let acceptR1 = false
  let acceptR2 = false

  for (const o of item.parent.locations) {
    if (!o.isAccepted) {
      continue
    }
    if (o.pairType === "first") {
      acceptR1 = true
    } else if (o.pairType === "second") {
      acceptR2 = true
    }
    if (acceptR1 && acceptR2) {
      break
    }
  }

  item.parent.isAcceptedR1 = acceptR1
  item.parent.isAcceptedR2 = acceptR2

  emit("content-changed", {})
}

// QUICK ACTIONS

const acceptAllMaxConfidenceReads = function() {

  let items = []

  for (const item of summaryTableData.value.items) {
    if (item.confidenceLevel === 5 && !item.isAccepted) {
      for (const l of item.locations) {
        l.isAccepted = true
        items.push(l)
      }
      item.isAcceptedR1 = true
      item.isAcceptedR2 = true
    }
  }

  acceptRecords(items)

  emit("content-changed", {})
}

// FILTERS

// SELECTION ACTIONS

const deselectAll = function() {
  for (const item of summaryTableData.value.items) {
    item.isSelected = false
  }
}

defineExpose({selectAndScrollToRead})

onMounted(() => {
  getSummaryTableData()
})
</script>

<style scoped>

</style>
