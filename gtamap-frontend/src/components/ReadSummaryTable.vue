<template>
  <!-- Filtering Options Section -->
  <div class="filter-section tw:mb-4 tw:p-4 tw:border tw:border-gray-300 tw:rounded-lg tw:bg-gray-50">
    <h3 class="tw:text-lg tw:font-medium tw:mb-3">Table Filters</h3>
    <div class="tw:flex tw:items-center tw:gap-4">
      <w-checkbox
        v-model="filterOptions.hideAcceptedRecords"
        label="Hide accepted records"
        @change="applyFilters"
      />
    </div>
  </div>

  <w-table
      :loading="tableData.loading"
      :headers="tableData.headers"
      :items="filteredItems"
      fixed-headers
      :pagination="pagination"
      :selectable-rows="1"
      :selected-rows="tableData.selectedRows"
      :expanded-rows="tableData.expandedRows"
      @row-expand="expandRowClickHandler"
      :sort="tableData.sortKey"
      :sort-function="customSort"
      @update:sort="handleSort"
      class="tw:text-xs"
      expandable-rows>

    <template #item-cell.readDetails="{ item }">
      <w-button @click.stop="openReadDetails(item)" xs>view</w-button>
    </template>

    <template #item-cell.isAccepted="{ item }">
      <w-checkbox
        :model-value="item.isAccepted"
        disabled
      />
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
</template>

<script setup>

import CustomTag from "./CustomTag.vue"

import {defineExpose, inject, nextTick, onMounted, ref, defineEmits, computed} from 'vue'

const emit = defineEmits(['open-read-details', 'content-changed'])

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
    {label: '', key: 'readDetails', sortable: false},
    {label: 'Accepted', key: 'isAccepted', sortable: false},
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

  const items = readSummaryTableData.value.items

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

const readSummaryTableData = ref({items: []})
const readMappingTableData = ref({items: []})

const filterOptions = ref({
  hideAcceptedRecords: false
})

const filteredItems = computed(() => {
  let items = readSummaryTableData.value.items

  if (filterOptions.value.hideAcceptedRecords) {
    items = items.filter(item => !item.isAccepted)
  }

  return items
})

const applyFilters = function() {
  // Filters are applied automatically through the computed property
  // This function can be used for additional filter logic if needed
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

let getReadSummaryTableData = function() {

  tableData.value.loading = true

  ApiService.get(props.url)
      .then(response => {

        if (response.status === 200) {

          readSummaryTableData.value.items = response.data
          pagination.value.total = filteredItems.value.length
          pagination.value.itemsPerPage = 15

          // add parent information to each location
          for (const item of readSummaryTableData.value.items) {
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

let openReadDetails = function(readItem) {
  emit("open-read-details", readItem)
}

const selectAndScrollToRead = async function(readInfo) {

  const targetRead = readSummaryTableData.value.items.find(item => item.qname === readInfo.qname)
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
      if (mapping.mappedBy[i] !== readInfo.mapperName) {
        continue
      }
      if (mapping.readIndices[i] !== readInfo.index) {
        continue
      }
      if (selectedMappingRow !== undefined) {
        console.warn("Multiple mappings found for the selected read, selecting the first one")
        break
      }
      selectedMappingRow = mapping
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

const acceptRecord = function(item) {
  ApiService.post("/api/acceptRecords", {
    recordIds: item.readIndices
  }).then(response => {
      console.log(response)
  }).catch(err => {
    console.error("Error accepting record:", err)
  })
}

const unacceptRecord = function(item) {
  ApiService.post("/api/unacceptRecords", {
    recordIds: item.readIndices
  }).then(response => {
      console.log(response)
  }).catch(err => {
    console.error("Error unaccepting record:", err)
  })
}

const updateItemAcceptance = function(item, isAccepted) {

  for (const o of item.parent.locations) {
    console.log(o)
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

  item.parent.isAccepted = acceptR1 && acceptR2

  emit("content-changed", {})
}

defineExpose({selectAndScrollToRead})

onMounted(() => {
  getReadSummaryTableData()
})
</script>

<style scoped>

</style>
