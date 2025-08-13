<template>
  <w-table
      :loading="tableData.loading"
      :headers="tableData.headers"
      :items="readSummaryTableData.items"
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

import {defineExpose, inject, nextTick, onMounted, ref, defineEmits} from 'vue'
import CustomTag from "./CustomTag.vue";

const emit = defineEmits(['open-read-details'])

const ApiService = inject("http")

const pagination = ref({
  total: 0,
  itemsPerPage: 30,
  itemsPerPageOptions: [{value: 15}, {value: 30}],
})

const randomLevel = () => Math.floor(Math.random() * 5) + 1

const tableData = ref({
  loading: true,
  sortKey: "+qname",
  headers: [
    {label: '', key: 'readDetails', sortable: false},
    {label: 'Confidence', key: 'confidence', sortable: true, type: 'number'},
    {label: 'Read Name', key: 'qname', sortable: true, type: 'string'},
    {label: 'Num Locations', key: 'numLocations', type: 'number'},
    {label: 'Num Mapped By', key: 'numMappedBy', type: 'number'},
    {label: 'Mapped By', key: 'mappedBy', sortable: false},
  ],
  expandedHeaders: [
    {label: 'Pair Type', key: 'pairType'},
    {label: 'Contig', key: 'contigName'},
    {label: 'Strand', key: 'isForwardStrand'},
    {label: 'Position', key: 'position'},
    {label: 'Cigar String', key: 'cigar'},
    {label: 'Num Mismatches', key: 'numMismatches', type: 'number'},
    {label: 'Num Gaps', key: 'numGaps', type: 'number'},
    {label: 'Num Mapped By', key: 'numMappedBy'},
    {label: 'Mapped By', key: 'mappedBy', sortable: false},
  ],
  selectedRows: [],
  expandedRows: [],
  selectedRowsInExpanded: []
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

  ApiService.get("/api/readSummaryTable")
      .then(response => {

        if (response.status === 200) {

          readSummaryTableData.value.items = response.data
          pagination.value.total = response.data.length
          pagination.value.itemsPerPage = 15

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

defineExpose({selectAndScrollToRead})

onMounted(() => {
  getReadSummaryTableData()
})
</script>

<style scoped>

</style>
