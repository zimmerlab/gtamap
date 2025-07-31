<template>
  <w-table
      :loading="tableData.loading"
      :headers="tableData.headers"
      :items="readSummaryTableData.items"
      fixed-headers
      :pagination="pagination"
      :selectable-rows="1"
      :selected-rows="tableData.selectedRows"
      v-model:expanded-rows="tableData.expandedRows"
      @row-select="selectedRead = $event"
      :sort="tableData.sortKey"
      :sort-function="customSort"
      class="tw:text-xs"
      expandable-rows>

<!--    <template #item="{ item, index }">-->
<!--      <tr :key="item.qname"-->
<!--          :id="item.qname">-->
<!--&lt;!&ndash;          :ref="el => rowRefs.value[index] = el">&ndash;&gt;-->
<!--        <td>-->
<!--          <w-button @click.stop="openReadDetails(item)" xs>view</w-button>-->
<!--        </td>-->
<!--        <td>{{ item.qname }}</td>-->
<!--        <td>{{ item.readLength }}</td>-->
<!--        <td>{{ item.numLocations }}</td>-->
<!--        <td>{{ item.numMappedBy }}</td>-->
<!--        <td>-->
<!--        <span v-if="Array.isArray(item.mappedBy)">-->
<!--          <w-tag-->
<!--              v-for="(val, idx) in item.mappedBy"-->
<!--              :key="idx"-->
<!--              class="mr4"-->
<!--              color="primary"-->
<!--          >{{ val }}</w-tag>-->
<!--        </span>-->
<!--          <span v-else>{{ item.mappedBy }}</span>-->
<!--        </td>-->
<!--      </tr>-->
<!--    </template>-->

    <template #item-cell.readDetails="{ item }">
      <w-button @click.stop="openReadDetails(item)" xs>view</w-button>
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
          :items="item.locations">

        <template #item-cell.pairType="{ item }">
            <span
                v-if="item.pairType === 'first'"
            >R1</span>
          <span
              v-else-if="item.pairType === 'second'"
          >R2</span>
          <span
              v-else
          >NA</span>
        </template>

        <template #item-cell.isForwardStrand="{ item }">
            <span
                v-if="item.isForwardStrand"
                color="primary"
            >+</span>
          <span
              v-else
              color="primary"
          >-</span>
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
      </w-table>
    </template>
  </w-table>
</template>

<script setup>

import {defineExpose, inject, nextTick, onMounted, ref} from 'vue'

const ApiService = inject("http")

const tableData = ref({
  loading: true,
  sortKey: "+qname",
  headers: [
    {label: '', key: 'readDetails', sortable: false},
    {label: 'Read Name', key: 'qname', sortable: true, type: 'string'},
    {label: 'Read Length', key: 'readLength', type: 'number'},
    {label: 'Num Locations', key: 'numLocations', type: 'number'},
    {label: 'Num Mapped By', key: 'numMappedBy', type: 'number'},
    {label: 'Mapped By', key: 'mappedBy', sortable: false},
  ],
  expandedHeaders: [
    {label: 'Pair Type', key: 'pairType'},
    {label: 'Contig', key: 'contigName'},
    {label: 'Strand', key: 'isForwardStrand'},
    {label: 'Position', key: 'position'},
    {label: 'Cigar String', key: 'cigarString'},
    {label: 'Num Mismatches', key: 'numMismatches'},
    {label: 'Num Gaps', key: 'numGaps'},
    {label: 'Num Mapped By', key: 'numMappedBy'},
    {label: 'Mapped By', key: 'mappedBy', sortable: false},
  ],
  selectedRows: [],
  expandedRows: [],
})

const customSort = function() {

  const items = readSummaryTableData.value.items
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
}

const pagination = ref({
  total: 0,
  itemsPerPage: 30,
  itemsPerPageOptions: [{value: 15}, {value: 30}],
})

const selectedRead = ref({})

const readSummaryTableData = ref({items: []})

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
  console.log("Opening read details for:", readItem)
}

const selectAndScrollToRead = async function(readInfo) {

  const targetRead = readSummaryTableData.value.items.find(item => item.qname === readInfo.qname)
  if (!targetRead) {
    console.warn("Read not found in the table data")
    return
  }

  pagination.value.page = Math.floor(targetRead._uid / pagination.value.itemsPerPage) + 1

  tableData.value.selectedRows = [targetRead._uid]
  tableData.value.expandedRows = [targetRead._uid]
  selectedRead.value = targetRead

  // await nextTick()
}

defineExpose({ selectAndScrollToRead})

onMounted(() => {
  getReadSummaryTableData()
})
</script>

<style scoped>

</style>