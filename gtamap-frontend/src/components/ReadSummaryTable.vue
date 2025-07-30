<template>
  <w-table
      :loading="tableData.loading"
      :headers="tableData.headers"
      :items="readSummaryTableData.items"
      fixed-headers
      :pagination="tableData.pagination"
      :selectable-rows="1"
      @row-select="selectedRead = $event"
      class="tw:text-xs"
      expandable-rows>

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

<script>

import {ref, inject, onMounted} from 'vue'

export default {
  name: 'ReadSummaryTable',
  setup() {

    const ApiService = inject("http")

    const tableData = ref({
      loading: true,
      headers: [
        {label: '', key: 'readDetails'},
        {label: 'Read Name', key: 'qname'},
        {label: 'Read Length', key: 'readLength'},
        {label: 'Num Locations', key: 'numLocations'},
        {label: 'Num Mapped By', key: 'numMappedBy'},
        {label: 'Mapped By', key: 'mappedBy'},
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
        {label: 'Mapped By', key: 'mappedBy'},
      ],
      pagination: {
        total: 10,
        itemsPerPage: 20,
        // itemsPerPageOptions: [{value: 10, label: '10'}, {value: 20, label: '20'}, {value: 50, label: '50'}],
        itemsPerPageOptions: [{value: 10}, {value: 20}, {value: 50}],
      }
    })

    const selectedRead = ref({})

    const readSummaryTableData = ref({items: []})

    let getReadSummaryTableData = function() {

      tableData.value.loading = true

      ApiService.get("/api/readSummaryTable")
          .then(response => {

            if (response.status === 200) {

              readSummaryTableData.value.items = response.data
              tableData.value.pagination.itemsPerPage = 10
              tableData.value.pagination.total = response.data.length

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

    onMounted(() => {
      getReadSummaryTableData()
    })

    return {
      tableData,
      selectedRead,
      readSummaryTableData,
      getReadSummaryTableData,
      openReadDetails
    }
  }
}
</script>

<style scoped>

</style>