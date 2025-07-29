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
      <template #row-expansion="{ item }">
        <w-table
            :headers="tableData.expandedHeaders"
            items=""
        ></w-table>
<!--        <w-tag class="mr4" color="primary">tag</w-tag>-->
      </template>
  </w-table>
</template>

<script>

import {ref, inject, onMounted} from 'vue'

export default {
  name: 'TablePage',
  setup() {

    const ApiService = inject("http")

    const tableData = ref({
      loading: true,
      headers: [
        {label: 'Read Name', key: 'qname'},
        {label: 'Read Length', key: 'readLength'},
        {label: 'Is Paired?', key: 'isPaired'},
        {label: 'Number Mismatches', key: 'numMismatches'},
        {label: 'Num Gaps', key: 'numGaps'},
        {label: 'Num Insertions', key: 'numInsertions'},
        {label: 'Num Matches GTAMap', key: 'numMatchesGtamap'},
        {label: 'Num Other Mappers Used', key: 'numOtherMappersUsed'}
      ],
      expandedHeaders: [
        {label: 'Contig', key: 'contigName'},
        {label: 'Position', key: 'position'},
        {label: 'Cigar String', key: 'cigarString'},
        {label: 'Num Mismatches', key: 'numMismatches'},
        {label: 'Num Gaps', key: 'numGaps'},
        {label: 'Score', key: 'score'},
        {label: 'Mapped By', key: 'mappedBy'},
      ],
      pagination: {
        itemsPerPage: 10,
        total: 0,
        itemsPerPageOptions: []
      }
    })

    const selectedRead = ref({})

    const readSummaryTableData = ref({items: []})

    let getReadSummaryTableData = function() {

      tableData.value.loading = true

      ApiService.get("/api/readSummaryTable")
          .then(response => {

            console.log(response)

            if (response.status === 200) {

              readSummaryTableData.value.items = response.data

              tableData.value.pagination.total = response.data.length

            } else {
              console.error("Failed to fetch read summary table data")
            }

            tableData.value.loading = false

          }).catch(err => {
        console.log(err);
      })
    }

    onMounted(() => {
      getReadSummaryTableData()
    })

    return {
      tableData,
      selectedRead,
      readSummaryTableData,
      getReadSummaryTableData
    }
  }
}
</script>

<style scoped>

</style>