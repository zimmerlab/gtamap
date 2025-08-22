<template>
  <div>
    <div
      class="filter-section tw:mb-4 tw:p-4 tw:border tw:border-gray-300 tw:rounded-lg tw:bg-gray-50"
    >
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
        <w-button
          @click="applyFilters"
          class="tw:mt-2"
        >
          <span class="tw-text-xs">Apply Filters</span>
        </w-button>
      </div>
    </div>

    <div
      class="filter-section tw:mb-4 tw:px-4 tw:py-2 tw:border tw:border-gray-300 tw:rounded-lg tw:bg-gray-50"
    >
      <h3 class="tw:text-md tw:font-medium tw:mb-3">Quick Actions</h3>
      <div class="tw:flex tw:items-center tw:gap-4">
        <w-button @click="acceptAllMaxConfidenceReads"
          ><span class="tw-text-xs">Accept All Max Confidence Reads</span></w-button
        >
      </div>
    </div>

    <div
      class="filter-section tw:mb-4 tw:px-4 tw:py-2 tw:border tw:border-gray-300 tw:rounded-lg tw:bg-gray-50"
    >
      <h3 class="tw:text-md tw:font-medium tw:mb-3">Selection</h3>
      <div class="tw:flex tw:items-center tw:gap-4">
        <w-button @click="deselectAll"><span class="tw-text-xs">Unselect All</span> </w-button>
      </div>
      <div class="tw:flex tw:items-center tw:gap-4 tw:pt-2">
        <w-button><span class="tw-text-xs">Reset Acceptance for Selected</span></w-button>
        <w-button
          :disabled="!tableData.actions.accept"
          @click="acceptSelected"
          ><span class="tw-text-xs">Accept Selected</span>
        </w-button>
        <w-button><span class="tw-text-xs">Discard Selected</span></w-button>
        <w-button @click="resetSelected"><span class="tw-text-xs">Reset Selected</span> </w-button>
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
      expandable-rows
    >
      <template #item-cell.isSelected="{ item }">
        <w-checkbox
          :model-value="item.isSelected"
          @update:model-value="(value) => toggleItemSelection(item, value)"
        />
      </template>

      <template #item-cell.actions="{ item }">
        <w-button
          @click.stop="openReadDetails(item)"
          outline
          xs
        >
          <w-icon class="mr1">mdi mdi-open-in-new</w-icon>
          view details
        </w-button>
        <w-button
          v-if="!item.isDiscarded"
          fg-color="red"
          class="tw:ml-2"
          @click.stop="discardRead(item)"
          outline
          xs
          >discard read</w-button
        >
        <w-button
          v-else
          class="tw:ml-2"
          @click.stop="undiscardRead(item)"
          outline
          xs
          >reset</w-button
        >
      </template>

      <template #item-cell.state="{ item }">
        <w-tag
          v-if="item.isDiscarded"
          xs
          bg-color="red"
          >discarded</w-tag
        >
        <w-tag
          v-else-if="item.isAcceptedR1 && item.isAcceptedR2"
          xs
          bg-color="green"
          >accepted</w-tag
        >
        <w-tag
          v-else-if="item.isAcceptedR1 || item.isAcceptedR2"
          xs
          bg-color="yellow"
          >in progress</w-tag
        >
      </template>

      <template #item-cell.confidence="{ item }">
        <CustomTag
          :level="item.confidenceLevel"
          size="xs"
        ></CustomTag>
      </template>

      <template #item-cell.mappedBy="{ item }">
        <span v-if="Array.isArray(item.mappedBy)">
          <w-tag
            v-for="(val, idx) in item.mappedBy"
            :key="idx"
            class="mr4"
            color="primary"
            >{{ val }}</w-tag
          >
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

          <template #item-cell.targetRegionOverlap="{ item }">
            <w-tag
              v-if="item.targetRegionOverlap === 0"
              xs
              bg-color="red"
              >no overlap</w-tag
            >
            <w-tag
              v-else-if="item.targetRegionOverlap === 1"
              xs
              bg-color="green"
              >contained</w-tag
            >
            <w-tag
              v-else
              xs
              bg-color="yellow"
              >partially</w-tag
            >
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
                >{{ val }}</w-tag
              >
            </span>
            <span v-else>{{ item.mappedBy }}</span>
          </template>
        </w-table>
      </template>
    </w-table>
  </div>
</template>

<script setup>
  import CustomTag from './CustomTag.vue'

  import { defineExpose, inject, nextTick, onMounted, ref, defineEmits, computed } from 'vue'

  const emit = defineEmits([
    'open-read-details',
    'content-changed',
    'summary-table-update',
    'accepted-table-update',
  ])

  const ApiService = inject('http')

  const props = defineProps({
    url: {
      type: String,
      required: true,
    },
  })

  const pagination = ref({
    total: 0,
    itemsPerPage: 30,
    itemsPerPageOptions: [{ value: 15 }, { value: 30 }],
  })

  const tableData = ref({
    loading: true,
    sortKey: '+qname',
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
  })

  const handleSort = function (sortInfo) {
    tableData.value.sortKey = sortInfo[0] ? sortInfo[0] : ''
    customSort()
  }

  const customSort = function () {
    const items = summaryTableData.value.items

    if (!tableData.value.sortKey || tableData.value.sortKey === '') {
      return items
    }

    const desc = tableData.value.sortKey[0] === '-'
    const key = tableData.value.sortKey.slice(1)

    if (!key) {
      return items
    }

    const sortLikeNumber = function (a, b) {
      return a[key] < b[key] ? -1 : a[key] > b[key] ? 1 : 0
    }

    const sortLikeString = function (a, b) {
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
    filteredItems: [],
  })
  // content of the inner table (mapping location records per qname)
  const readMappingTableData = ref({ items: [] })

  const filterOptions = ref({
    hideAcceptedRecords: false,
    hideDiscardedRecords: false,
    hideNonInProgressRecords: false,
  })

  const applyFilters = function () {
    let items = summaryTableData.value.items

    if (filterOptions.value.hideAcceptedRecords) {
      items = items.filter((item) => !(item.isAcceptedR1 && item.isAcceptedR2))
    }
    if (filterOptions.value.hideDiscardedRecords) {
      items = items.filter((item) => !item.isDiscarded)
    }
    if (filterOptions.value.hideNonInProgressRecords) {
      items = items.filter(
        (item) =>
          (item.isAcceptedR1 || item.isAcceptedR2) && !(item.isAcceptedR1 && item.isAcceptedR2)
      )
    }

    // emit the table update only once the server was notified about the new filters
    updateSummaryFilters().then(() => {
      emit('summary-table-update')
    })

    deselectRows()

    summaryTableData.value.filteredItems = items
  }

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

  let getSummaryTableData = function () {
    tableData.value.loading = true

    ApiService.get(props.url)
      .then((response) => {
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
          console.error('Failed to fetch read summary table data')
        }

        tableData.value.loading = false
      })
      .catch((err) => {
        console.log(err)
      })
  }

  const updateSummaryFilters = function () {
    return ApiService.post('/api/summary/filterUpdate', {
      hideAccepted: filterOptions.value.hideAcceptedRecords,
      hideDiscarded: filterOptions.value.hideDiscardedRecords,
      hideNonInProgress: filterOptions.value.hideNonInProgressRecords,
    })
  }

  let openReadDetails = function (readItem) {
    emit('open-read-details', readItem)
  }

  const selectAndScrollToRead = async function (readInfo) {
    const targetRead = summaryTableData.value.items.find((item) => item.qname === readInfo.qname)
    if (!targetRead) {
      console.warn('Read not found in the table data')
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
          if (
            mapping.mappedBy[i] === readInfo.mapperNames[j] &&
            mapping.readIndices[i] === readInfo.samIndices[j]
          ) {
            if (selectedMappingRow !== undefined) {
              console.warn('Multiple mappings found for the selected read, selecting the first one')
              break
            }
            selectedMappingRow = mapping
            break
          }
        }
      }
    }

    if (!selectedMappingRow) {
      console.warn('No matching mapping found for the selected read')
      return
    }

    tableData.value.selectedRowsInExpanded = [selectedMappingRow._uid]
  }

  const isItemSelected = function (item) {
    return tableData.value.customSelectedItems.includes(item._uid)
  }

  const discardRead = function (item) {
    discardReads([item])
  }

  const discardReads = function (items) {
    for (const item of items) {
      item.isDiscarded = true
      item.isAcceptedR1 = false
      item.isAcceptedR2 = false
      for (const l of item.locations) {
        l.isAccepted = false
      }
    }

    const qnames = items.map((item) => item.qname).flat()
    ApiService.post('/api/summary/discardReads', {
      qnames: qnames,
    })
      .then((response) => {
        console.log(response)
      })
      .catch((err) => {
        console.error('Error discarding records:', err)
      })
  }

  const resetRead = function (item) {
    resetReads([item])
  }

  const resetReads = function (items) {
    for (const item of items) {
      item.isDiscarded = false
      item.isAcceptedR1 = false
      item.isAcceptedR2 = false

      for (const l of item.locations) {
        l.isAccepted = false
      }
    }

    emit('accepted-table-update')

    const qnames = items.map((item) => item.qname).flat()
    ApiService.post('/api/summary/resetReads', {
      qnames: qnames,
    })
      .then((response) => {
        console.log(response)
      })
      .catch((err) => {
        console.error('Error resetting discarded reads:', err)
      })
  }

  const acceptReads = function (reads) {
    const records = reads.map((item) => item.locations).flat()
    acceptRecords(records)
  }

  const acceptRecord = function (record) {
    acceptRecords([record])
  }

  const acceptRecords = function (records) {
    for (const r of records) {

      r.parent.isDiscarded = false

      r.isAccepted = true

      if (r.pairType === 'first') {
        r.parent.isAcceptedR1 = true
      } else if (r.pairType === 'second') {
        r.parent.isAcceptedR2 = true
      }
    }

    emit('accepted-table-update')

    const recordIds = records.map((r) => r.readIndices).flat()

    ApiService.post('/api/summary/acceptRecords', {
      recordIds: recordIds,
    })
      .then((response) => {
        console.log(response)
      })
      .catch((err) => {
        console.error('Error accepting records:', err)
      })
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

  const updateItemAcceptance = function (item, isAccepted) {
    const recordsToUnaccept = item.parent.locations.filter(
      (o) => o.pairType === item.pairType && o.isAccepted
    )

    unacceptRecords(recordsToUnaccept)

    if (isAccepted) {
      item.isAccepted = isAccepted
      acceptRecord(item)
    }
  }

  // QUICK ACTIONS

  const acceptAllMaxConfidenceReads = function () {
    const reads = summaryTableData.value.items.filter(
      (item) => item.confidenceLevel === 5 && !item.isAccepted
    )
    const records = reads.map((item) => item.locations).flat()
    acceptRecords(records)
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
    let items = []
    const selected = summaryTableData.value.filteredItems.filter((item) => item.isSelected)

    acceptReads(selected)
  }

  const resetSelected = function () {
    const selected = summaryTableData.value.filteredItems.filter((item) => item.isSelected)

    resetReads(selected)
  }

  defineExpose({ selectAndScrollToRead })

  onMounted(() => {
    getSummaryTableData()
  })
</script>

<style scoped></style>
