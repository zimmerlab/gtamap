<template>
  <div>
    <div class="tw:flex tw:justify-center tw:mb-4 tw:gap-8">
      <!-- Detail Level Group -->
      <div class="tw:flex tw:flex-col tw:items-center">
        <div class="tw:text-sm tw:font-medium tw:mb-2">Detail Level</div>
        <div class="tw:flex tw:gap-4">
          <w-checkbox v-model="checkboxes.read" label="Read" :disabled="true" />
          <w-checkbox v-model="checkboxes.contigPos" label="Contig + Position" :disabled="checkboxes.cigar" />
          <w-checkbox v-model="checkboxes.cigar" label="CIGAR" @update:model-value="onCigarChange" />
        </div>
      </div>

      <!-- Read Filter Group -->
      <div class="tw:flex tw:flex-col tw:items-center">
        <div class="tw:text-sm tw:font-medium tw:mb-2">Read Filter</div>
        <div class="tw:flex tw:gap-4">
          <w-checkbox v-model="checkboxes.targetRegion" label="Target Region Only" />
        </div>
      </div>
    </div>

    <UpsetPlot :upset-data="upsetData" :title="computedTitle" @onSetClick="onSetClick" @onSetHover="onSetHover">
    </UpsetPlot>
  </div>
</template>

<script setup>
import { ref, inject, onMounted, defineProps, defineEmits, computed, watch } from 'vue'
import UpsetPlot from './UpsetPlot.vue'

const props = defineProps({
  url: {
    type: String,
    default: '/api/upsetData',
  },
  title: {
    type: String,
    default: 'Upset Plot',
  },
})

const emit = defineEmits(['onSetClick', 'onSetHover'])

const DataService = inject('data')
const ApiService = inject('http')

const upsetData = ref([])

const checkboxes = ref({
  read: true,
  contigPos: false,
  cigar: false,
  targetRegion: false,
})

const computedTitle = computed(() => {
  const parts = []

  if (checkboxes.value.read) parts.push('Read')
  if (checkboxes.value.contigPos) parts.push('Contig + Position')
  if (checkboxes.value.cigar) parts.push('CIGAR')

  // return "Mapper Agreement (" + parts.join(" + ") + ")"
  return ''
})

const fetchData = function (url) {
  ApiService.get(props.url, {
    params: {
      onlyTargetRegion: checkboxes.value.targetRegion,
      usePosition: checkboxes.value.contigPos,
      useCigar: checkboxes.value.cigar,
    },
  })
    .then((response) => {
      if (response.status === 200) {
        upsetData.value = response.data
      } else {
        console.error('Failed to fetch upset data')
        console.error(response)
        upsetData.value = []
      }
    })
    .catch((err) => {
      console.log(err)
    })
}

watch(
  [
    () => checkboxes.value.targetRegion,
    () => checkboxes.value.contigPos,
    () => checkboxes.value.cigar,
  ],
  () => {
    fetchData()
  }
)

const onCigarChange = function (value) {
  if (value) {
    checkboxes.value.contigPos = true
  }
}

const onSetHover = function (set) {
  emit('onSetHover', set)
}

const onSetClick = function (set) {
  // emit('onSetClick', set)

  const qnames = set.elems.map((e) => {
    return e.name.split('::')[1]
  })

  DataService.setSummaryTableSpecificFilters(
    set.name,
    checkboxes.value.targetRegion,
    checkboxes.value.contigPos,
    checkboxes.value.cigar,
    qnames
  )
}

onMounted(() => {
  fetchData()
})
</script>

<style scoped></style>
