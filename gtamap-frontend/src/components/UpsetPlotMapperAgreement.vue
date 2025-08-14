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
    
    <UpsetPlot
      :upset-data="upsetData"
      :title="computedTitle"
      @on-set-click="onSetClick">
    </UpsetPlot>
  </div>
</template>

<script setup>

import { ref, inject, onMounted, defineProps, defineEmits, computed, watch } from 'vue'
import UpsetPlot from "./UpsetPlot.vue";

const props = defineProps({
  url: {
    type: String,
    default: "/api/upsetDataRecordPos"
  },
  title: {
    type: String,
    default: "Upset Plot"
  }
})

const emit = defineEmits(['onSetClick'])

const ApiService = inject("http")

const upsetData = ref([])

const checkboxes = ref({
  read: true,
  contigPos: false,
  cigar: false,
  targetRegion: false
})

const computedUrl = computed(() => {
  if (checkboxes.value.cigar) {
    return "/api/upsetDataRecordPosCigar"
  }
  if (checkboxes.value.contigPos) {
    return "/api/upsetDataRecordPos"
  }
  return "/api/upsetDataRead"
})

const computedTitle = computed(() => {
  const parts = []
  
  if (checkboxes.value.read) parts.push("Read")
  if (checkboxes.value.contigPos) parts.push("Contig + Position")
  if (checkboxes.value.cigar) parts.push("CIGAR")
  
  return "Mapper Agreement (" + parts.join(" + ") + ")"
})

const fetchData = function() {
  ApiService.get(computedUrl.value, {
    params: {
      onlyTargetRegion: checkboxes.value.targetRegion
    }
  })
    .then(response => {
      if (response.status === 200) {
        upsetData.value = response.data
      } else {
        console.error("Failed to fetch upset data")
        console.error(response)
        upsetData.value = []
      }
    }).catch(err => {
      console.log(err);
    })
}

watch([computedUrl, () => checkboxes.value.targetRegion], () => {
  fetchData()
})

const onCigarChange = function(value) {
  if (value) {
    checkboxes.value.contigPos = true
  }
}

const onSetClick = function(set) {
  emit("onSetClick", set)
}

onMounted(() => {
  fetchData()
})

</script>

<style scoped>

</style>