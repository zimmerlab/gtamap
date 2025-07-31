<template>
  <UpsetPlot
    :upset-data="upsetData"
    :title="props.title"
    @on-set-click="onSetClick">
  </UpsetPlot>
</template>

<script setup>

import { ref, inject, onMounted, defineProps, defineEmits } from 'vue'
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

const fetchData = function() {
  ApiService.get(props.url)
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

const onSetClick = function(set) {
  emit("onSetClick", set)
}

onMounted(() => {
  fetchData()
})

</script>

<style scoped>

</style>