<template>
  <div class="w-full h-full flex flex-col items-center justify-center">
    <MapperMultimappingParallelPlot></MapperMultimappingParallelPlot>
  </div>
</template>

<script>

import { ref, inject, onMounted } from 'vue'

// import MultidimensionalScalingPlot from "../components/MultidimensionalScalingPlot.vue";
import MapperMultimappingParallelPlot from "../components/MapperMultimappingParallelPlot.vue";

export default {
  name: 'OtherPage',
  components: {
    // MultidimensionalScalingPlot
    MapperMultimappingParallelPlot
  },
  setup() {

    onMounted(() => {
    })

    return {

    }
  }
}
</script>

<style scoped>

</style>