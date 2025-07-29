<template>
  <div class="w-full h-full flex flex-col items-center justify-center">
<!--    <MapperMultimappingParallelPlot></MapperMultimappingParallelPlot>-->
<!--    <UmapTestPlot></UmapTestPlot>-->
<!--    <MultidimensionalScalingPlot></MultidimensionalScalingPlot>-->
    <MapperSimilarityUpgmaPlot></MapperSimilarityUpgmaPlot>
  </div>
</template>

<script>

import { ref, inject, onMounted } from 'vue'

// import UmapTestPlot from "../components/UmapTestPlot.vue";

// import MultidimensionalScalingPlot from "../components/MultidimensionalScalingPlot.vue";
// import MapperMultimappingParallelPlot from "../components/MapperMultimappingParallelPlot.vue";
import MapperSimilarityUpgmaPlot from "../components/MapperSimilarityUpgmaPlot.vue";

export default {
  name: 'OtherPage',
  components: {
    // MultidimensionalScalingPlot,
    // MapperMultimappingParallelPlot,
    // UmapTestPlot,
    MapperSimilarityUpgmaPlot
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