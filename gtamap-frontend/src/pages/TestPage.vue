<template>
  <UpsetPlotReadPos title="test title" url="/api/upsetDataRecordPosCigar"></UpsetPlotReadPos>
</template>

<script>
import UpsetPlot from "../components/UpsetPlot.vue";
import UpsetPlotReadPos from "../components/UpsetPlotReadPos.vue";

export default {
  name: 'TestPage',
  components: {UpsetPlotReadPos, UpsetPlot},
  setup() {
    // This is where you can set up your component logic
    // For example, you can import and use the UpsetPlot component here
  }
}
</script>

<style scoped>

</style>