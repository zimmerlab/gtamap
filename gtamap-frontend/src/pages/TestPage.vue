<template>
  <UpsetPlotReadPos title="test title" url="/api/upsetDataRecordPosCigar"></UpsetPlotReadPos>
</template>

<script>
import UpsetPlot from "../components/UpsetPlot.vue";
import UpsetPlotMapperAgreement from "../components/UpsetPlotMapperAgreement.vue";

export default {
  name: 'TestPage',
  components: {UpsetPlotReadPos: UpsetPlotMapperAgreement, UpsetPlot},
  setup() {
    // This is where you can set up your component logic
    // For example, you can import and use the UpsetPlot component here
  }
}
</script>

<style scoped>

</style>