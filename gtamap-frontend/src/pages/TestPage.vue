<template>
  <div class="tw:m-20">
    <CustomTag :level="1" size="sm"></CustomTag>
  </div>
  <div class="tw:m-20">
    <CustomTag :level="2" size="xs"></CustomTag>
  </div>
  <div class="tw:m-20">
    <CustomTag :level="3"></CustomTag>
  </div>
  <div class="tw:m-20">
    <CustomTag :level="4"></CustomTag>
  </div>
  <div class="tw:m-20">
    <CustomTag :level="5" size="xs"></CustomTag>
    <CustomTag :level="5" size="sm"></CustomTag>
    <CustomTag :level="5" size="md"></CustomTag>
  </div>
</template>

<script>
import CustomTag from "../components/CustomTag.vue";

export default {
  name: 'TestPage',
  components: {CustomTag},
  setup() {
  }
}
</script>

<style scoped>

</style>
