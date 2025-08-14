<template>
  <div ref="plotDiv"></div>
</template>

<script setup>

import { ref, onBeforeUnmount, defineProps, watch, defineEmits } from 'vue';

import * as UpSetJS from '@upsetjs/bundle'

const props = defineProps({
  upsetData: {
    type: Array,
    default: () => []
  },
  title: {
    type: String,
    default: "Upset Plot"
  }
})

const emit = defineEmits(['onSetClick'])

const plotDiv = ref(null);

const upsetDataRecordPos = ref({
  sets: [],
  combinations: []
});

const processUpsetData = function(elems) {

  const sets = UpSetJS.extractSets(elems, elem => elem.sets);
  const combinations = UpSetJS.generateCombinations(sets, {
    type: 'distinctIntersection'
  });

  // sort the combinations by cardinality (set size) and name
  const nameToCombination = {};
  for (const combination of combinations) {
    nameToCombination[combination.name] = combination;
  }
  const sortSets = (a, b) => {
    const sizeA = nameToCombination[a.name] ? nameToCombination[a.name].cardinality : 0;
    const sizeB = nameToCombination[b.name] ? nameToCombination[b.name].cardinality : 0;
    const sizeDiff = sizeB - sizeA;
    return sizeDiff !== 0 ? sizeDiff : a.name.localeCompare(b.name);
  };
  const sortCombinations = (a, b) => {
    const sizeDiff = b.cardinality - a.cardinality;
    return sizeDiff !== 0 ? sizeDiff : a.name.localeCompare(b.name);
  };
  sets.sort(sortSets);
  combinations.sort(sortCombinations);
  // end sort by cardinality and name

  upsetDataRecordPos.value.sets = sets;
  upsetDataRecordPos.value.combinations = combinations;

  renderUpsetRecordPos();
}

let selectionUpsetRecordPos = ref(null);

function onHoverUpsetRecordPos(set) {
  selectionUpsetRecordPos.value = set;
  renderUpsetRecordPos();
}

function onClickUpsetRecordPos(set) {
  emit('onSetClick', set);
}

function renderUpsetRecordPos() {

  if (!plotDiv.value) {
    console.error("Plot div is not available");
    return;
  }

  let sets = upsetDataRecordPos.value.sets
  let combinations = upsetDataRecordPos.value.combinations
  const upsetConfig = {
    sets,
    combinations,
    width: 800,
    height: 300,
    selection: selectionUpsetRecordPos.value,
    onHover: onHoverUpsetRecordPos,
    onClick: onClickUpsetRecordPos,
    color: "darkorchid",
    selectionColor: "#65cc32",
    //hoverHintColor: "#cc9932",
    //hasSelectionColor: "#cc9932",
    //alternatingBackgroundColor: true,
    //hasSelectionOpacity: 0.5,
    title: props.title,
    fontSizes: {
      barLabel: "8pt",
      chartLabel: "8pt",
      axisTick: "8pt",
      setLabel: "8pt",
      title: "11pt",
    }
  }
  UpSetJS.render(plotDiv.value, upsetConfig);
}

watch(() => props.upsetData, (newData) => {
  if (newData && newData.length > 0) {
    processUpsetData(newData);
  }
}, { immediate: true });

onBeforeUnmount(() => {
  if (plotDiv.value) {
    plotDiv.value.innerHTML = ''
  }
})

</script>

<style scoped>

</style>