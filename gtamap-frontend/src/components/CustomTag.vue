<template>
  <div class="confidence-score tw:select-none" :class="['level-'+props.level, 'size-'+props.size]">
    <div class="signal-container">
      <div
          v-for="index in 5"
          :key="index"
          class="signal-bar"
          :style="signalStyle(index)"
      ></div>
    </div>
  </div>
</template>

<script setup>

import { ref, computed, defineProps } from 'vue'

const props = defineProps({
  level: {
    type: Number,
    default: 1
  },
  size: {
    type: String,
    default: 'md',
    validator: (value) => ['xs', 'sm', 'md', 'lg'].includes(value)
  }
})

const signalStyle = function(index) {
  let opacity = 0;
  if (index <= props.level) {
    opacity = 1;
  } else if (index === props.level + 1) {
    opacity = 0.3;
  } else {
    opacity = 0.1;
  }
  return "opacity: " + opacity + ";";
}

</script>

<style scoped>

.confidence-score {
  display: inline-flex;
  align-items: center;
  gap: 4px;
  font-weight: bold;
  position: relative;
  overflow: hidden;
  box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
}

/* Size variants */
.size-xs {
  padding: 5px 6px;
  font-size: 9px;
  border-radius: 2px;
}

.size-sm {
  padding: 4px 8px;
  font-size: 10px;
  border-radius: 3px;
}

.size-md {
  padding: 9px 20px;
  font-size: 12px;
  border-radius: 4px;
}

.size-lg {
  padding: 12px 24px;
  font-size: 14px;
  border-radius: 5px;
}

.signal-container {
  display: flex;
  align-items: flex-end;
}

.signal-bar {
  background: white;
  border-radius: 2px;
  margin: 0 1px;
}

/* Extra small signal bar sizes */
.size-xs .signal-bar {
  width: 4px;
  height: 16px;
}

.size-xs .signal-bar:nth-child(1) { height: 5px; }
.size-xs .signal-bar:nth-child(2) { height: 6px; }
.size-xs .signal-bar:nth-child(3) { height: 7px; }
.size-xs .signal-bar:nth-child(4) { height: 8px; }
.size-xs .signal-bar:nth-child(5) { height: 9px; }

/* Small signal bar sizes */
.size-sm .signal-bar {
  width: 5px;
  height: 12px;
}

.size-sm .signal-bar:nth-child(1) { height: 5px; }
.size-sm .signal-bar:nth-child(2) { height: 7px; }
.size-sm .signal-bar:nth-child(3) { height: 9px; }
.size-sm .signal-bar:nth-child(4) { height: 11px; }
.size-sm .signal-bar:nth-child(5) { height: 13px; }

/* Default (md) signal bar sizes */
.size-md .signal-bar {
  width: 8px;
  height: 20px;
}

.size-md .signal-bar:nth-child(1) { height: 8px; }
.size-md .signal-bar:nth-child(2) { height: 12px; }
.size-md .signal-bar:nth-child(3) { height: 16px; }
.size-md .signal-bar:nth-child(4) { height: 20px; }
.size-md .signal-bar:nth-child(5) { height: 24px; }

/* Large signal bar sizes */
.size-lg .signal-bar {
  width: 10px;
  height: 24px;
}

.size-lg .signal-bar:nth-child(1) { height: 10px; }
.size-lg .signal-bar:nth-child(2) { height: 14px; }
.size-lg .signal-bar:nth-child(3) { height: 18px; }
.size-lg .signal-bar:nth-child(4) { height: 22px; }
.size-lg .signal-bar:nth-child(5) { height: 26px; }

/* Level 1 - Very Low */
.level-1 {
  background: linear-gradient(135deg, #ff4444, #cc2222);
  color: white;
}

/* Level 2 - Low */
.level-2 {
  background: linear-gradient(135deg, #ff8800, #cc6600);
  color: white;
}



/* Level 3 - Medium (Bronze) */
.level-3 {
  background: linear-gradient(135deg, #cd7f32, #b8722c, #a0601f);
  background-size: 200% 200%;
  color: white;
  box-shadow:
      0 0 15px rgba(205, 127, 50, 0.4),
      0 6px 20px rgba(205, 127, 50, 0.2),
      inset 0 1px 0 rgba(255, 255, 255, 0.3);
}

/* Level 4 - High (Silver) */
.level-4 {
  background: linear-gradient(135deg, #e8e8e8, #d4d4d4, #c0c0c0, #b8b8b8);
  background-size: 200% 200%;
  color: #333;
  box-shadow:
      0 0 20px rgba(232, 232, 232, 0.6),
      0 8px 25px rgba(212, 212, 212, 0.4),
      inset 0 1px 0 rgba(255, 255, 255, 0.8),
      inset 0 -1px 0 rgba(0, 0, 0, 0.1);
}

.level-4::after {
  content: '';
  position: absolute;
  top: 0;
  left: 0;
  right: 0;
  bottom: 0;
  background: linear-gradient(45deg,
  transparent 30%,
  rgba(255, 255, 255, 0.3) 50%,
  transparent 70%);
  animation: shine 3.5s linear infinite;
  pointer-events: none;
}

/* Level 5 - Very High (Gold) */
.level-5 {
  background: linear-gradient(135deg, #ffd700, #ffb300, #ffa500, #ff8c00);
  background-size: 200% 200%;
  color: #1a1a1a;
  box-shadow:
      0 0 20px rgba(255, 215, 0, 0.5),
      0 8px 25px rgba(255, 215, 0, 0.3),
      inset 0 1px 0 rgba(255, 255, 255, 0.6);
  animation: shimmer 2s ease-in-out infinite alternate;
  position: relative;
}

.level-5::after {
  content: '';
  position: absolute;
  top: 0;
  left: 0;
  right: 0;
  bottom: 0;
  background: linear-gradient(45deg,
  transparent 30%,
  rgba(255, 255, 255, 0.3) 50%,
  transparent 70%);
  animation: shine 3.5s linear infinite;
  pointer-events: none;
}

@keyframes shimmer {
  0% {
    background-position: 0% 50%;
  }
  100% {
    background-position: 100% 50%;
  }
}

@keyframes shine {
  0% {
    transform: translateX(-100%) skewX(-15deg);
  }
  50% {
    transform: translateX(100%) skewX(-15deg);
  }
  100% {
    transform: translateX(100%) skewX(-15deg);
  }
}

.confidence-text {
  margin-right: 8px;
}

.confidence-level {
  font-size: 18px;
  font-weight: 900;
}
</style>