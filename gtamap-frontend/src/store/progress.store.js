import { defineStore } from "pinia"

export const useProgressStore = defineStore('progress', {
	state: () => ({
		readsTotal: 0,
		readsAccepted: 0,
		readsDiscarded: 0,
		readsInProgress: 0,
		readsTodo: 0,
	}),
	getters: {
		getReadsTotal: (state) => state.readsTotal,
		getReadsAccepted: (state) => state.readsAccepted,
		getReadsDiscarded: (state) => state.readsDiscarded,
		getReadsInProgress: (state) => state.readsInProgress,
		getReadsTodo: (state) => state.readsTodo,
	},
	actions: {
		setReadsTotal(value) {
			this.readsTotal = value
			this.readsTodo = value
		},
		incrementAccepted(count) {
			this.readsTodo -= count
			this.readsAccepted += count
		},
	},
})
