import { defineStore } from 'pinia'

export const useDataStore = defineStore('data', {
	state: () => ({
		reads: [],
		qnameToReadMap: {},
		summaryTableFilters: {
			hideAcceptedRecords: false,
			hideDiscardedRecords: false,
			hideNonInProgressRecords: false,
		},
		summaryTableSortKey: '',
		summaryIgv: null,
		acceptedIgv: null,
	}),
	getters: {
		getReads() {
			return this.reads
		},
		getRead: (state) => {
			return (qname) => state.qnameToReadMap[qname]
		},
		getSummaryTableReads() {
			return this.reads
				.filter((item) => {
					if (
						this.summaryTableFilters.hideAcceptedRecords &&
						item.isAcceptedR1 &&
						item.isAcceptedR2
					)
						return false
					if (this.summaryTableFilters.hideDiscardedRecords && item.isDiscarded)
						return false
					if (
						this.summaryTableFilters.hideNonInProgressRecords &&
						!(
							(item.isAcceptedR1 || item.isAcceptedR2) &&
							!(item.isAcceptedR1 && item.isAcceptedR2)
						)
					)
						return false
					return true
				})
				.sort((a, b) => {
					if (!this.summaryTableSortKey || this.summaryTableSortKey === '') {
						return 0
					}

					const desc = this.summaryTableSortKey[0] === '-'
					const key = this.summaryTableSortKey.slice(1)

					if (!key) {
						return 0
					}

					// Define which fields are numbers for proper sorting
					const numberFields = [
						'confidenceLevel',
						'numLocations',
						'numMappedBy',
					]
					const isNumber = numberFields.includes(key)

					const sortLikeNumber = (a, b) => {
						return a[key] < b[key] ? -1 : a[key] > b[key] ? 1 : 0
					}

					const sortLikeString = (a, b) => {
						return String(a[key] || '').localeCompare(String(b[key] || ''))
					}

					if (isNumber) {
						return sortLikeNumber(a, b) * (desc ? -1 : 1)
					} else {
						return sortLikeString(a, b) * (desc ? -1 : 1)
					}
				})
		},
		isSummaryTableFilterHideAcceptedRecords: (state) => {
			return state.summaryTableFilters.hideAcceptedRecords
		},
		isSummaryTableFilterHideDiscardedRecords: (state) => {
			return state.summaryTableFilters.hideDiscardedRecords
		},
		isSummaryTableFilterHideNonInProgressRecords: (state) => {
			return state.summaryTableFilters.hideNonInProgressRecords
		},
	},
	actions: {
		setReads(reads) {
			this.reads = reads

			for (const read of reads) {
				this.qnameToReadMap[read.qname] = read
			}
		},
		setSummaryIgv(igv) {
			this.summaryIgv = igv
		},
		setAcceptedIgv(igv) {
			this.acceptedIgv = igv
		},
		toggleSummaryTableFilterHideAcceptedRecords() {
			this.summaryTableFilters.hideAcceptedRecords =
				!this.summaryTableFilters.hideAcceptedRecords
		},
		toggleSummaryTableFilterHideDiscardedRecords() {
			this.summaryTableFilters.hideDiscardedRecords =
				!this.summaryTableFilters.hideDiscardedRecords
		},
		toggleSummaryTableFilterHideNonInProgressRecords() {
			this.summaryTableFilters.hideNonInProgressRecords =
				!this.summaryTableFilters.hideNonInProgressRecords
		},
		setSummaryTableFilterHideAcceptedRecords(value) {
			this.summaryTableFilters.hideAcceptedRecords = value
		},
		setSummaryTableFilterHideDiscardedRecords(value) {
			this.summaryTableFilters.hideDiscardedRecords = value
		},
		setSummaryTableFilterHideNonInProgressRecords(value) {
			this.summaryTableFilters.hideNonInProgressRecords = value
		},
		setSummaryTableSortKey(sortKey) {
			this.summaryTableSortKey = sortKey
		},
	},
})
