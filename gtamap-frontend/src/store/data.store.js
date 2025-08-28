import { defineStore } from 'pinia'

export const useDataStore = defineStore('data', {
	state: () => ({
		reads: [],
		qnameToReadMap: {},
		indexToRecordMap: {},
		summaryTableFilters: {
			hideAcceptedRecords: false,
			hideDiscardedRecords: false,
			hideNonInProgressRecords: false,
		},
		summaryTableSortKey: '',
		igvSummary: {
			div: null,
			browser: null,
			config: null,
		},
		igvAccepted: {
			div: null,
			browser: null,
			config: null,
		},
	}),
	getters: {
		getReads() {
			return this.reads
		},
		getReadByQname: (state) => {
			return (qname) => state.qnameToReadMap[qname]
		},

		getRecordByIndex: (state) => {
			return (index) => state.indexToRecordMap[index]
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
		getSummaryTableFilters: (state) => {
			return state.summaryTableFilters
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
		getIgvSummary: (state) => {
			return state.igvSummary
		},
		getIgvAccepted: (state) => {
			return state.igvAccepted
		},
	},
	actions: {
		setReads(reads) {
			this.reads = reads

			for (const read of reads) {
				this.qnameToReadMap[read.qname] = read

				for (const record of read.locations) {
					for (const index of record.readIndices) {
						this.indexToRecordMap[index] = record
					}
				}
			}
		},
		setIgvSummary(igvDiv, igvBrowser, config) {
			this.igvSummary.div = igvDiv
			this.igvSummary.browser = igvBrowser
			this.igvSummary.config = config
		},
		setIgvAccepted(igvDiv, igvBrowser, config) {
			this.igvAccepted.div = igvDiv
			this.igvAccepted.browser = igvBrowser
			this.igvAccepted.config = config
		},
		toggleSummaryTableFiltersHideAcceptedRecords() {
			this.summaryTableFilters.hideAcceptedRecords =
				!this.summaryTableFilters.hideAcceptedRecords
		},
		toggleSummaryTableFiltersHideDiscardedRecords() {
			this.summaryTableFilters.hideDiscardedRecords =
				!this.summaryTableFilters.hideDiscardedRecords
		},
		toggleSummaryTableFiltersHideNonInProgressRecords() {
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
