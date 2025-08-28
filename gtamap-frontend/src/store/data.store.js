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
		summaryTableSpecificFilters: {
			qnames: new Set(),
			name: '',
			targetRegion: false,
			contigPos: false,
			cigar: false,
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
					) {
						return false
					}
					if (
						this.summaryTableFilters.hideDiscardedRecords &&
						item.isDiscarded
					) {
						return false
					}
					if (
						this.summaryTableFilters.hideNonInProgressRecords &&
						!(
							(item.isAcceptedR1 || item.isAcceptedR2) &&
							!(item.isAcceptedR1 && item.isAcceptedR2)
						)
					) {
						return false
					}

					if (
						this.summaryTableSpecificFilters.qnames &&
						this.summaryTableSpecificFilters.qnames.size > 0 &&
						!this.summaryTableSpecificFilters.qnames.has(item.qname)
					) {
						return false
					}

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

		isProgressLoaded: (state) => {
			return state.reads && state.reads.length > 0
		},
		numReadsTotal: (state) => {
			return state.reads.length
		},
		numReadsAccepted: (state) => {
			return state.reads.filter(
				(read) => read.isAcceptedR1 && read.isAcceptedR2
			).length
		},
		numReadsDiscarded: (state) => {
			return state.reads.filter((read) => read.isDiscarded).length
		},
		numReadsInProgress: (state) => {
			return state.reads.filter(
				(read) =>
					!(read.isAcceptedR1 && read.isAcceptedR2) &&
					(read.isAcceptedR1 || read.isAcceptedR2)
			).length
		},
		numReadsTodo: (state) => {
			return (
				state.numReadsTotal - (state.numReadsAccepted + state.numReadsDiscarded)
			)
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

		getSummaryTableSpecificFilters: (state) => {
			return state.summaryTableSpecificFilters
		},

		getSummaryTableSpecificFilterSetName: (state) => {
			return state.summaryTableSpecificFilters.name
		},
		getSummaryTableSpecificFilterAsString: (state) => {
			if (
				!state.summaryTableSpecificFilters.qnames ||
				state.summaryTableSpecificFilters.qnames.size === 0
			) {
				return undefined
			}

			let s = 'Read Name'
			if (state.summaryTableSpecificFilters.contigPos) {
				s += ' and Contig + Position'
			}
			if (state.summaryTableSpecificFilters.cigar) {
				s += ' and CIGAR'
			}
			if (state.summaryTableSpecificFilters.targetRegion) {
				s += ' only in Target Region'
			}

			return s
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
		resetSummaryTableFilters() {
			this.summaryTableFilters.hideAcceptedRecords = false
			this.summaryTableFilters.hideDiscardedRecords = false
			this.summaryTableFilters.hideNonInProgressRecords = false
		},

		setSummaryTableSpecificFilters(name, targetRegion, contigPos, cigar, qnames) {
			this.summaryTableSpecificFilters.name = name
			this.summaryTableSpecificFilters.targetRegion = targetRegion
			this.summaryTableSpecificFilters.contigPos = contigPos
			this.summaryTableSpecificFilters.cigar = cigar
			this.summaryTableSpecificFilters.qnames = new Set(qnames)
		},
		resetSummaryTableSpecificFilters() {
			this.summaryTableSpecificFilters.name = ''
			this.summaryTableSpecificFilters.targetRegion = false
			this.summaryTableSpecificFilters.contigPos = false
			this.summaryTableSpecificFilters.cigar = false
			this.summaryTableSpecificFilters.qnames = new Set()
		},

		setSummaryTableSortKey(sortKey) {
			this.summaryTableSortKey = sortKey
		},
	},
})
