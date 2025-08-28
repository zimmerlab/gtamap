import { useDataStore } from '@/store/data.store'
import { useProgressStore } from '@/store/progress.store'

const DataService = {
	ApiService: undefined,
	DataStore: undefined,
	ProgressStore: undefined,

	init: function(ApiService) {
		this.ApiService = ApiService
		this.DataStore = useDataStore()
		this.ProgressStore = useProgressStore()
	},

	fetchGeneralInfo: async function() {
		const response = await this.ApiService.get('/api/general/info')
		if (!response || response.status !== 200 || !response.data) {
			console.error('Error fetching general info')
			console.error(response)
			return
		}

		console.log(response.data)
	},

	fetchReads: async function() {
		const response = await this.ApiService.get('/api/summary/table')
		if (response && response.status === 200 && response.data) {
			// set parents
			for (const item of response.data) {
				item.isSelected = false
				for (const l of item.locations) {
					l.parent = item
				}
			}
			this.DataStore.setReads(response.data)
		}
	},

	updateSummaryIgvTrack: function() {
		const igvSummary = this.DataStore.getIgvSummary
		const tracks = igvSummary.browser.findTracks()
		const track = tracks[tracks.length - 1]
		igvSummary.browser.removeTrack(track)
		igvSummary.browser.loadTrack(igvSummary.config.tracks[0])
	},

	updateAcceptedIgvTrack: function() {
		const igvAccepted = this.DataStore.getIgvAccepted
		const tracks = igvAccepted.browser.findTracks()
		const track = tracks[tracks.length - 1]
		igvAccepted.browser.removeTrack(track)
		igvAccepted.browser.loadTrack(igvAccepted.config.tracks[0])
	},

	acceptAllMaxConfidenceReads: function() {
		const records = this.DataStore.getReads
			.filter((read) => read.confidenceLevel === 5 && !read.isAccepted)
			.map((read) => read.locations)
			.flat()

		this.acceptRecords(records)

		this.updateSummaryIgvTrack()
		this.updateAcceptedIgvTrack()
	},

	updateRecordAcceptance: function(record, isAccepted) {
		// get records of the same read with the same pair type
		const recordsToUnaccept = record.parent.locations.filter(
			(o) => o.pairType === record.pairType && o.isAccepted
		)

		this.unacceptRecords(recordsToUnaccept)

		if (isAccepted) {
			this.acceptRecords([record])
		}

		this.updateSummaryIgvTrack()
		this.updateAcceptedIgvTrack()
	},

	updateRecordAcceptanceByIndex: function(index, isAccepted) {
		const record = this.DataStore.getRecordByIndex(index)
		this.updateRecordAcceptance(record, isAccepted)
	},
	discardReadByQname: function(qname) {
		const read = this.DataStore.getReadByQname(qname)
		this.discardReads([read])
	},
	resetReadByQname: function(qname) {
		const read = this.DataStore.getReadByQname(qname)
		this.resetReads([read])
	},

	acceptRecords: function(records) {
		for (const record of records) {
			record.isAccepted = true
			if (record.pairType === 'first') {
				record.parent.isAcceptedR1 = true
			} else {
				record.parent.isAcceptedR2 = true
			}
		}

		const recordIds = records.map((r) => r.readIndices).flat()

		this.ApiService.post('/api/summary/acceptRecords', {
			recordIds: recordIds,
		})
			.then((response) => {
				// console.log(response)
			})
			.catch((err) => {
				console.error('Error accepting record:', err)
			})
	},

	unacceptRecords: function(records) {
		for (const record of records) {
			record.isAccepted = false
		}

		for (const record of records) {
			const isAnyAccepted = record.parent.locations.find(
				(o) => o.pairType === record.pairType && o.isAccepted
			)

			if (record.pairType === 'first') {
				record.parent.isAcceptedR1 = isAnyAccepted ? true : false
			} else {
				record.parent.isAcceptedR2 = isAnyAccepted ? true : false
			}
		}

		const recordIds = records.map((r) => r.readIndices).flat()

		this.ApiService.post('/api/summary/unacceptRecords', {
			recordIds: recordIds,
		})
			.then((response) => {
				// console.log(response)
			})
			.catch((err) => {
				console.error('Error unaccepting record:', err)
			})
	},

	discardReads: function(reads) {
		for (const read of reads) {
			read.isDiscarded = true
			read.isAcceptedR1 = false
			read.isAcceptedR2 = false
			for (const l of read.locations) {
				l.isAccepted = false
			}
		}

		const qnames = reads.map((r) => r.qname).flat()

		this.ApiService.post('/api/summary/discardReads', {
			qnames: qnames,
		})
			.then((response) => {
				// console.log(response)
			})
			.catch((err) => {
				console.error('Error discarding records:', err)
			})
	},

	resetReads: function(reads) {
		for (const read of reads) {
			read.isDiscarded = false
			read.isAcceptedR1 = false
			read.isAcceptedR2 = false
			for (const l of read.locations) {
				l.isAccepted = false
			}
		}

		const qnames = reads.map((r) => r.qname).flat()

		this.ApiService.post('/api/summary/resetReads', {
			qnames: qnames,
		})
			.then((response) => {
				// console.log(response)
			})
			.catch((err) => {
				console.error('Error resetting records:', err)
			})
	},

	toggleSummaryTableFiltersHideAcceptedRecords: function() {
		this.DataStore.toggleSummaryTableFiltersHideAcceptedRecords()

		this.updateSummaryTableFilters()
		this.updateSummaryIgvTrack()
	},
	toggleSummaryTableFiltersHideDiscardedRecords: function() {
		this.DataStore.toggleSummaryTableFiltersHideDiscardedRecords()

		this.updateSummaryTableFilters()
		this.updateSummaryIgvTrack()
	},
	toggleSummaryTableFiltersHideNonInProgressRecords: function() {
		this.DataStore.toggleSummaryTableFiltersHideNonInProgressRecords()

		this.updateSummaryTableFilters()
		this.updateSummaryIgvTrack()
	},

	updateSummaryTableFilters: function() {
		const filterOptions = this.DataStore.getSummaryTableFilters

		this.ApiService.post('/api/summary/filterUpdate', {
			hideAccepted: filterOptions.hideAcceptedRecords,
			hideDiscarded: filterOptions.hideDiscardedRecords,
			hideNonInProgress: filterOptions.hideNonInProgressRecords,
		})
	},
}

export default {
	install: (app, options) => {
		DataService.init(app.config.globalProperties.$http)
		app.provide('data', DataService)
	},
}
