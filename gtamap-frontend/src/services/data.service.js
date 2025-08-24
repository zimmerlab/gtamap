
import { useDataStore } from '@/store/data.store'

const DataService = {
	ApiService: undefined,
	init: function(ApiService) {
		this.ApiService = ApiService
	},
}

export default {
	install: (app, options) => {
		DataService.init(app.config.globalProperties.$http)
		app.provide('data', DataService)
	},
}
