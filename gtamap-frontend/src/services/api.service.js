import axios from 'axios'

const ApiService = {
    init: function (baseURL) {
        console.log('initialized api service with base url: ' + baseURL)

        this.baseUrl = baseURL

        axios.defaults.baseURL = baseURL
    },
    getBaseUrl: function () {
        return this.baseUrl
    },
    get: function (url, params = {}) {
        return axios.get(url, params)
    },
    post: function (url, body) {
        return axios.post(url, body)
    }
}

export default {
    install: (app, options) => {
        ApiService.init(options.API_BASEURL)
        app.provide("http", ApiService)
    }
}