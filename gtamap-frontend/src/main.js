import { createApp } from 'vue'
import { createPinia } from 'pinia'

import WaveUI from 'wave-ui'
import 'wave-ui/dist/wave-ui.css'

import './style.css'
import '@mdi/font/css/materialdesignicons.min.css'

import App from './App.vue'

const app = createApp(App)

const pinia = createPinia()
app.use(pinia)

app.use(WaveUI, {})

import router from "./router.js"
app.use(router)

let API_BASEURL = location.protocol + '//' + location.hostname + ':8000'
//if (process.env.NODE_ENV === 'development') {
//    API_BASEURL = location.protocol + '//' + location.hostname + ':8080'
//}

import ApiService from "./services/api.service.js"
app.use(ApiService, {
    API_BASEURL
})

import DataService from "@/services/data.service.js"
app.use(DataService, {})

app.mount('#app')
