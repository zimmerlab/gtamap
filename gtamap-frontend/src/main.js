import { createApp } from 'vue'

import WaveUI from 'wave-ui'
import 'wave-ui/dist/wave-ui.css'

import './style.css'

import App from './App.vue'

const app = createApp(App)

app.use(WaveUI, {})

import router from "./router.js"
app.use(router)

let API_BASEURL = location.protocol + '//' + location.hostname + ':8000'
//if (process.env.NODE_ENV === 'development') {
//    API_BASEURL = location.protocol + '//' + location.hostname + ':8080'
//}

import ApiService from "./services/api.service.js";
app.use(ApiService, {
    API_BASEURL
})

app.mount('#app')