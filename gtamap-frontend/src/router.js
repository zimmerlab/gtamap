import { createRouter, createWebHistory } from 'vue-router'

import IndexPage from "./pages/IndexPage.vue";
import OtherPage from "./pages/OtherPage.vue";

const router = createRouter({
    history: createWebHistory(),
    base: '/',
    routes: [
        {
            path: '/',
            name: 'index',
            component: IndexPage
        },
        {
            path: '/other',
            name: 'other',
            component: OtherPage
        },
    ]
})

export default router