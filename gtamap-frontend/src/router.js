import { createRouter, createWebHistory } from 'vue-router'

import IndexPage from "./pages/IndexPage.vue";
import OtherPage from "./pages/OtherPage.vue";
import OverviewPage from "./pages/OverviewPage.vue";

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
        {
            path: '/overview',
            name: 'overview',
            component: OverviewPage
        },
    ]
})

export default router