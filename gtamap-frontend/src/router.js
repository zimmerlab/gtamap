import { createRouter, createWebHistory } from 'vue-router'

// import IndexPage from "./pages/IndexPage.vue"
// import OtherPage from "./pages/OtherPage.vue"
import OverviewPage from "./pages/OverviewPage.vue"
// import TablePage from "./pages/TablePage.vue"
import TestPage from "./pages/TestPage.vue"
import ReadDetailsPage from "./pages/ReadDetailsPage.vue"

const router = createRouter({
    history: createWebHistory(),
    base: '/',
    routes: [
			{
					path: '/',
					name: 'index',
					component: OverviewPage
			},
			{
					path: '/test',
					name: 'test',
					component: TestPage
			},
			{
				path: '/readDetails',
				name: 'readDetailsPage',
				component: ReadDetailsPage
			}
        // {
        //     path: '/overview',
        //     name: 'overview',
        //     component: OverviewPage
        // },
    ]
})

export default router
