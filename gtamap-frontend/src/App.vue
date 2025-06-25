<template>
    <div class="tw:h-screen">
        <header class="no-shrink">
            <w-toolbar class="tw:h-10" fixed>
                <div class="tw:text-md tw:font-bold tw:w-20 tw:mx-8">
                    <router-link to="/">
                        GTAMap
                    </router-link>
                </div>
                <w-list :items="routes" nav color="primary" class="tw:py-2 tw:w-full w-flex row">
                    <template #item="{ item }">
                        <span class="tw:mx-2">{{ item.label }}</span>
                    </template>
                </w-list>
            </w-toolbar>
        </header>
        <div class="content-wrap w-flex no-shrink tw:mt-10">
            <div class="tw:container tw:mx-auto w-flex column grow">
                <router-view class="grow"></router-view>
            </div>
        </div>
    </div>
</template>

<script>

    import { ref } from "vue"

    export default {
        name: "App",
        setup() {
            const routes = ref([
                {
                    label: "Home",
                    route: "/"
                },
                {
                    label: "Other Page",
                    route: "/other"
                },
            ])

            return { routes }
        }
    }

</script>

<style scoped>

</style>
