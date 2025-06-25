<template>
    <div>
        <div>Index Page</div>
        <div id="igv-div"></div>
    </div>
</template>

<script>

import { ref, inject } from "vue"
import igv from "../js/igv/dist/igv.esm.js"

export default {
    name: "IndexPage",
    components: {},
    methods: {},
    setup() {

        let igvInfo = ref({})
        let igvBrowser = ref(undefined)

        const ApiService = inject("http")

        const getIgvConfig = function() {

            ApiService.get("/api/genomeConfig")
                .then((response) => {
                    if (response.status === 200) {

                        console.log(response.data)

                        igvInfo = ref({
                            genome: response.data.genomeConfig,
                            tracks: response.data.tracks,
                        })

                        initializeIgv()
                    }

                }).catch((err) => {
                console.log(err)
            })
        }

        const initializeIgv = function () {

            const igvDiv = ref(document.getElementById("igv-div"))

            // the genome information is loaded from API
            const options = {
                genome: igvInfo.value.genome,
            }

            igv.createBrowser(igvDiv.value, options)
                .then(function (browser) {

                    igvBrowser = ref(browser)

                    // initialize tracks by track config loaded from API
                    for (const trackConfig of igvInfo.value.tracks) {
                        browser.loadTrack(trackConfig)
                    }
                })
        }

        return {
            getIgvConfig,
            initializeIgv
        }
    },
    created() {
        this.getIgvConfig()
    },
}

</script>

<style scoped>

</style>