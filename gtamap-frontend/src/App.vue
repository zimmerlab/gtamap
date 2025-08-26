<template>
  <div class="tw:min-h-screen app-background">
    <w-toolbar class="tw:h-10 high-z-index" fixed>
      <div class="tw:flex tw:items-center tw:justify-between tw:w-full">
        <!-- Left section -->
        <div class="tw:flex tw:items-center tw:gap-4">
          <div class="tw:text-md tw:font-bold">
            <router-link to="/"> GTAMap </router-link>
          </div>
          <w-list :items="routes" nav color="primary" class="tw:py-2 w-flex row">
            <template #item="{ item }">
              <span class="tw:mx-2">{{ item.label }}</span>
            </template>
          </w-list>
        </div>

        <!-- Center section -->
        <div class="tw:flex tw:items-center tw:text-sm tw:bg-gray-100 tw:px-4 tw:py-1 tw:rounded">
          <span class="tw:mr-4">Total: {{ progressStore.getReadsTotal }}</span>
          <span class="tw:mr-4 tw:text-green-600">Accepted: {{ progressStore.getReadsAccepted }}</span>
          <span class="tw:text-blue-600">Todo: {{ progressStore.getReadsTodo }}</span>
        </div>

        <!-- Right section -->
        <div class="tw:flex tw:items-center">
          <w-button sm @click="downloadAcceptedBam">
            <w-icon class="mr1">mdi mdi-download</w-icon>
            Download Accepted BAM
          </w-button>
        </div>
      </div>
    </w-toolbar>
    <div class="tw:pt-10">
      <router-view class="grow" v-slot="{ Component }">
        <keep-alive>
          <component :is="Component" />
        </keep-alive>
      </router-view>
    </div>
  </div>
</template>

<script>
import { ref, inject, onMounted } from 'vue'
import { useDataStore } from './store/data.store'
import { useProgressStore } from './store/progress.store'

export default {
  name: 'App',
  setup() {
    const dataStore = useDataStore()
    const progressStore = useProgressStore()

    const DataService = inject('data')
    const ApiService = inject('http')

    const routes = ref([
      {
        label: 'Home',
        route: '/',
      },
      {
        label: 'Test',
        route: '/test',
      },
    ])

    const downloadAcceptedBam = async function () {
      try {
        const response = await ApiService.get('/api/accepted/bam', {
          responseType: 'blob',
        })

        const blob = response.data
        const url = window.URL.createObjectURL(blob)

        const contentDisposition = response.headers['content-disposition']
        let filename = 'accepted_reads.bam'
        if (contentDisposition) {
          const filenameMatch = contentDisposition.match(/filename="(.+)"/)
          if (filenameMatch) {
            filename = filenameMatch[1]
          }
        }

        const link = document.createElement('a')
        link.href = url
        link.download = filename
        document.body.appendChild(link)
        link.click()

        document.body.removeChild(link)
        window.URL.revokeObjectURL(url)
      } catch (error) {
        console.error('Download failed:', error)
      }
    }

    const fetchReads = function () {
      ApiService.get("/api/summary/table")
        .then((response) => {
          if (response.status === 200) {

            // set parents
            for (const item of response.data) {
              item.isSelected = false
              for (const l of item.locations) {
                l.parent = item
              }
            }

            dataStore.setReads(response.data)

            // summaryTableData.value.items = response.data
            //
            // applyFilters()
            //
            // pagination.value.total = summaryTableData.value.filteredItems.length
            // pagination.value.itemsPerPage = 15
            //
            // // add parent information to each location
            // for (const item of summaryTableData.value.items) {
            //   item.isSelected = false
            //   for (const l of item.locations) {
            //     l.parent = item
            //   }
            // }
            //
            // customSort()
          } else {
            console.error('Failed to fetch read summary table data')
          }

          // tableData.value.loading = false
        })
        .catch((err) => {
          console.log(err)
        })
    }

    onMounted(() => {
      DataService.fetchReads()
    })

    return { routes, progressStore, downloadAcceptedBam }
  },
}
</script>

<style scoped>
.high-z-index {
  z-index: 9999 !important;
}

.app-background {
  background-color: #f8fafc;
  background-image: radial-gradient(circle, #e2e8f0 1px, transparent 1px);
  background-size: 20px 20px;
  background-attachment: fixed;
}

/* Alternative: subtle radial gradient */
/* .app-background {
  background: radial-gradient(ellipse at center, #f8fafc 0%, #e2e8f0 100%);
  min-height: 100vh;
} */

/* Alternative: subtle diagonal gradient */
/* .app-background {
  background: linear-gradient(135deg, #f8fafc 0%, #e2e8f0 50%, #f1f5f9 100%);
  min-height: 100vh;
} */
</style>
