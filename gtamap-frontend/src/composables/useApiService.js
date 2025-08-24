import { getCurrentInstance } from 'vue'

export function useApiService() {
	const { appContext } = getCurrentInstance()
	return appContext.config.globalProperties.$http
}
