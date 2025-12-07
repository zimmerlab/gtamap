local dap = require("dap")
-- Go configuration
dap.configurations.go = dap.configurations.go or {}

-- Debug current file
table.insert(dap.configurations.go, {
	type = "go",
	name = "Debug Current File",
	request = "launch",
	mode = "debug",
	program = "${file}",
})

table.insert(dap.configurations.go, {
	type = "go",
	name = "Kartmann",
	request = "launch",
	mode = "debug",
	program = "/home/sam/Projects/gtamap/src/main.go",
	console = "integratedTerminal",
	showLog = true,
	args = {
		"map",
		"--index",
		"/home/sam/Data/kartmann/gtamap/ENSSSCG00000015652.gtai",
		"--output",
		"/home/sam/Data/kartmann/gtamap/",
		"--read-origin",
		"dna",
		"--reads-r1",
		"/home/sam/Data/kartmann/read2_Vis_36-15.fastq.gz",
	},
})

table.insert(dap.configurations.go, {
	type = "go",
	name = "giab chr20",
	request = "launch",
	mode = "debug",
	program = "/home/sam/Projects/gtamap/src/main.go",
	console = "integratedTerminal",
	showLog = true,
	args = {
		"map",
		"--index",
		"/home/sam/Data/gtamap/evaluation/index/20:39659472-39675318.gtai",
		"--output",
		"/home/sam/Data/gtamap/evaluation/output/",
		"--read-origin",
		"dna",
		"--reads-r1",
		"/home/sam/Data/gtamap/evaluation/fastq/test.r1.fastq",
		-- "/home/sam/Data/gtamap/evaluation/fastq/2A2_TGACCA_L001_R1_001.fastq.gz",
		"--reads-r2",
		"/home/sam/Data/gtamap/evaluation/fastq/test.r2.fastq",
		-- "/home/sam/Data/gtamap/evaluation/fastq/2A2_TGACCA_L001_R2_001.fastq.gz",
	},
})
