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
	type = "delve",
	name = "debug simple",
	request = "launch",
	mode = "debug",
	program = "/home/sam/Projects/gtamap/src/runner.go",
})

table.insert(dap.configurations.go, {
	type = "go",
	name = "Kartmann",
	request = "launch",
	mode = "debug",
	program = "/home/sam/Projects/gtamap/src/main.go",
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
	type = "delve",
	name = "debugging gap",
	request = "launch",
	-- mode = "debug",
	program = "/home/sam/Projects/gtamap/src/main.go",
	outputMode = "remote",
	args = {
		"map",
		"--config",
		"/home/sam/Data/gtamap/evaluation/bugs/20251211_gaptest/gtamap_config.yaml",
		"--index",
		"/home/sam/Data/gtamap/evaluation/bugs/20251211_gaptest/20:18931791-18934204.gtai",
		"--output",
		"/home/sam/Data/gtamap/evaluation/bugs/20251211_gaptest/",
		"--read-origin",
		"dna",
		"--reads-r1",
		"/home/sam/Data/gtamap/evaluation/bugs/20251211_gaptest/test.r1.fastq",
		"--reads-r2",
		"/home/sam/Data/gtamap/evaluation/bugs/20251211_gaptest/test.r2.fastq",
	},
})

table.insert(dap.configurations.go, {
	type = "delve",
	name = "debugging pig",
	request = "launch",
	-- mode = "debug",
	program = "/home/users/klein/Projects/gtamap/src/main.go",
	outputMode = "remote",
	args = {
		"map",
		"--config",
		"/mnt/proj/projekte/expressionlab/dmd-pig/gtamap/gtamap_config.yaml",
		"--index",
		"/mnt/proj/projekte/expressionlab/dmd-pig/results/gtamap/78918581/dc092a89/ENSSSCG00000053426.gtai",
		"--output",
		"/home/users/klein/test/",
		"--threads",
		"1",
		"--read-origin",
		"rna",
		"--reads-r1",
		"/mnt/raidinput/input/own/SeqReads/Jaudas/01.RawData/Tri_Bi_14085_1.fq.gz",
		"--reads-r2",
		"/mnt/raidinput/input/own/SeqReads/Jaudas/01.RawData/Tri_Bi_14085_2.fq.gz",
	},
})

table.insert(dap.configurations.go, {
	type = "delve",
	name = "debugging weight",
	request = "launch",
	-- mode = "debug",
	program = "/home/sam/Projects/gtamap/src/main.go",
	outputMode = "remote",
	args = {
		"map",
		"--config",
		"/home/sam/Data/gtamap/evaluation/genes/dmd/gtamap_config.yaml",
		"--index",
		"/home/sam/Data/gtamap/evaluation/genes/dmd/index.kmer-10.v2.gtai",
		"--output",
		"/home/sam/Data/gtamap/evaluation/genes/dmd/",
		"--read-origin",
		"rna",
		"--reads-r1",
		"/home/sam/Data/gtamap/evaluation/genes/dmd/1916600.r1.fastq",
		"--reads-r2",
		"/home/sam/Data/gtamap/evaluation/genes/dmd/1916600.r2.fastq",
	},
})

table.insert(dap.configurations.go, {
	type = "go",
	name = "rs bug",
	request = "launch",
	mode = "debug",
	program = "/home/sam/Projects/gtamap/src/main.go",
	showLog = true,
	args = {
		"map",
		"--config",
		"/home/sam/Data/gtamap/evaluation/bugs/20251211_rs200630371/gtamap_config.yaml",
		"--index",
		"/home/sam/Data/gtamap/evaluation/bugs/20251211_rs200630371/20:18345927-18347017.gtai",
		"--output",
		"/home/sam/Data/gtamap/evaluation/bugs/",
		"--read-origin",
		"dna",
		"--reads-r1",
		"/home/sam/Data/gtamap/evaluation/bugs/20251211_rs200630371/rs200630371.novoalign.R1.fastq",
		"--reads-r2",
		"/home/sam/Data/gtamap/evaluation/bugs/20251211_rs200630371/rs200630371.novoalign.R2.fastq",
	},
})
