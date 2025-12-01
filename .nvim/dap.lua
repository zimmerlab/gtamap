local dap = require("dap")

-- Configure the Go adapter (delve) - THIS IS REQUIRED
dap.adapters.go = {
	type = "server",
	port = "${port}",
	executable = {
		command = vim.fn.stdpath("data") .. "/mason/bin/dlv",
		args = { "dap", "-l", "127.0.0.1:${port}" },
	},
}

-- Go configuration
-- dap.configurations.go = dap.configurations.go or {}

-- table.insert(dap.configurations.go, {
-- 	type = "go",
-- 	name = "testrunner genome",
-- 	request = "launch",
-- 	target = "debug",
-- 	-- program = "${workspaceFolder}/src/testrunnergenome.go",
-- 	program = "${file}",
-- })

-- table.insert(dap.configurations.go, {
-- 	type = "go",
-- 	name = "test single read",
-- 	request = "launch",
-- 	target = "debug",
-- 	console = "integratedTerminal",
-- 	program = "${workspaceFolder}/src/main.go",
-- 	args = {
-- 		"map",
-- 		"--index",
-- 		"~/Data/kartmann/gtamap/ENSSSCG00000015652.gtai",
-- 		"--output",
-- 		"~/Data/kartmann/gtamap/",
-- 		"--read-origin",
-- 		"dna",
-- 		"--reads-r1",
-- 		"~/Data/kartmann/read2_Vis_36-15.fastq.gz",
-- 	},
-- })

dap.configurations.go = {
	{
		type = "go",
		name = "test single read",
		request = "launch",
		-- target = "debug",
		console = "integratedTerminal",
		program = "${workspaceFolder}/src/main.go",
		args = {
			"map",
			"--index",
			"~/Data/kartmann/gtamap/ENSSSCG00000015652.gtai",
			"--output",
			"~/Data/kartmann/gtamap/",
			"--read-origin",
			"dna",
			"--reads-r1",
			"~/Data/kartmann/read2_Vis_36-15.fastq.gz",
		},
	},
}
