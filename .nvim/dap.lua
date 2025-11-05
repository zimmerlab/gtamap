local dap = require("dap")

-- Go configuration
dap.configurations.go = dap.configurations.go or {}

table.insert(dap.configurations.go, {
	type = "go",
	name = "testrunner genome",
	request = "launch",
	target = "debug",
	-- program = "${workspaceFolder}/src/testrunnergenome.go",
	program = "${file}",
})
