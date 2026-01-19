.PHONY: build clean install check-go

VERSION=$(shell grep 'toolVersion string' src/config/config.go | sed 's/.*"\(.*\)".*/\1/')
BINARY_NAME=gtamap_v$(VERSION)
SRC_DIR=./src
BUILD_DIR=./build

check-go:
	@which go > /dev/null 2>&1 || (echo "Error: Go is not installed or not in PATH" && exit 1)
	@echo "Go found: $$(go version)"

build: check-go
	@mkdir -p $(BUILD_DIR)
	go build -o $(BUILD_DIR)/$(BINARY_NAME) $(SRC_DIR)/main.go
	@echo "Built $(BUILD_DIR)/$(BINARY_NAME)"

clean:
	rm -rf $(BUILD_DIR)

install: build
	cp $(BUILD_DIR)/$(BINARY_NAME) /usr/local/bin/
	@echo "Installed $(BINARY_NAME) to /usr/local/bin/"
