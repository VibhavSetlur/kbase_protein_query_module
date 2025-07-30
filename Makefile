SERVICE = kbase_protein_query_module
SERVICE_CAPS = kbase_protein_query_module
SPEC_FILE = kbase_protein_query_module.spec
URL = https://kbase.us/services/kbase_protein_query_module
DIR = $(shell pwd)
LIB_DIR = lib
SCRIPTS_DIR = scripts
TEST_DIR = test
TEST_LOCAL_DIR = test_local
LBIN_DIR = bin
WORK_DIR = /kb/module/work/tmp
EXECUTABLE_SCRIPT_NAME = run_$(SERVICE_CAPS)_async_job.sh
STARTUP_SCRIPT_NAME = start_server.sh
TEST_SCRIPT_NAME = run_tests.sh

.PHONY: test

default: compile

all: compile build build-startup-script build-executable-script build-test-script

compile:
	kb-sdk compile $(SPEC_FILE) \
		--out $(LIB_DIR) \
		--pysrvname $(SERVICE_CAPS).$(SERVICE_CAPS)Server \
		--pyimplname $(SERVICE_CAPS).$(SERVICE_CAPS)Impl;

build:
	chmod +x $(SCRIPTS_DIR)/entrypoint.sh

build-executable-script:
	mkdir -p $(LBIN_DIR)
	echo '#!/bin/bash' > $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
	echo 'script_dir=$$(dirname "$$(readlink -f "$$0")")' >> $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
	echo 'export PYTHONPATH=$$script_dir/../$(LIB_DIR):$$PATH:$$PYTHONPATH' >> $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
	echo 'python3 -u $$script_dir/../$(LIB_DIR)/$(SERVICE_CAPS)/$(SERVICE_CAPS)Server.py $$1 $$2 $$3' >> $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)
	chmod +x $(LBIN_DIR)/$(EXECUTABLE_SCRIPT_NAME)

build-startup-script:
	mkdir -p $(LBIN_DIR)
	echo '#!/bin/bash' > $(SCRIPTS_DIR)/$(STARTUP_SCRIPT_NAME)
	echo 'script_dir=$$(dirname "$$(readlink -f "$$0")")' >> $(SCRIPTS_DIR)/$(STARTUP_SCRIPT_NAME)
	echo 'export KB_DEPLOYMENT_CONFIG=$$script_dir/../deploy.cfg' >> $(SCRIPTS_DIR)/$(STARTUP_SCRIPT_NAME)
	echo 'export PYTHONPATH=$$script_dir/../$(LIB_DIR):$$PATH:$$PYTHONPATH' >> $(SCRIPTS_DIR)/$(STARTUP_SCRIPT_NAME)
	echo 'uwsgi --master --processes 5 --threads 5 --http :5000 --wsgi-file $$script_dir/../$(LIB_DIR)/$(SERVICE_CAPS)/$(SERVICE_CAPS)Server.py' >> $(SCRIPTS_DIR)/$(STARTUP_SCRIPT_NAME)
	chmod +x $(SCRIPTS_DIR)/$(STARTUP_SCRIPT_NAME)

build-test-script:
	echo '#!/bin/bash' > $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	echo 'script_dir=$$(dirname "$$(readlink -f "$$0")")' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	echo 'export KB_DEPLOYMENT_CONFIG=$$script_dir/../deploy.cfg' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	echo 'export KB_AUTH_TOKEN=`cat /kb/module/work/token`' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	echo 'echo "Removing temp files..."' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	echo 'rm -rf $(WORK_DIR)/*' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	echo 'echo "...done removing temp files."' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	echo 'export PYTHONPATH=$$script_dir/../$(LIB_DIR):$$PATH:$$PYTHONPATH' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	echo 'cd $$script_dir/../$(TEST_DIR)' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	echo 'python3 -m nose --with-coverage --cover-package=$(SERVICE_CAPS) --cover-html --cover-html-dir=/kb/module/work/test_coverage --nocapture  --nologcapture .' >> $(TEST_DIR)/$(TEST_SCRIPT_NAME)
	chmod +x $(TEST_DIR)/$(TEST_SCRIPT_NAME)

test:
	if [ ! -f /kb/module/work/token ]; then echo -e '\nOutside a docker container please run "kb-sdk test" rather than "make test"\n' && exit 1; fi
	bash $(TEST_DIR)/$(TEST_SCRIPT_NAME)

# KBase SDK test command (proper way to test KBase modules)
test-kbase: compile
	@echo "Running KBase SDK tests..."
	@kb-sdk test

# Generate test data
generate-test-data:
	@echo "Generating test data..."
	@cd $(TEST_LOCAL_DIR) && python3 generate_test_data.py
	@echo "Test data generation completed."

# Run tests with proper data generation
test-with-data: generate-test-data test-kbase

clean:
	rm -rfv $(LBIN_DIR)

clean-test:
	@echo "Cleaning test files..."
	@rm -rf $(TEST_LOCAL_DIR)/data
	@echo "✅ Test files cleaned"

# Validate the module
validate:
	@echo "Validating KBase module..."
	@kb-sdk validate

# Install the module
install:
	@echo "Installing KBase module..."
	@kb-sdk install

# Help target
help:
	@echo "Available targets:"
	@echo "  compile          - Compile the KIDL specification"
	@echo "  build            - Build all scripts"
	@echo "  test             - Run tests inside Docker container"
	@echo "  test-kbase       - Run KBase SDK tests"
	@echo "  test-with-data   - Generate test data and run tests"
	@echo "  generate-test-data - Generate test data only"
	@echo "  validate         - Validate the module"
	@echo "  install          - Install the module"
	@echo "  clean            - Clean build artifacts"
	@echo "  clean-test       - Clean test data"
	@echo "  help             - Show this help"
