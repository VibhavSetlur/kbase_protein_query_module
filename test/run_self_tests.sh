#!/bin/bash
# Self-test runner for all module components
# Runs main() functions from Python files directly using bash
# Independent of KB SDK - can be run from anywhere

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# Get script directory and project root
script_dir=$(dirname "$(readlink -f "$0")")
project_root=$(dirname "$script_dir")

# Set up Python path
export PYTHONPATH=${PYTHONPATH:-}:$project_root/lib

# Check if we're in a virtual environment, if not, try to activate
if [[ -z "${VIRTUAL_ENV:-}" ]] && [[ -f "$project_root/venv/bin/activate" ]]; then
    echo "Activating virtual environment..."
    source "$project_root/venv/bin/activate"
fi

# Check Python version
PYTHON=${PYTHON:-python3}
if ! command -v $PYTHON &> /dev/null; then
    echo -e "${RED}Error: python3 not found${NC}" >&2
    exit 1
fi

echo -e "${BOLD}${BLUE}============================================================${NC}"
echo -e "${BOLD}${BLUE}Self-Test Runner${NC}"
echo -e "${BOLD}${BLUE}============================================================${NC}"
echo ""

# Track test results
PASSED=0
FAILED=0
TOTAL=0

# Function to run a test as a module
run_test_module() {
    local test_name=$1
    local module_path=$2
    local suppress_warnings=${3:-false}
    
    TOTAL=$((TOTAL + 1))
    
    echo -e "${BLUE}Running ${test_name}...${NC}"
    
    # Run the test as a Python module to handle relative imports correctly
    local output
    local exit_code
    if [ "$suppress_warnings" = "true" ]; then
        # Suppress stderr warnings (HTTP 400, RuntimeWarnings) but capture stdout
        output=$($PYTHON -W ignore -m "$module_path" 2>/dev/null)
        exit_code=$?
        if [ -n "$output" ]; then
            echo "$output"
        fi
    else
        # Show all output but filter out repetitive warnings
        output=$($PYTHON -W ignore::RuntimeWarning -m "$module_path" 2>&1 | grep -v "UniProt metadata query failed: HTTP 400" | grep -v "RuntimeWarning" || true)
        exit_code=${PIPESTATUS[0]}
        echo "$output"
    fi
    
    if [ $exit_code -eq 0 ]; then
        echo -e "${GREEN}✅ ${test_name}: PASSED${NC}"
        PASSED=$((PASSED + 1))
        return 0
    else
        echo -e "${RED}❌ ${test_name}: FAILED (exit code: $exit_code)${NC}"
        FAILED=$((FAILED + 1))
        return $exit_code
    fi
}

# Run all tests as modules
echo -e "${BOLD}Running core component tests...${NC}"
run_test_module "Analysis Config" "kbase_protein_query_module.src.analysis.config"
run_test_module "Analysis Manager" "kbase_protein_query_module.src.analysis.analysis_manager"
run_test_module "Input Manager" "kbase_protein_query_module.src.input.input_manager"
run_test_module "Output Manager" "kbase_protein_query_module.src.output.output_manager"
run_test_module "Workflow Orchestrator" "kbase_protein_query_module.src.core.workflow_orchestrator"

echo ""
echo -e "${BOLD}Running utility component tests...${NC}"
run_test_module "Embedding Generator" "kbase_protein_query_module.src.util.embeddings.generator"
run_test_module "Protein Storage" "kbase_protein_query_module.src.util.storage.storage"

run_test_module "UniProt API" "kbase_protein_query_module.src.util.uniprot.api" "true"

echo ""
echo -e "${BOLD}Running analysis component tests...${NC}"
run_test_module "Network Analysis" "kbase_protein_query_module.src.analysis.network_analysis.network_analysis" "true"

# Print summary
echo ""
echo -e "${BOLD}${BLUE}============================================================${NC}"
echo -e "${BOLD}${BLUE}Test Summary${NC}"
echo -e "${BOLD}${BLUE}============================================================${NC}"
echo ""
echo "Total tests: $TOTAL"
echo -e "${GREEN}✅ Passed: $PASSED${NC}"
if [ $FAILED -gt 0 ]; then
    echo -e "${RED}❌ Failed: $FAILED${NC}"
fi

# Exit with error if any tests failed
if [ $FAILED -gt 0 ]; then
    exit 1
else
    exit 0
fi
