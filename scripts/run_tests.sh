#!/usr/bin/env bash
set -e

TEST_DIR="$(dirname "$0")/../test"
cd "$TEST_DIR"

# Add lib/ to PYTHONPATH for absolute imports to work
export PYTHONPATH="$(cd "$(dirname "$0")/../lib" && pwd):$PYTHONPATH"

if [ $# -eq 0 ]; then
    echo "Running all tests in $TEST_DIR..."
    for t in test_*.py; do
        echo "\n===== Running $t ====="
        python3 "$t"
    done
    echo "\nAll tests completed."
else
    TEST_FILE="test_$1.py"
    if [ -f "$TEST_FILE" ]; then
        echo "Running $TEST_FILE..."
        python3 "$TEST_FILE"
    else
        echo "Test $TEST_FILE not found in $TEST_DIR."
        exit 1
    fi
fi 