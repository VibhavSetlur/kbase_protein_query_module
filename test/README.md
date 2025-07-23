# test/

This directory contains scripts and files needed to test your module's code.

## Structure

- `unit_tests/`: Unit tests for each main module in the toolkit.
- `kbase_protein_network_analysis_toolkit_server_test.py`: Integration and functional tests for the main KBase server implementation.

## Running Tests

1. Ensure your KBase developer token is set in `test_local/test.cfg`.
2. Run all tests with:
   ```bash
   kb-sdk test
   ```

## Documentation

- All test files include clear docstrings and comments for each test case.
- Follow modular, easy-to-understand commenting practices for all new tests.
 
