# KBase Protein Query Module - Unified Test Suite

This directory contains a unified test suite for the KBase Protein Query Module with consolidated outputs and a single test runner.

## Test Structure

```
test/
├── README.md                           # This file
├── test.cfg                            # Test configuration
├── run_tests.py                        # Unified test runner (single entry point)
├── clean_and_run.py                    # Clean and run script
├── unit_tests/                         # Unit tests
│   ├── core/                           # Core module tests
│   ├── data/                           # Data model tests
│   ├── processing/                     # Processing pipeline tests
│   ├── analysis/                       # Analysis component tests
│   ├── storage/                        # Storage layer tests
│   ├── reports/                        # Report generation tests
│   ├── stages/                         # Pipeline stage tests
│   ├── utils/                          # Utility function tests
│   └── workflows/                      # Workflow orchestration tests
├── integration_tests/                   # Integration tests
│   └── kbase_protein_query_module_query_server_test.py
├── scripts/                            # Test scripts
│   └── run_pipeline.py                 # Pipeline runner
├── outputs/                            # Consolidated test outputs
│   ├── html_reports/                   # HTML reports
│   ├── exports/                        # Data exports
│   ├── reports/                        # Test reports
│   ├── data/                           # Test data outputs
│   └── *.log                           # Log files
└── data/                               # Test data
```

## Test Categories

### 1. Unit Tests (`unit_tests/`)
- **Core functionality testing** for each module
- **Isolated testing** without external dependencies
- **Mock data** for consistent results
- **Fast execution** for development feedback

### 2. Integration Tests (`integration_tests/`)
- **Workspace integration** testing
- **KBase service** functionality
- **Real data** testing with actual protein data
- **End-to-end** workflow testing

### 3. HTML Report Tests (`unit_tests/test_html_report_generator.py`)
- **Report generation** functionality
- **Sequence analysis** integration
- **Visualization** components
- **KBase report** integration

### 4. KBase SDK Tests
- **SDK compliance** testing
- **Deployment** validation
- **Service registration** testing

## Running Tests

### Quick Start (Recommended)
```bash
# Run unified test suite with cleanup
cd test
python3 clean_and_run.py

# Or run the unified test runner directly
python3 run_tests.py
```

### Manual Test Execution
```bash
# Unified test suite (recommended)
cd test
python3 run_tests.py

# Clean and run (includes cleanup of old outputs)
python3 clean_and_run.py

# Individual test categories (if needed)
python3 -m pytest unit_tests/ -v
python3 -m pytest integration_tests/ -v
```

### Output Consolidation
All test outputs are now consolidated in the `test/outputs/` directory:
- **HTML Reports**: `outputs/html_reports/`
- **Data Exports**: `outputs/exports/`
- **Test Reports**: `outputs/reports/`
- **Log Files**: `outputs/*.log`
- **Test Results**: `outputs/test_results.json`

### KBase SDK Testing
```bash
# Standard KBase SDK test
kb-sdk test

# Quick validation
kb-sdk test --quick
```

## Test Configuration

The `test.cfg` file contains configuration for:
- **Data paths** and directories
- **Test parameters** (protein IDs, sequences, etc.)
- **Model configuration** (ESM-2 parameters)
- **Network parameters** (similarity thresholds)
- **Timeout settings** for long-running tests
- **Error handling** configuration

## Test Data Requirements

### Required Data Structure
```
data/
├── families/                    # Protein family data
├── indexes/                     # FAISS indexes
├── metadata/                    # Protein metadata
└── family_centroids/           # Family centroid data
    └── family_centroids_binary.npz
```

### Data Generation
```bash
# Generate test data with automatic index creation
make generate-test-data

# Create indexes for existing data
make create-indexes
```

## Test Reports

### Generated Reports
- **test_report.txt**: Comprehensive test results
- **test.log**: Detailed test logs
- **Coverage reports**: Code coverage analysis

### Report Contents
- **Test summary** with pass/fail statistics
- **Performance metrics** and timing
- **Data availability** status
- **Error details** and recommendations
- **KBase compliance** validation

## Error Handling

### Test Failures
1. **Check data availability** - Ensure required data files exist
2. **Verify dependencies** - Check Python packages and KBase SDK
3. **Review logs** - Check `test.log` for detailed error information
4. **Environment issues** - Verify KBase environment setup

### Common Issues
- **Missing data files**: Run `make generate-test-data`
- **SDK issues**: Run `kb-sdk validate`
- **Import errors**: Check `PYTHONPATH` and dependencies
- **Workspace errors**: Verify KBase authentication and workspace access

## KBase Testing Guidelines

### Best Practices
- **Comprehensive coverage** of all module functions
- **Proper error handling** and graceful degradation
- **Real data testing** with actual protein sequences
- **Performance testing** for large datasets
- **Documentation** of test cases and expected results

### Quality Standards
- **80%+ success rate** for production readiness
- **Complete integration** with KBase services
- **Proper logging** and error reporting
- **Consistent results** across different environments

## Development Workflow

### Adding New Tests
1. **Create test file** in appropriate directory
2. **Follow naming convention** `test_*.py`
3. **Include proper documentation** and docstrings
4. **Add to test runner** if needed
5. **Update configuration** if new parameters required

### Test Maintenance
- **Regular execution** of full test suite
- **Update test data** when module changes
- **Review coverage** reports for gaps
- **Validate KBase compliance** regularly

## References

- [KBase SDK Documentation](https://kbase.github.io/kb_sdk_docs/)
- [KBase Testing Guidelines](https://kbase.github.io/kb_sdk_docs/tutorial/7_implement.html)
- [Python Testing Best Practices](https://docs.python.org/3/library/unittest.html)
- [Module Documentation](https://github.com/kbaseapps/kbase_protein_query_module)

## Support

For issues with tests:
1. Check the test logs in `test.log`
2. Review the test report in `test_report.txt`
3. Verify data availability and configuration
4. Contact the development team with specific error details

---

**Author**: Vibhav Setlur  
**Contact**: https://kbase.us/contact-us/  
**Last Updated**: 2024
 
