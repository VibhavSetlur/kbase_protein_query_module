# KBase Protein Query Module - Test Suite

This directory contains comprehensive tests for the KBase Protein Query Module following KBase testing guidelines and best practices.

## Test Structure

```
test/
├── README.md                           # This file
├── test.cfg                            # Test configuration
├── run_tests.py                        # Comprehensive test runner
├── unit_tests/                         # Unit tests
│   ├── test_assign_protein_family_query.py
│   ├── test_check_existence.py
│   ├── test_embedding_generator.py
│   ├── test_network_builder.py
│   ├── test_similarity_index.py
│   ├── test_storage.py
│   ├── test_workflow_orchestrator.py
│   └── test_html_report_generator.py
├── integration_tests/                   # Integration tests
│   └── kbase_protein_query_module_query_server_test.py
├── scripts/                            # Test scripts
│   └── run_pipeline.py                 # Pipeline runner
├── outputs/                            # Test outputs
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

### Quick Start
```bash
# Run all tests
make test-comprehensive

# Run specific test categories
make test-unit
make test-integration
make test-html
make test-kbase
```

### Manual Test Execution
```bash
# Comprehensive test suite
cd test
python3 run_tests.py

# Unit tests only
python3 -m pytest unit_tests/ -v

# Integration tests
python3 -m pytest integration_tests/kbase_protein_query_module_query_server_test.py -v

# HTML report tests
python3 unit_tests/test_html_report_generator.py
```

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
 
