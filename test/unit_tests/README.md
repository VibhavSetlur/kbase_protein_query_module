# Unit Tests for KBase Protein Query Module

This directory contains comprehensive unit tests for all components of the KBase Protein Query Module.

## Test Organization

### Core Module Tests
- `test_input_parser.py` - Tests for input parsing functionality
- `test_embedding_generator.py` - Tests for protein embedding generation
- `test_similarity_index.py` - Tests for similarity search and indexing
- `test_assign_protein_family.py` - Tests for protein family assignment
- `test_network_builder.py` - Tests for network construction and analysis
- `test_sequence_analyzer.py` - Tests for sequence analysis features
- `test_storage.py` - Tests for data storage and retrieval
- `test_html_report_generator.py` - Tests for HTML report generation

### Workflow Tests
- `test_workflow_orchestrator.py` - Tests for workflow orchestration
- `test_unified_workflow.py` - Tests for unified workflow orchestrator
- `test_workflow_integration.py` - Tests for complete workflow integration
- `test_stages.py` - Tests for all workflow stages (input, processing, analysis, output)

### Query Type Tests
- `test_query_types.py` - Tests for all supported query types:
  - Protein existence queries
  - Similarity search queries
  - Family assignment queries
  - Network analysis queries

### Edge Case Tests
- `test_check_existence.py` - Tests for protein existence checking
- Additional edge case tests in each module

## Test Coverage

### Input Variations
- FASTA files with different header formats
- CSV files with various column configurations
- JSON files with different structures
- Sequence strings (single and multiple)
- Protein ID lists
- Workspace objects
- Batch processing of multiple input sources

### Query Types
- Single protein queries
- Multiple protein queries
- Sequence-based queries
- Similarity threshold queries
- Filtered queries with metadata
- Batch queries

### Configuration Options
- Different embedding models (ESM2 variants)
- Various pooling methods
- Network construction methods
- Similarity thresholds
- K-neighbors parameters

### Edge Cases
- Empty files and inputs
- Invalid sequences and formats
- Very short and very long sequences
- Missing data scenarios
- Error handling and recovery

### Performance Tests
- Large dataset processing
- Memory usage optimization
- Execution time validation
- Caching mechanisms

## Running Tests

### Individual Test Files
```bash
# Run specific test file
python -m pytest test_input_parser.py -v

# Run with coverage
python -m pytest test_input_parser.py --cov=kbase_protein_query_module.src.input_parser
```

### All Unit Tests
```bash
# Run all unit tests
python -m pytest test/unit_tests/ -v

# Run with coverage report
python -m pytest test/unit_tests/ --cov=kbase_protein_query_module.src --cov-report=html
```

### Comprehensive Test Suite
```bash
# Run comprehensive test suite (includes all test types)
python test/run_comprehensive_tests.py
```

## Test Data

Test data is generated dynamically for each test to ensure:
- Isolation between tests
- No dependency on external files
- Consistent test results
- Easy cleanup

## Mock Usage

Tests use mocks for:
- External API calls
- File system operations
- Network requests
- Heavy computational tasks (embeddings, similarity calculations)

This ensures:
- Fast test execution
- Reliable test results
- No external dependencies
- Focused unit testing

## Test Results

Test results are stored in:
- `test_report.txt` - Basic test results
- `comprehensive_test_report.txt` - Detailed test report
- `comprehensive_test.log` - Test execution log

## Continuous Integration

These tests are designed to run in CI/CD pipelines:
- No external dependencies
- Fast execution
- Clear pass/fail criteria
- Comprehensive coverage reporting

## Best Practices

1. **Test Isolation**: Each test is independent and cleans up after itself
2. **Descriptive Names**: Test methods have clear, descriptive names
3. **Comprehensive Coverage**: Tests cover normal cases, edge cases, and error conditions
4. **Mock External Dependencies**: External services and heavy computations are mocked
5. **Clear Assertions**: Each test has clear, specific assertions
6. **Documentation**: Tests include docstrings explaining what they test
7. **Performance**: Tests are optimized for speed while maintaining coverage 