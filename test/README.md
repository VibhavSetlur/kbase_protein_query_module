# KBase Protein Query Module - Test Suite

This directory contains comprehensive tests for the KBase Protein Query Module, organized by test type and functionality.

## Test Structure

```
test/
├── README.md                           # This file
├── requirements-test.txt               # Testing dependencies
├── pytest.ini                         # Pytest configuration
├── run_comprehensive_tests.py         # Main test runner script
├── conftest.py                        # Shared test fixtures
├── kbase_protein_query_module_server_test.py  # Original KBase server tests
├── unit_tests/                        # Unit tests by module
│   ├── conftest.py                    # Unit test fixtures
│   ├── test_runner.py                 # Unit test runner
│   ├── input/                         # Input module tests
│   │   ├── test_input_manager.py
│   │   ├── test_protein_sequence_processor.py
│   │   ├── test_uniprot_processor.py
│   │   └── test_workspace_processor.py
│   ├── analysis/                      # Analysis module tests
│   │   ├── test_analysis_manager.py
│   │   └── test_analysis_config.py
│   ├── output/                        # Output module tests
│   │   └── (to be implemented)
│   ├── core/                          # Core module tests
│   │   ├── test_pipeline_config.py
│   │   └── test_workflow_orchestrator.py
│   ├── util/                          # Utility module tests
│   │   ├── test_protein_embedding_generator.py
│   │   └── (additional utility tests)
│   ├── integration/                   # Integration tests
│   │   ├── test_input_analysis_integration.py
│   │   └── (additional integration tests)
│   ├── functional/                    # Functional tests
│   │   ├── test_end_to_end_workflow.py
│   │   └── (additional functional tests)
│   └── kbase_compliance/              # KBase compliance tests
│       └── test_kbase_standards.py
└── test_outputs/                      # Test output directory
```

## Test Categories

### Unit Tests
- **Input Tests**: Test input processing modules (protein sequences, UniProt IDs, workspace objects)
- **Analysis Tests**: Test analysis management and configuration
- **Output Tests**: Test output generation and management
- **Core Tests**: Test core workflow orchestration and configuration
- **Utility Tests**: Test utility modules (embeddings, storage, similarity search)

### Integration Tests
- **Input-Analysis Integration**: Test complete flow from input processing to analysis execution
- **Analysis-Output Integration**: Test analysis result processing and output generation
- **Cross-Module Integration**: Test interactions between different modules

### Functional Tests
- **End-to-End Workflows**: Test complete workflows from input to output
- **KBase Integration**: Test integration with KBase services and standards
- **Performance Tests**: Test performance and scalability
- **Error Handling**: Test error scenarios and recovery

### KBase Compliance Tests
- **Standards Compliance**: Test adherence to KBase development standards
- **API Compliance**: Test API structure and behavior compliance
- **Integration Compliance**: Test KBase service integration compliance

## Running Tests

### Quick Start

```bash
# Install test dependencies
pip install -r requirements-test.txt

# Run all tests
python run_comprehensive_tests.py

# Run with verbose output
python run_comprehensive_tests.py --verbose

# Run specific test categories
python run_comprehensive_tests.py --categories unit_tests integration_tests

# Run only KBase server tests
python run_comprehensive_tests.py --server-tests

# Run linting tests
python run_comprehensive_tests.py --linting
```

### Using pytest directly

```bash
# Run all unit tests
pytest unit_tests/

# Run specific test category
pytest unit_tests/input/

# Run specific test file
pytest unit_tests/input/test_input_manager.py

# Run with coverage
pytest --cov=../../lib unit_tests/

# Run with verbose output
pytest -v unit_tests/

# Run specific test method
pytest unit_tests/input/test_input_manager.py::TestInputManager::test_initialization
```

### Using unittest

```bash
# Run KBase server tests
python -m unittest discover -s . -p "*server_test.py"

# Run specific test module
python -m unittest unit_tests.input.test_input_manager
```

## Test Configuration

### Environment Variables

```bash
# KBase authentication (for integration tests)
export KB_AUTH_TOKEN="your_kbase_token"
export KB_WORKSPACE_URL="https://appdev.kbase.us/services/ws"

# Test configuration
export TEST_OUTPUT_DIR="/tmp/test_outputs"
export TEST_VERBOSE="true"
```

### Pytest Configuration

The `pytest.ini` file configures:
- Test discovery patterns
- Output formatting
- Markers for test categorization
- Timeout settings
- Warning filters

### Test Fixtures

The `conftest.py` files provide shared fixtures for:
- Mock KBase clients
- Sample test data
- Temporary directories
- Configuration objects

## Writing Tests

### Test Structure

```python
import pytest
from unittest.mock import Mock, patch
from your_module import YourClass

class TestYourClass:
    """Test cases for YourClass."""
    
    def test_initialization(self, test_config):
        """Test class initialization."""
        instance = YourClass(test_config)
        assert instance.config == test_config
    
    def test_method_with_mock(self, mock_kb_util):
        """Test method with mocked dependencies."""
        instance = YourClass({}, mock_kb_util)
        result = instance.some_method()
        assert result is not None
```

### Test Markers

Use pytest markers to categorize tests:

```python
@pytest.mark.unit
def test_unit_functionality():
    """Unit test for specific functionality."""
    pass

@pytest.mark.integration
def test_integration_flow():
    """Integration test for component interactions."""
    pass

@pytest.mark.slow
def test_performance():
    """Performance test that takes time."""
    pass
```

### Mock Usage

```python
# Mock external dependencies
@patch('your_module.ExternalService')
def test_with_external_mock(mock_service):
    mock_service.return_value.get_data.return_value = "mocked_data"
    result = your_function()
    assert result == "expected_result"

# Use fixtures for common mocks
def test_with_fixture_mock(mock_workspace_client):
    result = process_workspace_data(mock_workspace_client)
    assert result is not None
```

## Test Data

### Sample Data Fixtures

Tests use fixtures for consistent sample data:
- `sample_protein_sequences`: Protein sequence strings
- `sample_uniprot_ids`: UniProt protein identifiers
- `sample_embeddings`: Protein embedding arrays
- `sample_metadata_df`: Protein metadata DataFrames

### Mock Data

Tests use mocks for:
- KBase workspace clients
- KBase utility libraries
- External API services
- File system operations

## Continuous Integration

### GitHub Actions

```yaml
name: Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v3
    - name: Set up Python
      uses: actions/setup-python@v4
      with:
        python-version: '3.8'
    - name: Install dependencies
      run: |
        pip install -r test/requirements-test.txt
    - name: Run tests
      run: |
        python test/run_comprehensive_tests.py --verbose
```

### Local CI

```bash
#!/bin/bash
# Run all tests locally
python test/run_comprehensive_tests.py --verbose --failfast
```

## Troubleshooting

### Common Issues

1. **Import Errors**: Ensure the lib directory is in the Python path
2. **Mock Failures**: Check that mocks are properly configured
3. **Timeout Issues**: Increase timeout in pytest.ini
4. **Permission Errors**: Check file permissions for test output directories

### Debug Mode

```bash
# Run with debug output
pytest -v -s unit_tests/input/test_input_manager.py

# Run single test with debugging
pytest -v -s unit_tests/input/test_input_manager.py::TestInputManager::test_initialization
```

### Test Reports

Test reports are generated in:
- `test_outputs/comprehensive_test_report.txt` (human-readable)
- `test_outputs/comprehensive_test_report.json` (machine-readable)

## Contributing

When adding new tests:

1. Follow the existing test structure and naming conventions
2. Use appropriate test markers for categorization
3. Include comprehensive docstrings
4. Mock external dependencies appropriately
5. Add fixtures for reusable test data
6. Update this README if adding new test categories

## KBase Compliance

Tests ensure compliance with:
- KBase API standards
- KBase service integration patterns
- KBase error handling requirements
- KBase performance standards
- KBase security requirements

Run compliance tests with:
```bash
python test/run_comprehensive_tests.py --categories kbase_compliance
```

