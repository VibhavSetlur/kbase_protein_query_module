# Testing Guide

## Test Structure

```
test/
├── unit_tests/                    # Unit tests for backend components
│   ├── core/                     # Core component tests
│   ├── analysis/                 # Analysis module tests
│   ├── input/                    # Input processing tests
│   ├── outputs/                  # Output generation tests
│   └── util/                     # Utility tests
├── integration_tests/            # Integration tests
└── kbase_protein_query_module_server_test.py  # KBase integration tests
```

## Running Tests



## KBase Integration Testing

### Using kb-sdk

Run server tests through the KBase SDK to align with KBase documentation and conventions.

```bash
# Run full KBase validation and server tests
kb-sdk test

# Run specific steps as needed
kb-sdk validate
kb-sdk build
```

Environment notes:
- Ensure `KB_AUTH_TOKEN`, `KB_WORKSPACE_URL`, and `SDK_CALLBACK_URL` are set when running in environments that require them.
- Tests mock KBase clients (`Workspace`, `KBaseReport`, `DataFileUtil`) to keep execution deterministic.

## Test Types

### Unit Tests

Test individual components in isolation:

```python
def test_protein_storage_add_protein():
    storage = ProteinStorage()
    protein = {"id": "P12345", "sequence": "MKFLVNVALVFMVVYISYIYA"}
    
    result = storage.add_protein(protein)
    
    assert result is True
    assert storage.get_protein("P12345") == protein
```

### Integration Tests

Test component interactions:

```python
def test_workflow_orchestrator_integration():
    orchestrator = WorkflowOrchestrator()
    input_data = {"input_type": "uniprot_ids", "input_data": ["P12345"]}
    
    result = orchestrator.run_workflow(input_data, "/tmp/test")
    
    assert result.success is True
    assert len(result.analyses_completed) > 0
```

### Mock Testing

Use mocks for external dependencies:

```python
@patch('requests.get')
def test_protein_existence_checker(mock_get):
    mock_get.return_value.json.return_value = {"entry": {"id": "P12345"}}
    
    checker = ProteinExistenceChecker()
    exists, metadata = checker.check_protein_exists("P12345")
    
    assert exists is True
    assert metadata["id"] == "P12345"
```

## Test Configuration

### pytest.ini

```ini
[tool:pytest]
testpaths = test
python_files = test_*.py
python_classes = Test*
python_functions = test_*
addopts = -v --tb=short
markers =
    unit: Unit tests
    integration: Integration tests
    slow: Slow running tests
```

### conftest.py

```python
import pytest
import tempfile

@pytest.fixture
def temp_dir():
    """Create temporary directory for tests."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir

@pytest.fixture
def sample_proteins():
    """Sample protein data for tests."""
    return [
        {"id": "P12345", "sequence": "MKFLVNVALVFMVVYISYIYA"},
        {"id": "P67890", "sequence": "MKFLVNVALVFMVVYISYIYA"}
    ]
```

## Writing Tests

### Test Structure

```python
class TestYourComponent:
    """Test cases for YourComponent."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.component = YourComponent()
        self.sample_data = {"key": "value"}
    
    def test_initialization(self):
        """Test component initialization."""
        assert self.component is not None
        assert self.component.name == "expected_name"
    
    def test_method_success(self):
        """Test successful method execution."""
        result = self.component.method(self.sample_data)
        
        assert result is not None
        assert result["status"] == "success"
    
    def test_method_failure(self):
        """Test method failure handling."""
        with pytest.raises(ValueError):
            self.component.method(None)
```

### Test Best Practices

1. **Clear Names**: Use descriptive test names
2. **Single Responsibility**: Test one thing per test
3. **Arrange-Act-Assert**: Structure tests clearly
4. **Mock Dependencies**: Mock external services
5. **Edge Cases**: Test boundary conditions

### Assertions

```python
# Basic assertions
assert result is not None
assert result["status"] == "success"
assert len(result["data"]) > 0

# Exception testing
with pytest.raises(ValueError):
    component.method(invalid_input)

# Approximate equality
assert result["score"] == pytest.approx(0.95, rel=1e-2)
```

## Debugging Tests

### Verbose Output

```bash
# Maximum verbosity
pytest -vvv

# Show local variables on failure
pytest -l

# Show print statements
pytest -s
```

### Debug Mode

```python
# Add breakpoints
import pdb; pdb.set_trace()

# Use pytest debugging
pytest --pdb
```

## Performance Testing

### Benchmark Tests

```python
import time

def test_performance():
    """Test component performance."""
    start_time = time.time()
    
    result = component.expensive_operation()
    
    end_time = time.time()
    execution_time = end_time - start_time
    
    assert execution_time < 5.0  # Should complete in under 5 seconds
    assert result is not None
```

## Best Practices

### Test Organization

1. **Group Related Tests**: Use test classes for related functionality
2. **Clear Naming**: Use descriptive test and method names
3. **Single Responsibility**: Each test should verify one thing
4. **Independent Tests**: Tests should not depend on each other

### Test Maintenance

1. **Keep Tests Updated**: Update tests when code changes
2. **Remove Obsolete Tests**: Delete tests for removed functionality
3. **Refactor Tests**: Improve test quality over time
4. **Document Tests**: Add comments for complex test logic

### Performance

1. **Fast Tests**: Keep unit tests fast
2. **Parallel Execution**: Use parallel test execution
3. **Selective Testing**: Run only relevant tests during development
4. **CI Optimization**: Optimize tests for continuous integration