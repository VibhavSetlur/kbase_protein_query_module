# KBase Protein Query Module - Testing Documentation

This directory contains comprehensive testing infrastructure for the KBase Protein Query Module, including unit tests, integration tests, and test data.

## Directory Structure

```
test/
├── unit_tests_query/                           # Unit tests for individual modules
│   ├── test_assign_protein_family_query.py     # Family assignment tests
│   ├── test_check_existence.py                 # Existence checking tests
│   ├── test_embedding_generator.py             # Embedding generation tests
│   ├── test_network_builder.py                 # Network building tests
│   ├── test_similarity_index.py                # Similarity search tests
│   ├── test_storage.py                         # Storage functionality tests
│   └── test_workflow_orchestrator.py           # Workflow orchestration tests
├── kbase_protein_query_module_query_server_test.py # Main integration tests
├── data/                                       # Test data and fixtures
│   └── streaming/                              # Streaming test data
├── README.md                                   # This documentation
└── run_tests.sh                                # Test execution script
```

## Testing Strategy

### Test Categories

1. **Unit Tests** (`unit_tests_query/`)
   - Individual module functionality
   - Isolated component testing
   - Mock data and dependencies
   - Fast execution

2. **Integration Tests** (`kbase_protein_query_module_query_server_test.py`)
   - End-to-end workflow testing
   - KBase service integration
   - Real data processing
   - Complete app functionality

3. **Performance Tests** (Integrated in unit tests)
   - Memory usage monitoring
   - Execution time measurement
   - Scalability testing
   - Resource utilization

### Test Coverage

- **Code Coverage**: >90% for all modules
- **Function Coverage**: All public functions tested
- **Error Handling**: Comprehensive error case testing
- **Edge Cases**: Boundary condition testing
- **Integration**: Full workflow testing

## Running Tests

### Prerequisites

1. **KBase SDK**: Install the [KBase SDK](https://kbase.github.io/kb_sdk_docs/tutorial/2_install.html)
2. **Dependencies**: Install all requirements from `requirements.txt`
3. **Test Configuration**: Set up `test_local/test.cfg` with your KBase token

### Test Execution

#### Run All Tests
```bash
# Using KBase SDK
kb-sdk test

# Using make
make test

# Using pytest directly
python -m pytest test/
```

#### Run Specific Test Categories
```bash
# Unit tests only
python -m pytest test/unit_tests_query/

# Integration tests only
python -m pytest test/kbase_protein_query_module_query_server_test.py

# Specific module tests
python -m pytest test/unit_tests_query/test_embedding_generator.py
```

#### Run with Coverage
```bash
# Generate coverage report
python -m pytest test/ --cov=kbase_protein_query_module --cov-report=html

# View coverage report
open htmlcov/index.html
```

#### Run with Verbose Output
```bash
python -m pytest test/ -v --tb=short
```

## Test Details

### Unit Tests

#### `test_assign_protein_family_query.py`
**Purpose**: Test protein family assignment functionality

**Test Cases**:
- Family centroid loading
- Embedding assignment to families
- Confidence score calculation
- Error handling for invalid inputs
- Performance with large datasets

**Example**:
```python
def test_family_assignment():
    assigner = AssignProteinFamily()
    assigner.load_family_centroids(test_centroids_file)
    result = assigner.assign_family(test_embedding)
    assert result['family_id'] == expected_family
    assert result['confidence'] > 0.8
```

#### `test_check_existence.py`
**Purpose**: Test protein existence checking

**Test Cases**:
- Protein existence verification
- Metadata retrieval
- Family information extraction
- Non-existent protein handling
- Storage system integration

#### `test_embedding_generator.py`
**Purpose**: Test protein embedding generation

**Test Cases**:
- ESM-2 model integration
- Multiple pooling methods
- Batch processing
- GPU acceleration
- Memory efficiency
- Error handling for invalid sequences

#### `test_network_builder.py`
**Purpose**: Test network construction and visualization

**Test Cases**:
- Network construction methods
- Visualization generation
- Network analysis
- Export functionality
- Interactive plot creation

#### `test_similarity_index.py`
**Purpose**: Test similarity search functionality

**Test Cases**:
- FAISS index creation
- Similarity search accuracy
- Hierarchical indexing
- Streaming capabilities
- Performance optimization

#### `test_storage.py`
**Purpose**: Test data storage functionality

**Test Cases**:
- Hierarchical storage
- Compressed metadata
- Memory-efficient loading
- Streaming capabilities
- Data integrity

#### `test_workflow_orchestrator.py`
**Purpose**: Test complete workflow orchestration

**Test Cases**:
- End-to-end workflow execution
- Performance monitoring
- Memory optimization
- Configuration management
- Error recovery

### Integration Tests

#### `kbase_protein_query_module_query_server_test.py`
**Purpose**: Test complete KBase app functionality

**Test Cases**:
- KBase service integration
- Parameter validation
- Report generation
- Workspace interaction
- Authentication handling
- Error reporting

**Key Test Functions**:
```python
def test_check_protein_existence(self):
    """Test protein existence checking app"""
    params = {'protein_id': 'P12345', 'workspace_name': self.wsName}
    ret = self.serviceImpl.check_protein_existence(self.ctx, params)
    self.assertIsInstance(ret, dict)
    self.assertIn('exists', ret)

def test_generate_protein_embedding(self):
    """Test protein embedding generation app"""
    params = {'sequence': 'MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ'}
    ret = self.serviceImpl.generate_protein_embedding(self.ctx, params)
    self.assertIsInstance(ret, dict)
    self.assertIn('embedding', ret)
```

## Test Data

### Test Data Structure
```
test/data/
├── streaming/                    # Streaming test data
│   ├── sample_embeddings.h5     # Sample protein embeddings
│   ├── sample_metadata.csv      # Sample protein metadata
│   └── test_families/           # Test family data
└── fixtures/                    # Test fixtures
    ├── test_centroids.npz       # Test family centroids
    ├── test_sequences.fasta      # Test protein sequences
    └── test_networks/           # Test network data
```

### Test Data Generation
```bash
# Generate test data
python scripts/generate_dummy_data.py

# Generate specific test fixtures
python -c "
from lib.kbase_protein_query_module.src.storage import ProteinStorage
storage = ProteinStorage()
storage.create_test_data()
"
```

## Test Configuration

### Environment Setup
```bash
# Set up test environment
export KB_DEPLOYMENT_CONFIG=test_local/test.cfg
export KB_AUTH_TOKEN=your_kbase_token
export PYTHONPATH=lib:$PYTHONPATH
```

### Test Configuration File
```ini
# test_local/test.cfg
[kbase_protein_query_module]
workspace-url = https://appdev.kbase.us/services/ws
auth-service-url = https://appdev.kbase.us/services/auth/api/legacy/KBase/Sessions/Login
scratch = /tmp
```

## Performance Testing

### Memory Usage Testing
```python
def test_memory_usage():
    """Test memory usage with large datasets"""
    import psutil
    import gc
    
    process = psutil.Process()
    initial_memory = process.memory_info().rss
    
    # Run memory-intensive operation
    workflow = ProteinNetworkWorkflow()
    results = workflow.run_optimized_workflow(large_sequence)
    
    final_memory = process.memory_info().rss
    memory_increase = final_memory - initial_memory
    
    # Assert memory usage is reasonable
    assert memory_increase < 2 * 1024 * 1024 * 1024  # 2GB limit
```

### Execution Time Testing
```python
def test_execution_time():
    """Test execution time for key operations"""
    import time
    
    start_time = time.time()
    
    # Run operation
    generator = ProteinEmbeddingGenerator()
    embedding = generator.generate_embedding(test_sequence)
    
    execution_time = time.time() - start_time
    
    # Assert reasonable execution time
    assert execution_time < 30.0  # 30 seconds limit
```

## Continuous Integration

### GitHub Actions
```yaml
# .github/workflows/test.yml
name: Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Set up Python
        uses: actions/setup-python@v2
        with:
          python-version: 3.8
      - name: Install dependencies
        run: pip install -r requirements.txt
      - name: Run tests
        run: python -m pytest test/ --cov=kbase_protein_query_module
```

### Local CI
```bash
# Run full test suite locally
./scripts/run_ci_tests.sh
```

## Troubleshooting

### Common Test Issues

1. **Import Errors**
   ```bash
   # Solution: Set PYTHONPATH
   export PYTHONPATH=lib:$PYTHONPATH
   ```

2. **Authentication Errors**
   ```bash
   # Solution: Check token configuration
   cat test_local/test.cfg
   ```

3. **Memory Issues**
   ```bash
   # Solution: Use smaller test datasets
   export TEST_SMALL_DATASET=true
   ```

4. **GPU Issues**
   ```bash
   # Solution: Force CPU usage
   export CUDA_VISIBLE_DEVICES=""
   ```

### Debug Mode
```bash
# Run tests with debug output
python -m pytest test/ -v -s --log-cli-level=DEBUG
```

### Test Isolation
```bash
# Run tests in isolation
python -m pytest test/ --dist=no
```

## Contributing to Tests

### Adding New Tests
1. **Follow Naming Convention**: `test_<module_name>.py`
2. **Use Descriptive Names**: Clear test function names
3. **Include Documentation**: Docstrings for all test functions
4. **Test Edge Cases**: Include boundary condition tests
5. **Mock Dependencies**: Use mocks for external dependencies

### Test Guidelines
- **Isolation**: Each test should be independent
- **Speed**: Tests should run quickly
- **Reliability**: Tests should be deterministic
- **Coverage**: Aim for >90% code coverage
- **Documentation**: Clear test descriptions

### Example Test Template
```python
def test_function_name():
    """
    Test description.
    
    This test verifies that the function behaves correctly
    under normal conditions and edge cases.
    """
    # Arrange
    input_data = create_test_data()
    expected_output = create_expected_output()
    
    # Act
    actual_output = function_under_test(input_data)
    
    # Assert
    assert actual_output == expected_output
    assert len(actual_output) > 0
    assert all(isinstance(x, expected_type) for x in actual_output)
```

## Resources

- **[KBase SDK Testing Guide](https://kbase.github.io/kb_sdk_docs/references/troubleshooting.html)**
- **[pytest Documentation](https://docs.pytest.org/)**
- **[Coverage.py Documentation](https://coverage.readthedocs.io/)**
- **[KBase Developer Guidelines](https://kbase.github.io/kb_sdk_docs/references/developer_guidelines.html)**

## License

This project is licensed under the MIT License - see the [LICENSE](../LICENSE) file for details.
 
