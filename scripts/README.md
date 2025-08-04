# KBase Protein Query Module - Scripts

This directory contains utility scripts for the KBase Protein Query Module.

## Scripts Overview

### Core Scripts

#### `create_faiss_indexes.py`
**Purpose**: Create FAISS indexes for protein family data  
**Usage**:
```bash
# Create indexes for existing data
python scripts/create_faiss_indexes.py

# Check data availability only
python scripts/create_faiss_indexes.py --check-only

# Force recreate all indexes
python scripts/create_faiss_indexes.py --force
```

**Features**:
- Creates FAISS IVF indexes for protein similarity search
- Generates both float32 and binary indexes
- Creates metadata files for each index
- Handles large-scale protein family data
- Includes data validation and error handling

### Deployment Scripts

#### `entrypoint.sh`
**Purpose**: Main entry point for Docker container  
**Usage**: Automatically called by Docker container

**Features**:
- Sets up environment variables
- Handles different execution modes (test, async, init, bash)
- Configures deployment settings
- Manages authentication tokens

#### `prepare_deploy_cfg.py`
**Purpose**: Prepare deployment configuration  
**Usage**: Automatically called by entrypoint.sh

**Features**:
- Processes deployment configuration templates
- Sets up KBase service endpoints
- Configures authentication and workspace URLs
- Handles environment-specific settings

#### `run_async.sh`
**Purpose**: Run asynchronous jobs  
**Usage**: Called by entrypoint.sh for async job execution

**Features**:
- Executes async job processing
- Handles input/output JSON files
- Manages authentication tokens
- Provides error handling and logging

## Script Categories

### 1. Data Processing Scripts
- **Index Creation**: Generate FAISS indexes for similarity search
- **Data Validation**: Check data availability and integrity
- **Format Conversion**: Convert between different data formats

### 2. Deployment Scripts
- **Container Management**: Docker entry points and configuration
- **Service Setup**: KBase service registration and configuration
- **Environment Management**: Environment variable handling

### 3. Utility Scripts
- **Configuration**: Setup and validation scripts
- **Testing**: Test data generation and validation
- **Documentation**: Automated documentation generation

## Usage Examples

### Creating FAISS Indexes
```bash
# Basic usage
python scripts/create_faiss_indexes.py

# Check if data is available
python scripts/create_faiss_indexes.py --check-only

# Force recreate all indexes
python scripts/create_faiss_indexes.py --force
```

### Docker Deployment
```bash
# Build and run container
docker build -t kbase_protein_query_module .
docker run -it kbase_protein_query_module

# Run tests in container
docker run -it kbase_protein_query_module test

# Run async job
docker run -it kbase_protein_query_module async
```

## Script Requirements

### Dependencies
- **Python 3.8+**: All Python scripts
- **FAISS**: For index creation (`pip install faiss-cpu`)
- **NumPy**: For numerical operations
- **H5Py**: For HDF5 file handling
- **Pathlib**: For file path operations

### Data Requirements
- **Family Data**: H5 files in `data/families/`
- **Model Files**: ESM-2 model in `data/esm2_t6_8M_UR50D_local/`
- **Configuration**: Proper deployment configuration

## Error Handling

### Common Issues
1. **Missing Dependencies**: Install required Python packages
2. **Data Not Found**: Ensure data files are in correct locations
3. **Permission Issues**: Check file permissions and ownership
4. **Memory Issues**: Ensure sufficient RAM for large datasets

### Debugging
- **Logs**: Check script output for error messages
- **Data Validation**: Use `--check-only` flag to validate data
- **Configuration**: Verify deployment configuration files
- **Environment**: Check environment variables and paths

## Best Practices

### Script Development
- **Error Handling**: Always include proper error handling
- **Logging**: Use consistent logging format
- **Documentation**: Include docstrings and usage examples
- **Validation**: Validate inputs and check data availability

### Deployment
- **Configuration**: Use environment-specific configurations
- **Security**: Handle authentication tokens securely
- **Monitoring**: Include proper logging and monitoring
- **Testing**: Test scripts in containerized environment

## Integration with KBase

### Service Integration
- **Workspace Client**: Proper workspace object handling
- **Authentication**: Secure token management
- **Error Reporting**: KBase-compatible error messages
- **Logging**: KBase logging standards

### Data Flow
1. **Input Validation**: Check input parameters and data
2. **Processing**: Execute core functionality
3. **Output Generation**: Create KBase-compatible outputs
4. **Error Handling**: Provide meaningful error messages

## References

- [KBase SDK Documentation](https://kbase.github.io/kb_sdk_docs/)
- [FAISS Documentation](https://github.com/facebookresearch/faiss)
- [Docker Best Practices](https://docs.docker.com/develop/dev-best-practices/)
- [Python Scripting Guidelines](https://docs.python.org/3/tutorial/)

---

**Author**: Vibhav Setlur  
**Contact**: https://kbase.us/contact-us/  
**Last Updated**: 2024 