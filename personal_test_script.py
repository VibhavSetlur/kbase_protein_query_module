#!/usr/bin/env python3
"""
Personal Test Script for KBase Protein Query Module

This script tests the protein query analysis workflow outside of the KBase environment
by bypassing authentication and using direct protein sequence/UniProt ID inputs.

Usage:
    python personal_test_script.py
"""

import os
import sys
import json
import time
import logging
import tempfile
from pathlib import Path

# Add the lib directory to the path so we can import the modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'lib'))

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Mock KBase environment variables to avoid authentication
os.environ['PYTEST_CURRENT_TEST'] = '1'
os.environ['KPQM_TEST_FAST'] = '1'

def create_test_data():
    """Create test protein sequences and UniProt IDs for testing."""
    return {
        'protein_sequences': [
            {
                'id': 'test_protein_1',
                'sequence': 'MKLLILTCLVAVALARPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNEL',
                'description': 'Test protein 1 - Ras-related protein'
            },
            {
                'id': 'test_protein_2', 
                'sequence': 'MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT',
                'description': 'Test protein 2 - Insulin-like protein'
            }
        ],
        'uniprot_ids': [
            {
                'id': 'P01116',
                'description': 'HRAS - GTPase HRas'
            },
            {
                'id': 'P01308',
                'description': 'INS - Insulin'
            }
        ]
    }

def mock_kbase_clients():
    """Mock KBase clients to avoid authentication issues."""
    try:
        from unittest.mock import MagicMock, patch
        
        # Mock the KBase clients
        with patch('installed_clients.KBaseReportClient.KBaseReport') as mock_report, \
             patch('installed_clients.WorkspaceClient.Workspace') as mock_workspace, \
             patch('installed_clients.DataFileUtilClient.DataFileUtil') as mock_dfu, \
             patch('installed_clients.KBUtilLibClient.KBUtilLib') as mock_kb_util:
            
            # Configure mock return values
            mock_report.return_value.create_extended_report.return_value = {
                'name': 'test_report',
                'ref': 'test_report_ref'
            }
            
            mock_dfu.return_value.file_to_shock.return_value = {
                'shock_id': 'test_shock_id',
                'node_file_name': 'test_archive.zip',
                'shock_url': 'https://shock.test/node/test_shock_id'
            }
            
            return True
    except ImportError:
        logger.warning("unittest.mock not available, using fallback mocking")
        return False

def run_protein_sequence_test():
    """Test with direct protein sequences."""
    logger.info("=== Testing with Protein Sequences ===")
    
    test_data = create_test_data()
    
    # Prepare input data for protein sequences - convert to the expected format
    protein_input = []
    for seq_data in test_data['protein_sequences']:
        protein_input.append(f">{seq_data['id']}\n{seq_data['sequence']}")
    
    input_data = {
        'input_type': 'protein_input',
        'protein_input': protein_input
    }
    
    return run_analysis(input_data, 'protein_sequence_test')

def run_uniprot_test():
    """Test with UniProt IDs."""
    logger.info("=== Testing with UniProt IDs ===")
    
    test_data = create_test_data()
    
    # Prepare input data for UniProt IDs
    input_data = {
        'input_type': 'uniprot_ids',
        'uniprot_ids': [item['id'] for item in test_data['uniprot_ids']]
    }
    
    return run_analysis(input_data, 'uniprot_test')

def run_analysis(input_data, test_name):
    """Run the analysis with given input data."""
    try:
        # Import the workflow orchestrator
        from kbase_protein_query_module.src.core.workflow_orchestrator import WorkflowOrchestrator
        from kbase_protein_query_module.src.core.pipeline_config import PipelineConfig
        
        # Create output directory
        output_dir = os.path.join('personal_test_output', test_name)
        os.makedirs(output_dir, exist_ok=True)
        
        # Create pipeline configuration with selected analyses
        config = {
            'output_dir': output_dir,
            'workspace_name': 'test_workspace',
            'analysis_name': test_name,
            'run_id': f'test_{int(time.time())}',
            'enabled_input_types': ['protein_input', 'uniprot_ids', 'workspace_object'],
            'enabled_stages': ['input_processing', 'analysis', 'output_generation'],
            'selected_analyses': ['network_analysis']  # Enable network analysis
        }
        
        pipeline_config = PipelineConfig()
        # Update config attributes
        for key, value in config.items():
            if hasattr(pipeline_config, key):
                setattr(pipeline_config, key, value)
        
        # Create workflow orchestrator
        workflow = WorkflowOrchestrator(config=pipeline_config, kb_util=None)
        
        logger.info(f"Starting analysis: {test_name}")
        logger.info(f"Input data: {json.dumps(input_data, indent=2)}")
        
        # Execute workflow
        result = workflow.execute(input_data)
        
        logger.info(f"Analysis completed: {result.success}")
        if result.success:
            logger.info(f"Execution time: {result.execution_time:.2f} seconds")
            logger.info(f"Analyses completed: {result.analyses_completed}")
            logger.info(f"Output directory: {result.output_directory}")
            
            # Check if output files were created
            check_output_files(result.output_directory)
            
            # Try to generate zip file manually for testing
            try_generate_zip_file(result.output_directory, test_name)
            
            return result
        else:
            logger.error(f"Analysis failed: {result.error_message}")
            return None
            
    except Exception as e:
        logger.error(f"Error running analysis: {e}")
        import traceback
        traceback.print_exc()
        return None

def check_output_files(output_dir):
    """Check if expected output files were created."""
    logger.info(f"Checking output files in: {output_dir}")
    
    if not os.path.exists(output_dir):
        logger.error(f"Output directory does not exist: {output_dir}")
        return
    
    # List all files in the output directory
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            file_path = os.path.join(root, file)
            file_size = os.path.getsize(file_path)
            logger.info(f"  {file_path} ({file_size} bytes)")
    
    # Look for zip files specifically
    zip_files = []
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            if file.endswith('.zip'):
                zip_files.append(os.path.join(root, file))
    
    if zip_files:
        logger.info(f"Found {len(zip_files)} zip file(s):")
        for zip_file in zip_files:
            zip_size = os.path.getsize(zip_file)
            logger.info(f"  {zip_file} ({zip_size} bytes)")
    else:
        logger.warning("No zip files found in output")

def try_generate_zip_file(output_dir, test_name):
    """Try to manually generate a zip file for testing purposes."""
    try:
        import zipfile
        
        # Create zip file in the personal test output directory
        zip_path = os.path.join('personal_test_output', f'{test_name}_output.zip')
        
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, dirs, files in os.walk(output_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    # Calculate relative path for zip archive
                    arc_path = os.path.relpath(file_path, output_dir)
                    zipf.write(file_path, arc_path)
                    logger.info(f"Added to zip: {arc_path}")
        
        zip_size = os.path.getsize(zip_path)
        logger.info(f"✓ Generated zip file: {zip_path} ({zip_size} bytes)")
        
    except Exception as e:
        logger.error(f"Failed to generate zip file: {e}")

def main():
    """Main test function."""
    logger.info("Starting Personal Test Script for KBase Protein Query Module")
    logger.info("=" * 60)
    
    # Create personal test output directory
    os.makedirs('personal_test_output', exist_ok=True)
    
    # Mock KBase clients if possible
    mock_kbase_clients()
    
    # Run tests
    results = []
    
    # Test 1: Protein sequences
    result1 = run_protein_sequence_test()
    results.append(('protein_sequences', result1))
    
    # Test 2: UniProt IDs
    result2 = run_uniprot_test()
    results.append(('uniprot_ids', result2))
    
    # Summary
    logger.info("=" * 60)
    logger.info("TEST SUMMARY")
    logger.info("=" * 60)
    
    for test_name, result in results:
        if result and result.success:
            logger.info(f"✓ {test_name}: PASSED")
        else:
            logger.info(f"✗ {test_name}: FAILED")
    
    logger.info(f"\nOutput directory: {os.path.abspath('personal_test_output')}")
    logger.info("Test completed!")

if __name__ == '__main__':
    main()
