#!/usr/bin/env python3
"""
Comprehensive Test Runner for Reorganized KBase Protein Query Module

This script runs all tests for the reorganized module structure, ensuring
that all components work correctly with the new organization.
"""

import os
import sys
import time
import logging
import subprocess
import importlib
from pathlib import Path
from typing import List, Dict, Any
from datetime import datetime

# Add the lib directory to the Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))

def setup_logging():
    """Setup logging for test execution."""
    # Create output directory first
    os.makedirs('test_outputs', exist_ok=True)
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('test_outputs/reorganized_test_run.log'),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

def run_pytest_tests(test_path: str, test_name: str = None) -> Dict[str, Any]:
    """Run pytest tests for a specific module."""
    logger = logging.getLogger(__name__)
    
    cmd = [
        'python3', '-m', 'pytest', 
        test_path,
        '-v',
        '--tb=short',
        '--capture=no'
    ]
    
    if test_name:
        cmd.extend(['-k', test_name])
    
    logger.info(f"Running tests for {test_path}")
    logger.info(f"Command: {' '.join(cmd)}")
    
    start_time = time.time()
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=os.path.dirname(__file__)
        )
        
        execution_time = time.time() - start_time
        
        return {
            'success': result.returncode == 0,
            'stdout': result.stdout,
            'stderr': result.stderr,
            'returncode': result.returncode,
            'execution_time': execution_time
        }
        
    except Exception as e:
        execution_time = time.time() - start_time
        logger.error(f"Failed to run tests for {test_path}: {str(e)}")
        
        return {
            'success': False,
            'stdout': '',
            'stderr': str(e),
            'returncode': -1,
            'execution_time': execution_time
        }

def test_module_imports() -> Dict[str, Any]:
    """Test that all modules can be imported correctly."""
    logger = logging.getLogger(__name__)
    logger.info("Testing module imports...")
    
    modules_to_test = [
        'kbase_protein_query_module.src.core',
        'kbase_protein_query_module.src.embeddings',
        'kbase_protein_query_module.src.storage',
        'kbase_protein_query_module.src.similarity',
        'kbase_protein_query_module.src.networks',
        'kbase_protein_query_module.src.analysis',
        'kbase_protein_query_module.src.reports',
        'kbase_protein_query_module.src.utils',
        'kbase_protein_query_module.src.workflows',
        'kbase_protein_query_module.src.storage',
        'kbase_protein_query_module.src.stages'
    ]
    
    results = {}
    failed_imports = []
    
    for module_name in modules_to_test:
        try:
            module = importlib.import_module(module_name)
            logger.info(f"✓ Successfully imported {module_name}")
            results[module_name] = {'success': True, 'error': None}
        except Exception as e:
            logger.error(f"✗ Failed to import {module_name}: {str(e)}")
            results[module_name] = {'success': False, 'error': str(e)}
            failed_imports.append(module_name)
    
    return {
        'success': len(failed_imports) == 0,
        'results': results,
        'failed_imports': failed_imports
    }

def test_workflow_orchestrator() -> Dict[str, Any]:
    """Test the new workflow orchestrator."""
    logger = logging.getLogger(__name__)
    logger.info("Testing workflow orchestrator...")
    
    try:
        from kbase_protein_query_module.src.workflows import ProteinQueryWorkflow, WorkflowResult
        from kbase_protein_query_module.src.core import PipelineConfig
        
        # Test basic initialization
        config = PipelineConfig(
            input_proteins=['P12345'],
            perform_embedding_generation=True,
            perform_family_assignment=True
        )
        
        workflow = ProteinQueryWorkflow(config)
        
        logger.info("✓ Successfully created workflow orchestrator")
        
        return {
            'success': True,
            'error': None
        }
        
    except Exception as e:
        logger.error(f"✗ Failed to test workflow orchestrator: {str(e)}")
        return {
            'success': False,
            'error': str(e)
        }

def test_stage_registry() -> Dict[str, Any]:
    """Test the stage registry and dependencies."""
    logger = logging.getLogger(__name__)
    logger.info("Testing stage registry...")
    
    try:
        from kbase_protein_query_module.src.stages import (
            STAGE_REGISTRY, STAGE_DEPENDENCIES,
            get_stage_class, get_stage_dependencies, get_all_stages
        )
        
        # Test stage registry
        all_stages = get_all_stages()
        logger.info(f"✓ Found {len(all_stages)} stages in registry")
        
        # Test stage class retrieval
        for stage_name in all_stages:
            stage_class = get_stage_class(stage_name)
            dependencies = get_stage_dependencies(stage_name)
            logger.info(f"✓ Stage {stage_name}: {len(dependencies)} dependencies")
        
        return {
            'success': True,
            'stages_count': len(all_stages),
            'error': None
        }
        
    except Exception as e:
        logger.error(f"✗ Failed to test stage registry: {str(e)}")
        return {
            'success': False,
            'error': str(e)
        }

def run_unit_tests() -> Dict[str, Any]:
    """Run all unit tests for the reorganized modules."""
    logger = logging.getLogger(__name__)
    logger.info("Running unit tests...")
    
    test_modules = [
        ('test/unit_tests/core', 'Core Module Tests'),
        ('test/unit_tests/embeddings', 'Embeddings Module Tests'),
        ('test/unit_tests/storage', 'Storage Module Tests'),
        ('test/unit_tests/similarity', 'Similarity Module Tests'),
        ('test/unit_tests/networks', 'Networks Module Tests'),
        ('test/unit_tests/analysis', 'Analysis Module Tests'),
        ('test/unit_tests/reports', 'Reports Module Tests'),
        ('test/unit_tests/utils', 'Utils Module Tests'),
        ('test/unit_tests/workflows', 'Workflows Module Tests'),
        ('test/unit_tests/storage', 'Storage Module Tests'),
        ('test/unit_tests/stages', 'Stages Module Tests')
    ]
    
    results = {}
    total_success = 0
    total_failed = 0
    
    for test_path, test_name in test_modules:
        if os.path.exists(test_path):
            result = run_pytest_tests(test_path)
            results[test_name] = result
            
            if result['success']:
                total_success += 1
                logger.info(f"✓ {test_name} passed")
            else:
                total_failed += 1
                logger.error(f"✗ {test_name} failed")
        else:
            logger.warning(f"Test path {test_path} does not exist")
    
    return {
        'success': total_failed == 0,
        'results': results,
        'total_success': total_success,
        'total_failed': total_failed
    }

def generate_test_report(results: Dict[str, Any]) -> str:
    """Generate a comprehensive test report."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    report = f"""
# KBase Protein Query Module - Reorganized Test Report

**Generated:** {timestamp}

## Test Summary

### Module Imports
- **Status:** {'✓ PASSED' if results['imports']['success'] else '✗ FAILED'}
- **Failed Imports:** {len(results['imports']['failed_imports'])}
- **Details:** {', '.join(results['imports']['failed_imports']) if results['imports']['failed_imports'] else 'All imports successful'}

### Workflow Orchestrator
- **Status:** {'✓ PASSED' if results['workflow']['success'] else '✗ FAILED'}
- **Error:** {results['workflow']['error'] or 'None'}

### Stage Registry
- **Status:** {'✓ PASSED' if results['registry']['success'] else '✗ FAILED'}
- **Stages Found:** {results['registry'].get('stages_count', 0)}
- **Error:** {results['registry'].get('error', 'None')}

### Unit Tests
- **Status:** {'✓ PASSED' if results['unit_tests']['success'] else '✗ FAILED'}
- **Passed:** {results['unit_tests']['total_success']}
- **Failed:** {results['unit_tests']['total_failed']}

## Detailed Results

### Unit Test Results
"""
    
    for test_name, result in results['unit_tests']['results'].items():
        status = '✓ PASSED' if result['success'] else '✗ FAILED'
        report += f"- **{test_name}:** {status} ({result['execution_time']:.2f}s)\n"
    
    report += f"""
## Overall Status

**Overall Result:** {'✓ ALL TESTS PASSED' if all([
    results['imports']['success'],
    results['workflow']['success'],
    results['registry']['success'],
    results['unit_tests']['success']
]) else '✗ SOME TESTS FAILED'}

## Recommendations

"""
    
    if not results['imports']['success']:
        report += "- Fix module import issues\n"
    
    if not results['workflow']['success']:
        report += "- Fix workflow orchestrator issues\n"
    
    if not results['registry']['success']:
        report += "- Fix stage registry issues\n"
    
    if not results['unit_tests']['success']:
        report += "- Fix failing unit tests\n"
    
    if all([
        results['imports']['success'],
        results['workflow']['success'],
        results['registry']['success'],
        results['unit_tests']['success']
    ]):
        report += "- All tests passed! The reorganization is complete and working correctly.\n"
    
    return report

def main():
    """Main test execution function."""
    logger = setup_logging()
    
    # Create output directory
    os.makedirs('test_outputs', exist_ok=True)
    
    logger.info("Starting comprehensive test run for reorganized KBase Protein Query Module")
    logger.info("=" * 80)
    
    start_time = time.time()
    
    # Run all tests
    results = {
        'imports': test_module_imports(),
        'workflow': test_workflow_orchestrator(),
        'registry': test_stage_registry(),
        'unit_tests': run_unit_tests()
    }
    
    total_time = time.time() - start_time
    
    # Generate report
    report = generate_test_report(results)
    
    # Save report
    report_file = f'test_outputs/reorganized_test_report_{datetime.now().strftime("%Y%m%d_%H%M%S")}.md'
    with open(report_file, 'w') as f:
        f.write(report)
    
    # Print summary
    logger.info("=" * 80)
    logger.info("TEST EXECUTION COMPLETE")
    logger.info(f"Total execution time: {total_time:.2f} seconds")
    
    overall_success = all([
        results['imports']['success'],
        results['workflow']['success'],
        results['registry']['success'],
        results['unit_tests']['success']
    ])
    
    if overall_success:
        logger.info("✓ ALL TESTS PASSED - Reorganization is complete and working!")
    else:
        logger.error("✗ SOME TESTS FAILED - Please review the report for details")
    
    logger.info(f"Detailed report saved to: {report_file}")
    
    return 0 if overall_success else 1

if __name__ == '__main__':
    sys.exit(main())
