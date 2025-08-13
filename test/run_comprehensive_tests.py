#!/usr/bin/env python3
"""
Comprehensive Test Runner for KBase Protein Query Module

This script runs all tests for the reorganized codebase, including:
- Unit tests for all modules
- Integration tests
- Performance tests
- End-to-end workflow tests

Updated to reflect the new professional src structure:
- core: Core abstractions and configurations
- data: Data models and reference data
- processing: All processing logic (embeddings, similarity, networks)
- analysis: Analysis components and algorithms
- storage: Data storage and persistence layer
- reports: Report generation and visualization
- stages: Pipeline stages (input, processing, output)
- utils: Utility functions and helpers
- workflows: Workflow orchestration and management
"""

import os
import sys
import unittest
import time
import tempfile
import shutil
from pathlib import Path
from datetime import datetime
import json

# Add the lib directory to the Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))

def run_comprehensive_tests():
    """Run all comprehensive tests for the reorganized codebase."""
    
    print("=" * 80)
    print("COMPREHENSIVE TEST SUITE - KBase Protein Query Module")
    print("=" * 80)
    print(f"Test run started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Create test output directory
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    test_output_dir = f"test_outputs/comprehensive_test_run_{timestamp}"
    os.makedirs(test_output_dir, exist_ok=True)
    
    # Test results storage
    test_results = {
        'timestamp': timestamp,
        'modules': {},
        'summary': {
            'total_tests': 0,
            'passed': 0,
            'failed': 0,
            'errors': 0,
            'skipped': 0
        },
        'performance': {},
        'errors': []
    }
    
    # Define test modules with their descriptions
    test_modules = [
        ('test/unit_tests/core', 'Core Module Tests'),
        ('test/unit_tests/data', 'Data Module Tests'),
        ('test/unit_tests/processing', 'Processing Module Tests'),
        ('test/unit_tests/analysis', 'Analysis Module Tests'),
        ('test/unit_tests/storage', 'Storage Module Tests'),
        ('test/unit_tests/reports', 'Reports Module Tests'),
        ('test/unit_tests/stages', 'Pipeline Stages Tests'),
        ('test/unit_tests/utils', 'Utils Module Tests'),
        ('test/unit_tests/workflows', 'Workflow Tests'),
        ('test/integration_tests', 'Integration Tests')
    ]
    
    # Run tests for each module
    for test_dir, module_name in test_modules:
        print(f"\n{'='*60}")
        print(f"Testing: {module_name}")
        print(f"Directory: {test_dir}")
        print(f"{'='*60}")
        
        if not os.path.exists(test_dir):
            print(f"WARNING: Test directory {test_dir} does not exist. Skipping.")
            test_results['modules'][module_name] = {
                'status': 'skipped',
                'reason': 'Directory not found',
                'tests_run': 0,
                'failures': 0,
                'errors': 0,
                'skipped': 0
            }
            continue
        
        # Create test suite for this module
        loader = unittest.TestLoader()
        suite = loader.discover(test_dir, pattern='test_*.py')
        
        # Run tests
        runner = unittest.TextTestRunner(verbosity=2, stream=sys.stdout)
        start_time = time.time()
        
        try:
            result = runner.run(suite)
            end_time = time.time()
            
            # Store results
            module_results = {
                'status': 'completed',
                'tests_run': result.testsRun,
                'failures': len(result.failures),
                'errors': len(result.errors),
                'skipped': 0,
                'execution_time': end_time - start_time
            }
            
            # Update summary
            test_results['summary']['total_tests'] += result.testsRun
            test_results['summary']['failed'] += len(result.failures)
            test_results['summary']['errors'] += len(result.errors)
            
            # Store detailed results
            if result.failures:
                module_results['failure_details'] = [
                    {'test': str(failure[0]), 'traceback': failure[1]} 
                    for failure in result.failures
                ]
            
            if result.errors:
                module_results['error_details'] = [
                    {'test': str(error[0]), 'traceback': error[1]} 
                    for error in result.errors
                ]
            
            test_results['modules'][module_name] = module_results
            test_results['performance'][module_name] = end_time - start_time
            
            # Print module summary
            print(f"\n{module_name} Results:")
            print(f"  Tests run: {result.testsRun}")
            print(f"  Failures: {len(result.failures)}")
            print(f"  Errors: {len(result.errors)}")
            print(f"  Execution time: {end_time - start_time:.2f} seconds")
            
        except Exception as e:
            print(f"ERROR: Failed to run tests for {module_name}: {e}")
            test_results['modules'][module_name] = {
                'status': 'error',
                'reason': str(e),
                'tests_run': 0,
                'failures': 0,
                'errors': 1,
                'skipped': 0
            }
            test_results['summary']['errors'] += 1
            test_results['errors'].append({
                'module': module_name,
                'error': str(e)
            })
    
    # Calculate passed tests
    test_results['summary']['passed'] = (
        test_results['summary']['total_tests'] - 
        test_results['summary']['failed'] - 
        test_results['summary']['errors'] - 
        test_results['summary']['skipped']
    )
    
    # Print comprehensive summary
    print("\n" + "="*80)
    print("COMPREHENSIVE TEST SUMMARY")
    print("="*80)
    print(f"Total tests run: {test_results['summary']['total_tests']}")
    print(f"Passed: {test_results['summary']['passed']}")
    print(f"Failed: {test_results['summary']['failed']}")
    print(f"Errors: {test_results['summary']['errors']}")
    print(f"Skipped: {test_results['summary']['skipped']}")
    
    # Calculate success rate
    if test_results['summary']['total_tests'] > 0:
        success_rate = (test_results['summary']['passed'] / test_results['summary']['total_tests']) * 100
        print(f"Success rate: {success_rate:.1f}%")
    
    # Print performance summary
    print("\nPerformance Summary:")
    total_time = sum(test_results['performance'].values())
    print(f"Total execution time: {total_time:.2f} seconds")
    
    for module, exec_time in test_results['performance'].items():
        print(f"  {module}: {exec_time:.2f} seconds")
    
    # Print detailed results for failed modules
    failed_modules = [
        (name, results) for name, results in test_results['modules'].items()
        if results.get('failures', 0) > 0 or results.get('errors', 0) > 0
    ]
    
    if failed_modules:
        print("\nFailed Modules:")
        for module_name, module_results in failed_modules:
            print(f"  {module_name}:")
            print(f"    Failures: {module_results.get('failures', 0)}")
            print(f"    Errors: {module_results.get('errors', 0)}")
            
            if 'failure_details' in module_results:
                print("    Failure details:")
                for failure in module_results['failure_details'][:3]:  # Show first 3
                    print(f"      - {failure['test']}")
            
            if 'error_details' in module_results:
                print("    Error details:")
                for error in module_results['error_details'][:3]:  # Show first 3
                    print(f"      - {error['test']}")
    
    # Save detailed results
    results_file = os.path.join(test_output_dir, 'comprehensive_test_results.json')
    with open(results_file, 'w') as f:
        json.dump(test_results, f, indent=2)
    
    # Save summary report
    summary_file = os.path.join(test_output_dir, 'test_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("COMPREHENSIVE TEST SUMMARY\n")
        f.write("=" * 50 + "\n")
        f.write(f"Test run: {timestamp}\n")
        f.write(f"Total tests: {test_results['summary']['total_tests']}\n")
        f.write(f"Passed: {test_results['summary']['passed']}\n")
        f.write(f"Failed: {test_results['summary']['failed']}\n")
        f.write(f"Errors: {test_results['summary']['errors']}\n")
        f.write(f"Skipped: {test_results['summary']['skipped']}\n")
        f.write(f"Total time: {total_time:.2f} seconds\n")
        
        if test_results['summary']['total_tests'] > 0:
            success_rate = (test_results['summary']['passed'] / test_results['summary']['total_tests']) * 100
            f.write(f"Success rate: {success_rate:.1f}%\n")
    
    print(f"\nDetailed results saved to: {test_output_dir}")
    print(f"Results file: {results_file}")
    print(f"Summary file: {summary_file}")
    
    # Return overall success
    overall_success = (
        test_results['summary']['failed'] == 0 and 
        test_results['summary']['errors'] == 0
    )
    
    if overall_success:
        print("\n🎉 ALL TESTS PASSED! 🎉")
        return 0
    else:
        print("\n❌ SOME TESTS FAILED ❌")
        return 1

def run_specific_module_tests(module_name):
    """Run tests for a specific module."""
    module_test_dirs = {
        'core': 'test/unit_tests/core',
        'data': 'test/unit_tests/data',
        'processing': 'test/unit_tests/processing',
        'analysis': 'test/unit_tests/analysis',
        'storage': 'test/unit_tests/storage',
        'reports': 'test/unit_tests/reports',
        'stages': 'test/unit_tests/stages',
        'utils': 'test/unit_tests/utils',
        'workflows': 'test/unit_tests/workflows',
        'integration': 'test/integration_tests'
    }
    
    if module_name not in module_test_dirs:
        print(f"Unknown module: {module_name}")
        print(f"Available modules: {', '.join(module_test_dirs.keys())}")
        return 1
    
    test_dir = module_test_dirs[module_name]
    print(f"Running tests for module: {module_name}")
    print(f"Test directory: {test_dir}")
    
    if not os.path.exists(test_dir):
        print(f"Test directory {test_dir} does not exist.")
        return 1
    
    loader = unittest.TestLoader()
    suite = loader.discover(test_dir, pattern='test_*.py')
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return 0 if result.wasSuccessful() else 1

if __name__ == '__main__':
    if len(sys.argv) > 1:
        # Run specific module tests
        module_name = sys.argv[1].lower()
        exit_code = run_specific_module_tests(module_name)
    else:
        # Run all comprehensive tests
        exit_code = run_comprehensive_tests()
    
    sys.exit(exit_code)
