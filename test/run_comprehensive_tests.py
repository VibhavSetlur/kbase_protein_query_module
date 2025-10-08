#!/usr/bin/env python3
"""
Comprehensive test runner for KBase Protein Query Module.

This script runs all tests including unit tests, integration tests,
functional tests, and KBase compliance tests with detailed reporting.
"""

import os
import sys
import subprocess
import argparse
import json
import time
from pathlib import Path
from typing import Dict, List, Any
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class ComprehensiveTestRunner:
    """Comprehensive test runner for the KBase Protein Query Module."""
    
    def __init__(self, project_root: str = None):
        """
        Initialize the comprehensive test runner.
        
        Args:
            project_root: Root directory of the project
        """
        self.project_root = project_root or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.test_dir = os.path.join(self.project_root, 'test')
        self.results = {
            'start_time': time.time(),
            'test_suites': {},
            'overall_success': True,
            'summary': {
                'total_tests': 0,
                'passed': 0,
                'failed': 0,
                'skipped': 0,
                'errors': 0
            }
        }
        
        logger.info(f"Initialized test runner for project: {self.project_root}")
    
    def run_all_tests(self, verbose: bool = False, failfast: bool = False, 
                     categories: List[str] = None) -> Dict[str, Any]:
        """
        Run all test categories.
        
        Args:
            verbose: Enable verbose output
            failfast: Stop on first failure
            categories: Specific categories to run (None for all)
            
        Returns:
            Combined test results
        """
        logger.info("Starting comprehensive test execution...")
        
        test_categories = categories or [
            'unit_tests',
            'integration_tests', 
            'functional_tests',
            'kbase_compliance'
        ]
        
        for category in test_categories:
            logger.info(f"Running {category}...")
            category_result = self._run_test_category(category, verbose, failfast)
            self.results['test_suites'][category] = category_result
            
            # Update overall summary
            self._update_summary(category_result)
            
            # Check for failures
            if not category_result.get('success', False):
                self.results['overall_success'] = False
                if failfast:
                    logger.error(f"Test suite {category} failed. Stopping due to --failfast.")
                    break
        
        self.results['end_time'] = time.time()
        self.results['total_execution_time'] = self.results['end_time'] - self.results['start_time']
        
        logger.info(f"Test execution completed in {self.results['total_execution_time']:.2f} seconds")
        logger.info(f"Overall success: {self.results['overall_success']}")
        
        return self.results
    
    def _run_test_category(self, category: str, verbose: bool, failfast: bool) -> Dict[str, Any]:
        """Run a specific test category."""
        category_dir = os.path.join(self.test_dir, category)
        
        if not os.path.exists(category_dir):
            logger.warning(f"Test category directory not found: {category_dir}")
            return {
                'success': False,
                'error': f'Category directory not found: {category_dir}',
                'tests_run': 0,
                'tests_passed': 0,
                'tests_failed': 0,
                'tests_skipped': 0,
                'tests_errors': 0
            }
        
        # Use pytest for running tests
        cmd = ['python', '-m', 'pytest', category_dir]
        
        if verbose:
            cmd.append('-v')
        
        if failfast:
            cmd.append('-x')
        
        # Add coverage if available
        cmd.extend(['--tb=short', '--strict-markers'])
        
        logger.info(f"Running command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=self.project_root,
                timeout=300  # 5 minute timeout per category
            )
            
            # Parse pytest output
            return self._parse_pytest_output(result, category)
            
        except subprocess.TimeoutExpired:
            logger.error(f"Test category {category} timed out")
            return {
                'success': False,
                'error': 'Test execution timed out',
                'tests_run': 0,
                'tests_passed': 0,
                'tests_failed': 0,
                'tests_skipped': 0,
                'tests_errors': 0
            }
        except Exception as e:
            logger.error(f"Error running test category {category}: {e}")
            return {
                'success': False,
                'error': str(e),
                'tests_run': 0,
                'tests_passed': 0,
                'tests_failed': 0,
                'tests_skipped': 0,
                'tests_errors': 0
            }
    
    def _parse_pytest_output(self, result: subprocess.CompletedProcess, category: str) -> Dict[str, Any]:
        """Parse pytest output to extract test results."""
        success = result.returncode == 0
        
        # Basic parsing of pytest output
        output_lines = result.stdout.split('\n')
        error_lines = result.stderr.split('\n')
        
        # Extract test counts from output
        tests_run = 0
        tests_passed = 0
        tests_failed = 0
        tests_skipped = 0
        tests_errors = 0
        
        for line in output_lines:
            if 'failed' in line and 'passed' in line:
                # Parse summary line like "5 failed, 10 passed in 2.34s"
                parts = line.split(',')
                for part in parts:
                    part = part.strip()
                    if 'failed' in part:
                        tests_failed = int(part.split()[0])
                    elif 'passed' in part:
                        tests_passed = int(part.split()[0])
                    elif 'skipped' in part:
                        tests_skipped = int(part.split()[0])
                    elif 'error' in part:
                        tests_errors = int(part.split()[0])
        
        tests_run = tests_passed + tests_failed + tests_errors
        
        return {
            'success': success,
            'tests_run': tests_run,
            'tests_passed': tests_passed,
            'tests_failed': tests_failed,
            'tests_skipped': tests_skipped,
            'tests_errors': tests_errors,
            'stdout': result.stdout,
            'stderr': result.stderr,
            'return_code': result.returncode
        }
    
    def _update_summary(self, category_result: Dict[str, Any]):
        """Update overall summary with category results."""
        self.results['summary']['total_tests'] += category_result.get('tests_run', 0)
        self.results['summary']['passed'] += category_result.get('tests_passed', 0)
        self.results['summary']['failed'] += category_result.get('tests_failed', 0)
        self.results['summary']['skipped'] += category_result.get('tests_skipped', 0)
        self.results['summary']['errors'] += category_result.get('tests_errors', 0)
    
    def generate_comprehensive_report(self, output_file: str = None) -> str:
        """
        Generate a comprehensive test report.
        
        Args:
            output_file: Output file path
            
        Returns:
            Path to generated report file
        """
        if output_file is None:
            output_file = os.path.join(self.test_dir, 'comprehensive_test_report.json')
        
        # Generate JSON report
        with open(output_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        
        # Generate human-readable report
        readable_report = output_file.replace('.json', '.txt')
        with open(readable_report, 'w') as f:
            f.write("KBase Protein Query Module - Comprehensive Test Report\n")
            f.write("=" * 60 + "\n\n")
            
            f.write(f"Test Execution Time: {self.results['total_execution_time']:.2f} seconds\n")
            f.write(f"Overall Success: {self.results['overall_success']}\n\n")
            
            f.write("Summary:\n")
            f.write(f"  Total Tests: {self.results['summary']['total_tests']}\n")
            f.write(f"  Passed: {self.results['summary']['passed']}\n")
            f.write(f"  Failed: {self.results['summary']['failed']}\n")
            f.write(f"  Skipped: {self.results['summary']['skipped']}\n")
            f.write(f"  Errors: {self.results['summary']['errors']}\n\n")
            
            f.write("Test Suite Results:\n")
            f.write("-" * 40 + "\n")
            for category, result in self.results['test_suites'].items():
                f.write(f"{category}:\n")
                f.write(f"  Success: {result['success']}\n")
                f.write(f"  Tests Run: {result['tests_run']}\n")
                f.write(f"  Passed: {result['tests_passed']}\n")
                f.write(f"  Failed: {result['tests_failed']}\n")
                f.write(f"  Skipped: {result['tests_skipped']}\n")
                f.write(f"  Errors: {result['tests_errors']}\n")
                
                if not result['success'] and 'error' in result:
                    f.write(f"  Error: {result['error']}\n")
                
                f.write("\n")
        
        logger.info(f"Comprehensive test report generated: {readable_report}")
        logger.info(f"JSON report generated: {output_file}")
        
        return readable_report
    
    def run_kbase_server_tests(self, verbose: bool = False) -> Dict[str, Any]:
        """
        Run the original KBase server tests.
        
        Args:
            verbose: Enable verbose output
            
        Returns:
            Server test results
        """
        logger.info("Running KBase server tests...")
        
        server_test_file = os.path.join(self.test_dir, 'kbase_protein_query_module_server_test.py')
        
        if not os.path.exists(server_test_file):
            logger.warning("KBase server test file not found")
            return {
                'success': False,
                'error': 'Server test file not found',
                'tests_run': 0,
                'tests_passed': 0,
                'tests_failed': 0
            }
        
        # Use unittest for server tests
        cmd = ['python', '-m', 'unittest', 'discover', '-s', self.test_dir, '-p', '*server_test.py']
        
        if verbose:
            cmd.append('-v')
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=self.project_root,
                timeout=180  # 3 minute timeout
            )
            
            success = result.returncode == 0
            
            return {
                'success': success,
                'stdout': result.stdout,
                'stderr': result.stderr,
                'return_code': result.returncode
            }
            
        except Exception as e:
            logger.error(f"Error running server tests: {e}")
            return {
                'success': False,
                'error': str(e)
            }
    
    def run_linting_tests(self) -> Dict[str, Any]:
        """Run code quality and linting tests."""
        logger.info("Running linting tests...")
        
        # Check for common linting tools
        linting_results = {}
        
        # Python syntax check
        try:
            result = subprocess.run(
                ['python', '-m', 'py_compile'] + 
                [os.path.join(self.project_root, 'lib', 'kbase_protein_query_module', 'kbase_protein_query_moduleImpl.py'),
                 os.path.join(self.project_root, 'lib', 'kbase_protein_query_module', 'kbase_protein_query_moduleServer.py')],
                capture_output=True,
                text=True,
                cwd=self.project_root
            )
            
            linting_results['syntax_check'] = {
                'success': result.returncode == 0,
                'stdout': result.stdout,
                'stderr': result.stderr
            }
            
        except Exception as e:
            linting_results['syntax_check'] = {
                'success': False,
                'error': str(e)
            }
        
        return linting_results


def main():
    """Main entry point for comprehensive test runner."""
    parser = argparse.ArgumentParser(description='Comprehensive Test Runner for KBase Protein Query Module')
    parser.add_argument('--project-root', default=None, help='Project root directory')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    parser.add_argument('--failfast', '-x', action='store_true', help='Stop on first failure')
    parser.add_argument('--categories', nargs='+', 
                       choices=['unit_tests', 'integration_tests', 'functional_tests', 'kbase_compliance'],
                       help='Specific test categories to run')
    parser.add_argument('--server-tests', action='store_true', help='Run KBase server tests')
    parser.add_argument('--linting', action='store_true', help='Run linting tests')
    parser.add_argument('--output-file', default=None, help='Output file for test report')
    
    args = parser.parse_args()
    
    # Initialize test runner
    runner = ComprehensiveTestRunner(args.project_root)
    
    # Run tests based on arguments
    if args.server_tests:
        # Run only server tests
        result = runner.run_kbase_server_tests(args.verbose)
        print("Server Test Results:")
        print(f"Success: {result['success']}")
        if not result['success']:
            print(f"Error: {result.get('error', 'Unknown error')}")
            print(f"Output: {result.get('stdout', '')}")
            print(f"Errors: {result.get('stderr', '')}")
        return 0 if result['success'] else 1
    
    if args.linting:
        # Run only linting tests
        results = runner.run_linting_tests()
        print("Linting Results:")
        for check, result in results.items():
            print(f"{check}: {'PASS' if result['success'] else 'FAIL'}")
            if not result['success']:
                print(f"  Error: {result.get('error', result.get('stderr', 'Unknown error'))}")
        return 0 if all(r['success'] for r in results.values()) else 1
    
    # Run comprehensive tests
    results = runner.run_all_tests(args.verbose, args.failfast, args.categories)
    
    # Generate report
    runner.generate_comprehensive_report(args.output_file)
    
    # Print summary
    print("\n" + "="*60)
    print("COMPREHENSIVE TEST SUMMARY")
    print("="*60)
    print(f"Overall Success: {results['overall_success']}")
    print(f"Total Tests: {results['summary']['total_tests']}")
    print(f"Passed: {results['summary']['passed']}")
    print(f"Failed: {results['summary']['failed']}")
    print(f"Skipped: {results['summary']['skipped']}")
    print(f"Errors: {results['summary']['errors']}")
    print(f"Execution Time: {results['total_execution_time']:.2f} seconds")
    
    print("\nTest Suite Results:")
    for category, result in results['test_suites'].items():
        status = "PASS" if result['success'] else "FAIL"
        print(f"  {category}: {status} ({result['tests_passed']}/{result['tests_run']} passed)")
    
    # Return appropriate exit code
    return 0 if results['overall_success'] else 1


if __name__ == '__main__':
    sys.exit(main())

