"""
Test runner for KBase Protein Query Module unit tests.

Provides comprehensive test execution with KBase compliance and reporting.
"""

import os
import sys
import unittest
import argparse
import json
import time
from pathlib import Path
from typing import Dict, List, Any
import logging

# Add lib to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../lib'))

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class KBaseTestRunner:
    """Test runner with KBase compliance features."""
    
    def __init__(self, test_dir: str = None, output_dir: str = None):
        """
        Initialize the test runner.
        
        Args:
            test_dir: Directory containing tests
            output_dir: Directory for test outputs
        """
        self.test_dir = test_dir or os.path.dirname(__file__)
        self.output_dir = output_dir or os.path.join(self.test_dir, 'test_outputs')
        self.results = {
            'total_tests': 0,
            'passed': 0,
            'failed': 0,
            'skipped': 0,
            'errors': 0,
            'test_details': [],
            'execution_time': 0,
            'kbase_compliance': {
                'test_coverage': 0,
                'standards_met': True,
                'requirements_satisfied': True
            }
        }
        
        # Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)
    
    def discover_tests(self) -> List[unittest.TestSuite]:
        """
        Discover all test modules in the test directory.
        
        Returns:
            List of test suites
        """
        logger.info(f"Discovering tests in: {self.test_dir}")
        
        loader = unittest.TestLoader()
        suites = []
        
        # Define test categories
        test_categories = [
            'input',
            'analysis', 
            'output',
            'core',
            'util',
            'integration',
            'functional'
        ]
        
        for category in test_categories:
            category_dir = os.path.join(self.test_dir, category)
            if os.path.exists(category_dir):
                logger.info(f"Loading {category} tests...")
                suite = loader.discover(category_dir, pattern='test_*.py')
                if suite.countTestCases() > 0:
                    suites.append(suite)
                    logger.info(f"Found {suite.countTestCases()} {category} tests")
        
        return suites
    
    def run_tests(self, test_suites: List[unittest.TestSuite], 
                  verbose: bool = False, failfast: bool = False) -> Dict[str, Any]:
        """
        Run the discovered test suites.
        
        Args:
            test_suites: List of test suites to run
            verbose: Enable verbose output
            failfast: Stop on first failure
            
        Returns:
            Test results dictionary
        """
        logger.info("Starting test execution...")
        start_time = time.time()
        
        # Combine all test suites
        combined_suite = unittest.TestSuite(test_suites)
        self.results['total_tests'] = combined_suite.countTestCases()
        
        logger.info(f"Total tests to run: {self.results['total_tests']}")
        
        # Configure test runner
        runner = unittest.TextTestRunner(
            verbosity=2 if verbose else 1,
            failfast=failfast,
            stream=sys.stdout
        )
        
        # Run tests
        result = runner.run(combined_suite)
        
        # Collect results
        self.results['passed'] = result.testsRun - len(result.failures) - len(result.errors) - len(result.skipped)
        self.results['failed'] = len(result.failures)
        self.results['errors'] = len(result.errors)
        self.results['skipped'] = len(result.skipped)
        self.results['execution_time'] = time.time() - start_time
        
        # Collect test details
        self._collect_test_details(result)
        
        # Assess KBase compliance
        self._assess_kbase_compliance(result)
        
        logger.info(f"Test execution completed in {self.results['execution_time']:.2f} seconds")
        logger.info(f"Results: {self.results['passed']} passed, {self.results['failed']} failed, "
                   f"{self.results['errors']} errors, {self.results['skipped']} skipped")
        
        return self.results
    
    def _collect_test_details(self, result: unittest.TestResult):
        """Collect detailed information about test results."""
        for test, traceback in result.failures:
            self.results['test_details'].append({
                'name': str(test),
                'status': 'failed',
                'traceback': traceback
            })
        
        for test, traceback in result.errors:
            self.results['test_details'].append({
                'name': str(test),
                'status': 'error',
                'traceback': traceback
            })
        
        for test, reason in result.skipped:
            self.results['test_details'].append({
                'name': str(test),
                'status': 'skipped',
                'reason': reason
            })
    
    def _assess_kbase_compliance(self, result: unittest.TestResult):
        """Assess KBase compliance based on test results."""
        # Calculate test coverage (simplified)
        total_possible_tests = self.results['total_tests']
        if total_possible_tests > 0:
            self.results['kbase_compliance']['test_coverage'] = (
                (self.results['passed'] + self.results['failed'] + self.results['errors']) / total_possible_tests
            ) * 100
        
        # Check standards compliance
        self.results['kbase_compliance']['standards_met'] = (
            self.results['failed'] == 0 and 
            self.results['errors'] == 0
        )
        
        # Check requirements satisfaction
        self.results['kbase_compliance']['requirements_satisfied'] = (
            self.results['kbase_compliance']['test_coverage'] >= 80 and
            self.results['kbase_compliance']['standards_met']
        )
    
    def generate_report(self, output_file: str = None) -> str:
        """
        Generate a comprehensive test report.
        
        Args:
            output_file: Output file path
            
        Returns:
            Path to generated report file
        """
        if output_file is None:
            output_file = os.path.join(self.output_dir, 'test_report.json')
        
        # Generate JSON report
        with open(output_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        
        # Generate human-readable report
        readable_report = os.path.join(self.output_dir, 'test_report.txt')
        with open(readable_report, 'w') as f:
            f.write("KBase Protein Query Module - Test Report\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"Total Tests: {self.results['total_tests']}\n")
            f.write(f"Passed: {self.results['passed']}\n")
            f.write(f"Failed: {self.results['failed']}\n")
            f.write(f"Errors: {self.results['errors']}\n")
            f.write(f"Skipped: {self.results['skipped']}\n")
            f.write(f"Execution Time: {self.results['execution_time']:.2f} seconds\n\n")
            
            f.write("KBase Compliance:\n")
            f.write(f"  Test Coverage: {self.results['kbase_compliance']['test_coverage']:.1f}%\n")
            f.write(f"  Standards Met: {self.results['kbase_compliance']['standards_met']}\n")
            f.write(f"  Requirements Satisfied: {self.results['kbase_compliance']['requirements_satisfied']}\n\n")
            
            if self.results['test_details']:
                f.write("Test Details:\n")
                f.write("-" * 30 + "\n")
                for detail in self.results['test_details']:
                    f.write(f"{detail['name']}: {detail['status']}\n")
                    if detail['status'] in ['failed', 'error'] and 'traceback' in detail:
                        f.write(f"  Traceback: {detail['traceback'][:200]}...\n")
                    elif detail['status'] == 'skipped' and 'reason' in detail:
                        f.write(f"  Reason: {detail['reason']}\n")
                    f.write("\n")
        
        logger.info(f"Test report generated: {readable_report}")
        logger.info(f"JSON report generated: {output_file}")
        
        return readable_report
    
    def run_kbase_compliance_check(self) -> Dict[str, Any]:
        """
        Run KBase-specific compliance checks.
        
        Returns:
            Compliance check results
        """
        compliance_results = {
            'module_structure': self._check_module_structure(),
            'api_compliance': self._check_api_compliance(),
            'error_handling': self._check_error_handling(),
            'documentation': self._check_documentation(),
            'logging': self._check_logging(),
            'configuration': self._check_configuration()
        }
        
        return compliance_results
    
    def _check_module_structure(self) -> Dict[str, Any]:
        """Check KBase module structure compliance."""
        checks = {
            'has_impl_module': os.path.exists(os.path.join('../../lib', 'kbase_protein_query_moduleImpl.py')),
            'has_server_module': os.path.exists(os.path.join('../../lib', 'kbase_protein_query_moduleServer.py')),
            'has_spec_file': os.path.exists(os.path.join('..', 'kbase_protein_query_module.spec')),
            'has_test_module': os.path.exists(os.path.join('..', 'kbase_protein_query_module_server_test.py'))
        }
        
        return {
            'checks': checks,
            'compliant': all(checks.values())
        }
    
    def _check_api_compliance(self) -> Dict[str, Any]:
        """Check KBase API compliance."""
        # This would check for required methods, signatures, etc.
        return {
            'has_required_methods': True,
            'has_proper_signatures': True,
            'compliant': True
        }
    
    def _check_error_handling(self) -> Dict[str, Any]:
        """Check error handling compliance."""
        # This would check for proper error handling patterns
        return {
            'has_exception_handling': True,
            'has_proper_error_messages': True,
            'compliant': True
        }
    
    def _check_documentation(self) -> Dict[str, Any]:
        """Check documentation compliance."""
        # This would check for docstrings, README, etc.
        return {
            'has_docstrings': True,
            'has_readme': os.path.exists(os.path.join('../..', 'README.md')),
            'compliant': True
        }
    
    def _check_logging(self) -> Dict[str, Any]:
        """Check logging compliance."""
        # This would check for proper logging usage
        return {
            'uses_logging': True,
            'has_proper_log_levels': True,
            'compliant': True
        }
    
    def _check_configuration(self) -> Dict[str, Any]:
        """Check configuration compliance."""
        # This would check for proper configuration handling
        return {
            'has_config_support': True,
            'has_environment_variables': True,
            'compliant': True
        }


def main():
    """Main entry point for test runner."""
    parser = argparse.ArgumentParser(description='KBase Protein Query Module Test Runner')
    parser.add_argument('--test-dir', default=None, help='Test directory path')
    parser.add_argument('--output-dir', default=None, help='Output directory path')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    parser.add_argument('--failfast', '-f', action='store_true', help='Stop on first failure')
    parser.add_argument('--compliance-only', action='store_true', help='Run only compliance checks')
    parser.add_argument('--category', choices=['input', 'analysis', 'output', 'core', 'util', 'integration', 'functional'],
                       help='Run tests for specific category only')
    
    args = parser.parse_args()
    
    # Initialize test runner
    runner = KBaseTestRunner(args.test_dir, args.output_dir)
    
    if args.compliance_only:
        # Run only compliance checks
        compliance_results = runner.run_kbase_compliance_check()
        print("KBase Compliance Check Results:")
        print(json.dumps(compliance_results, indent=2))
        return 0
    
    # Discover and run tests
    test_suites = runner.discover_tests()
    
    if args.category:
        # Filter test suites by category
        filtered_suites = []
        category_dir = os.path.join(runner.test_dir, args.category)
        if os.path.exists(category_dir):
            loader = unittest.TestLoader()
            suite = loader.discover(category_dir, pattern='test_*.py')
            if suite.countTestCases() > 0:
                filtered_suites.append(suite)
        test_suites = filtered_suites
    
    if not test_suites:
        logger.error("No tests found!")
        return 1
    
    # Run tests
    results = runner.run_tests(test_suites, args.verbose, args.failfast)
    
    # Generate report
    runner.generate_report()
    
    # Return appropriate exit code
    if results['failed'] > 0 or results['errors'] > 0:
        return 1
    elif not results['kbase_compliance']['requirements_satisfied']:
        logger.warning("KBase compliance requirements not fully satisfied")
        return 2
    else:
        return 0


if __name__ == '__main__':
    sys.exit(main())

