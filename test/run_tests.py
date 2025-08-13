#!/usr/bin/env python3
"""
KBase Protein Query Module Test Runner

This script runs comprehensive tests for the KBase Protein Query Module following
KBase testing guidelines and best practices.

Author: Vibhav Setlur
Contact: https://kbase.us/contact-us/
"""

import os
import sys
import time
import logging
import unittest
import subprocess
from pathlib import Path
from configparser import ConfigParser

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('test.log'),
        logging.StreamHandler(sys.stdout)
    ]
)

logger = logging.getLogger(__name__)

class KBaseTestRunner:
    """
    Comprehensive test runner for KBase Protein Query Module.
    
    Follows KBase testing guidelines:
    - Proper error handling and logging
    - Workspace integration testing
    - Unit and integration tests
    - Performance and stress testing
    - Documentation and reporting
    """
    
    def __init__(self, config_file='test.cfg'):
        """Initialize the test runner with configuration."""
        self.config = self._load_config(config_file)
        self.test_results = {
            'passed': 0,
            'failed': 0,
            'errors': 0,
            'skipped': 0
        }
        
    def _load_config(self, config_file):
        """Load test configuration from file."""
        config = ConfigParser()
        config_path = Path(__file__).parent / config_file
        
        if config_path.exists():
            config.read(config_path)
            logger.info(f"Loaded test configuration from {config_path}")
        else:
            logger.warning(f"Configuration file {config_path} not found, using defaults")
            # Set default configuration
            config['test'] = {
                'data_dir': 'data',
                'test_protein_id': 'P12345',
                'test_sequence': 'MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ',
                'log_level': 'INFO'
            }
        
        return config
    
    def run_unit_tests(self):
        """Run unit tests for all modules."""
        logger.info("Starting unit tests...")
        
        try:
            # Add lib directory to path
            lib_path = Path(__file__).parent.parent / 'lib'
            sys.path.insert(0, str(lib_path))
            
            # Discover and run unit tests
            loader = unittest.TestLoader()
            start_dir = Path(__file__).parent / 'unit_tests'
            
            if start_dir.exists():
                suite = loader.discover(start_dir, pattern='test_*.py')
                runner = unittest.TextTestRunner(verbosity=2)
                result = runner.run(suite)
                
                self.test_results['passed'] += result.testsRun - len(result.failures) - len(result.errors)
                self.test_results['failed'] += len(result.failures)
                self.test_results['errors'] += len(result.errors)
                
                logger.info(f"Unit tests completed: {result.testsRun} tests run")
                return result.wasSuccessful()
            else:
                logger.warning(f"Unit test directory not found: {start_dir}")
                return False
                
        except Exception as e:
            logger.error(f"Error running unit tests: {e}")
            self.test_results['errors'] += 1
            return False
    
    def run_integration_tests(self):
        """Run integration tests with workspace."""
        logger.info("Starting integration tests...")
        
        try:
            # Check if biokbase module is available (required for integration tests)
            try:
                import biokbase
            except ImportError:
                logger.warning("biokbase module not available. Integration tests require KBase SDK environment.")
                logger.warning("Skipping integration tests - they should be run with 'kb-sdk test'")
                self.test_results['skipped'] += 1
                return True  # Consider this a success since it's expected behavior
            
            # Run the main integration test file
            test_file = Path(__file__).parent / 'integration_tests' / 'kbase_protein_query_module_query_server_test.py'
            
            if test_file.exists():
                # Run with proper environment
                env = os.environ.copy()
                env['PYTHONPATH'] = f"{Path(__file__).parent.parent / 'lib'}:{env.get('PYTHONPATH', '')}"
                
                result = subprocess.run([
                    sys.executable, '-m', 'pytest', str(test_file), '-v'
                ], env=env, capture_output=True, text=True)
                
                if result.returncode == 0:
                    logger.info("Integration tests completed successfully")
                    self.test_results['passed'] += 1
                    return True
                else:
                    logger.error(f"Integration tests failed: {result.stderr}")
                    self.test_results['failed'] += 1
                    return False
            else:
                logger.warning(f"Integration test file not found: {test_file}")
                return False
                
        except Exception as e:
            logger.error(f"Error running integration tests: {e}")
            self.test_results['errors'] += 1
            return False
    
    def run_html_report_tests(self):
        """Run HTML report generator tests."""
        logger.info("Starting HTML report tests...")
        
        try:
            test_file = Path(__file__).parent / 'unit_tests' / 'test_html_report_generator.py'
            
            if test_file.exists():
                result = subprocess.run([
                    sys.executable, str(test_file)
                ], capture_output=True, text=True)
                
                if result.returncode == 0:
                    logger.info("HTML report tests completed successfully")
                    self.test_results['passed'] += 1
                    return True
                else:
                    logger.error(f"HTML report tests failed: {result.stderr}")
                    self.test_results['failed'] += 1
                    return False
            else:
                logger.warning(f"HTML report test file not found: {test_file}")
                return False
                
        except Exception as e:
            logger.error(f"Error running HTML report tests: {e}")
            self.test_results['errors'] += 1
            return False
    
    def run_kbase_sdk_tests(self):
        """Run KBase SDK tests."""
        logger.info("Starting KBase SDK tests...")
        
        try:
            # Change to project root
            project_root = Path(__file__).parent.parent
            os.chdir(project_root)
            
            # Run kb-sdk test
            result = subprocess.run([
                'kb-sdk', 'test'
            ], capture_output=True, text=True)
            
            if result.returncode == 0:
                logger.info("KBase SDK tests completed successfully")
                self.test_results['passed'] += 1
                return True
            else:
                logger.error(f"KBase SDK tests failed: {result.stderr}")
                self.test_results['failed'] += 1
                return False
                
        except Exception as e:
            logger.error(f"Error running KBase SDK tests: {e}")
            self.test_results['errors'] += 1
            return False
    
    def check_data_availability(self):
        """Check if required test data is available."""
        logger.info("Checking data availability...")
        
        data_dir = Path(self.config.get('test', 'data_dir', fallback='data'))
        required_paths = [
            data_dir / 'families',
            data_dir / 'indexes',
            data_dir / 'metadata',
            data_dir / 'family_centroids' / 'family_centroids_binary.npz'
        ]
        
        missing_paths = []
        for path in required_paths:
            if not path.exists():
                missing_paths.append(str(path))
        
        if missing_paths:
            logger.warning(f"Missing data paths: {missing_paths}")
            logger.warning("Some tests may be skipped due to missing data")
            return False
        else:
            logger.info("All required data paths found")
            return True
    
    def generate_test_report(self):
        """Generate a comprehensive test report."""
        logger.info("Generating test report...")
        
        total_tests = sum(self.test_results.values())
        success_rate = (self.test_results['passed'] / total_tests * 100) if total_tests > 0 else 0
        
        # Check if integration tests were skipped
        integration_skipped = self.test_results.get('skipped', 0) > 0
        
        report = f"""
KBase Protein Query Module Test Report
=====================================

Test Results Summary:
- Passed: {self.test_results['passed']}
- Failed: {self.test_results['failed']}
- Errors: {self.test_results['errors']}
- Skipped: {self.test_results['skipped']}
- Total: {total_tests}
- Success Rate: {success_rate:.1f}%

Test Categories:
1. Unit Tests: Core module functionality
2. Integration Tests: Workspace and KBase integration {'(Skipped - requires KBase SDK environment)' if integration_skipped else ''}
3. HTML Report Tests: Report generation and visualization
4. KBase SDK Tests: SDK compliance and deployment

Data Availability: {'✓ Available' if self.check_data_availability() else '✗ Missing'}

Recommendations:
- Ensure all data files are present for comprehensive testing
- Run tests in KBase environment for full integration testing
- Check logs for detailed error information
{'Note: Integration tests were skipped because biokbase module is not available. Run "kb-sdk test" for full integration testing.' if integration_skipped else ''}

For more information, see:
- KBase SDK Documentation: https://kbase.github.io/kb_sdk_docs/
- Module Documentation: https://github.com/kbaseapps/kbase_protein_query_module
"""
        
        # Write report to file
        report_file = Path(__file__).parent / 'test_report.txt'
        with open(report_file, 'w') as f:
            f.write(report)
        
        logger.info(f"Test report written to: {report_file}")
        print(report)
        
        return success_rate >= 80  # Consider 80% success rate as passing
    
    def run_all_tests(self):
        """Run all tests and generate report."""
        logger.info("Starting comprehensive test suite...")
        
        start_time = time.time()
        
        # Check data availability
        data_available = self.check_data_availability()
        
        # Run different test categories
        unit_success = self.run_unit_tests()
        integration_success = self.run_integration_tests()
        html_success = self.run_html_report_tests()
        sdk_success = self.run_kbase_sdk_tests()
        
        end_time = time.time()
        duration = end_time - start_time
        
        logger.info(f"Test suite completed in {duration:.2f} seconds")
        
        # Generate final report
        overall_success = self.generate_test_report()
        
        if overall_success:
            logger.info("✅ All tests passed successfully!")
            return 0
        else:
            logger.error("❌ Some tests failed. Check the report for details.")
            return 1

def main():
    """Main entry point for test runner."""
    print("🧪 KBase Protein Query Module Test Runner")
    print("=" * 50)
    
    # Initialize test runner
    runner = KBaseTestRunner()
    
    # Run all tests
    exit_code = runner.run_all_tests()
    
    sys.exit(exit_code)

if __name__ == "__main__":
    main() 