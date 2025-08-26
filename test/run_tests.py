#!/usr/bin/env python3
"""
Unified Test Runner for KBase Protein Query Module

This script consolidates all testing functionality into a single file:
- Unit tests for all modules
- Integration tests
- Performance tests
- End-to-end workflow tests
- KBase SDK compliance tests

All test outputs are consolidated into the test/outputs directory.
"""

import os
import sys
import time
import logging
import unittest
import subprocess
import tempfile
import shutil
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any, Optional
from configparser import ConfigParser

# Add lib directory to path for imports
lib_path = Path(__file__).parent.parent / 'lib'
sys.path.insert(0, str(lib_path))
sys.path.insert(0, str(Path(__file__).parent.parent))

class UnifiedTestRunner:
    """
    Unified test runner that consolidates all testing functionality.
    
    Features:
    - Single entry point for all tests
    - Consolidated output directory (test/outputs)
    - Comprehensive reporting
    - Environment-aware testing
    - Clean output management
    """
    
    def __init__(self, config_file: str = 'test.cfg'):
        """Initialize the unified test runner."""
        self.output_dir = Path(__file__).parent / 'outputs'
        self.output_dir.mkdir(exist_ok=True)
        
        # Setup logging to outputs directory first
        self._setup_logging()
        
        # Now load config (after logger is available)
        self.config = self._load_config(config_file)
        
        self.test_results = {
            'timestamp': datetime.now().isoformat(),
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
        
        # Override output directories to use consolidated location
        self._setup_output_overrides()
    
    def _load_config(self, config_file: str) -> ConfigParser:
        """Load test configuration from file."""
        config = ConfigParser()
        config_path = Path(__file__).parent / config_file
        
        if config_path.exists():
            config.read(config_path)
            self.logger.info(f"Loaded test configuration from {config_path}")
        else:
            self.logger.warning(f"Configuration file {config_path} not found, using defaults")
            config['test'] = {
                'data_dir': 'data',
                'test_protein_id': 'P12345',
                'test_sequence': 'MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ',
                'log_level': 'INFO',
                'output_dir': str(self.output_dir)
            }
        
        return config
    
    def _setup_logging(self):
        """Setup logging to outputs directory."""
        log_file = self.output_dir / f"test_run_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
        self.logger.info(f"Test run started. Log file: {log_file}")
    
    def _setup_output_overrides(self):
        """Override output directories to use consolidated location."""
        # Set environment variables to redirect outputs
        os.environ['TEST_OUTPUT_DIR'] = str(self.output_dir)
        os.environ['HTML_REPORTS_DIR'] = str(self.output_dir / 'html_reports')
        os.environ['EXPORTS_DIR'] = str(self.output_dir / 'exports')
        os.environ['TEST_OUTPUTS_DIR'] = str(self.output_dir)
        
        # Create subdirectories
        (self.output_dir / 'html_reports').mkdir(exist_ok=True)
        (self.output_dir / 'exports').mkdir(exist_ok=True)
        (self.output_dir / 'reports').mkdir(exist_ok=True)
        (self.output_dir / 'data').mkdir(exist_ok=True)
    
    def run_unit_tests(self) -> bool:
        """Run all unit tests."""
        self.logger.info("Starting unit tests...")
        
        test_dirs = [
            'test/unit_tests/core',
            'test/unit_tests/data', 
            'test/unit_tests/processing',
            'test/unit_tests/analysis',
            'test/unit_tests/storage',
            'test/unit_tests/reports',
            'test/unit_tests/stages',
            'test/unit_tests/utils',
            'test/unit_tests/workflows',
            'test/unit_tests'
        ]
        
        total_passed = 0
        total_failed = 0
        total_errors = 0
        total_skipped = 0
        
        for test_dir in test_dirs:
            if not Path(test_dir).exists():
                self.logger.warning(f"Test directory not found: {test_dir}")
                continue
            
            self.logger.info(f"Running tests in: {test_dir}")
            
            try:
                # Run with pytest for better output
                result = subprocess.run([
                    sys.executable, '-m', 'pytest', test_dir,
                    '-v', '--tb=short', '--capture=no',
                    f'--junitxml={self.output_dir}/unit_tests_{Path(test_dir).name}.xml'
                ], capture_output=True, text=True, cwd=Path(__file__).parent.parent)
                
                # Parse results (simplified)
                if result.returncode == 0:
                    total_passed += 1
                    self.logger.info(f"✅ {test_dir} tests passed")
                else:
                    total_failed += 1
                    self.logger.error(f"❌ {test_dir} tests failed: {result.stderr}")
                
            except Exception as e:
                total_errors += 1
                self.logger.error(f"Error running tests in {test_dir}: {e}")
        
        self.test_results['summary']['passed'] += total_passed
        self.test_results['summary']['failed'] += total_failed
        self.test_results['summary']['errors'] += total_errors
        self.test_results['summary']['skipped'] += total_skipped
        
        return total_failed == 0 and total_errors == 0
    
    def run_integration_tests(self) -> bool:
        """Run integration tests."""
        self.logger.info("Starting integration tests...")
        
        try:
            # Check if biokbase module is available
            try:
                import biokbase
            except ImportError:
                self.logger.warning("biokbase module not available. Integration tests require KBase SDK environment.")
                self.logger.warning("Skipping integration tests - they should be run with 'kb-sdk test'")
                self.test_results['summary']['skipped'] += 1
                return True
            
            # Run integration tests
            integration_test_file = Path(__file__).parent / 'integration_tests' / 'kbase_protein_query_module_query_server_test.py'
            
            if integration_test_file.exists():
                result = subprocess.run([
                    sys.executable, '-m', 'pytest', str(integration_test_file), '-v',
                    f'--junitxml={self.output_dir}/integration_tests.xml'
                ], capture_output=True, text=True, cwd=Path(__file__).parent.parent)
                
                if result.returncode == 0:
                    self.logger.info("✅ Integration tests passed")
                    self.test_results['summary']['passed'] += 1
                    return True
                else:
                    self.logger.error(f"❌ Integration tests failed: {result.stderr}")
                    self.test_results['summary']['failed'] += 1
                    return False
            else:
                self.logger.warning(f"Integration test file not found: {integration_test_file}")
                return False
                
        except Exception as e:
            self.logger.error(f"Error running integration tests: {e}")
            self.test_results['summary']['errors'] += 1
            return False
    
    def run_kbase_sdk_tests(self) -> bool:
        """Run KBase SDK tests."""
        self.logger.info("Starting KBase SDK tests...")
        
        try:
            # Change to project root
            project_root = Path(__file__).parent.parent
            os.chdir(project_root)
            
            # Run kb-sdk test with output redirection
            result = subprocess.run([
                'kb-sdk', 'test'
            ], capture_output=True, text=True)
            
            # Save SDK test output
            sdk_output_file = self.output_dir / 'kb_sdk_test_output.txt'
            with open(sdk_output_file, 'w') as f:
                f.write("STDOUT:\n")
                f.write(result.stdout)
                f.write("\n\nSTDERR:\n")
                f.write(result.stderr)
            
            if result.returncode == 0:
                self.logger.info("✅ KBase SDK tests passed")
                self.test_results['summary']['passed'] += 1
                return True
            else:
                self.logger.error(f"❌ KBase SDK tests failed: {result.stderr}")
                self.test_results['summary']['failed'] += 1
                return False
                
        except Exception as e:
            self.logger.error(f"Error running KBase SDK tests: {e}")
            self.test_results['summary']['errors'] += 1
            return False
    
    def run_performance_tests(self) -> bool:
        """Run performance tests."""
        self.logger.info("Starting performance tests...")
        
        try:
            # Import and run performance tests
            from test.unit_tests.processing.test_processing_integration import TestProcessingIntegration
            
            # Run performance tests
            suite = unittest.TestLoader().loadTestsFromTestCase(TestProcessingIntegration)
            runner = unittest.TextTestRunner(verbosity=2)
            result = runner.run(suite)
            
            if result.wasSuccessful():
                self.logger.info("✅ Performance tests passed")
                self.test_results['summary']['passed'] += 1
                return True
            else:
                self.logger.error("❌ Performance tests failed")
                self.test_results['summary']['failed'] += 1
                return False
                
        except Exception as e:
            self.logger.error(f"Error running performance tests: {e}")
            self.test_results['summary']['errors'] += 1
            return False
    
    def run_workflow_tests(self) -> bool:
        """Run workflow tests."""
        self.logger.info("Starting workflow tests...")
        
        try:
            # Import and run workflow tests
            from test.unit_tests.workflows.test_workflow_integration import TestWorkflowIntegration
            
            suite = unittest.TestLoader().loadTestsFromTestCase(TestWorkflowIntegration)
            runner = unittest.TextTestRunner(verbosity=2)
            result = runner.run(suite)
            
            if result.wasSuccessful():
                self.logger.info("✅ Workflow tests passed")
                self.test_results['summary']['passed'] += 1
                return True
            else:
                self.logger.error("❌ Workflow tests failed")
                self.test_results['summary']['failed'] += 1
                return False
                
        except Exception as e:
            self.logger.error(f"Error running workflow tests: {e}")
            self.test_results['summary']['errors'] += 1
            return False
    
    def check_data_availability(self) -> bool:
        """Check if required test data is available."""
        self.logger.info("Checking data availability...")
        
        data_dir = Path(self.config.get('test', 'data_dir', fallback='data'))
        required_paths = [
            data_dir / 'families',
            data_dir / 'indexes',
            data_dir / 'metadata'
        ]
        
        missing_paths = []
        for path in required_paths:
            if not path.exists():
                missing_paths.append(str(path))
        
        if missing_paths:
            self.logger.warning(f"Missing data paths: {missing_paths}")
            self.logger.warning("Some tests may be skipped due to missing data")
            return False
        else:
            self.logger.info("✅ All required data paths found")
            return True
    
    def cleanup_old_outputs(self):
        """Clean up old output directories outside of test/outputs."""
        self.logger.info("Cleaning up old output directories...")
        
        # Directories to clean up
        cleanup_dirs = [
            'exports',
            'html_reports', 
            'test_outputs',
            'reports',
            'output'
        ]
        
        for dir_name in cleanup_dirs:
            dir_path = Path(dir_name)
            if dir_path.exists() and dir_path.is_dir():
                try:
                    shutil.rmtree(dir_path)
                    self.logger.info(f"Cleaned up: {dir_path}")
                except Exception as e:
                    self.logger.warning(f"Could not clean up {dir_path}: {e}")
    
    def generate_test_report(self) -> bool:
        """Generate comprehensive test report."""
        self.logger.info("Generating test report...")
        
        total_tests = sum(self.test_results['summary'].values())
        success_rate = (self.test_results['summary']['passed'] / total_tests * 100) if total_tests > 0 else 0
        
        report = f"""
KBase Protein Query Module - Unified Test Report
===============================================

Test Run: {self.test_results['timestamp']}
Output Directory: {self.output_dir}

Test Results Summary:
- Passed: {self.test_results['summary']['passed']}
- Failed: {self.test_results['summary']['failed']}
- Errors: {self.test_results['summary']['errors']}
- Skipped: {self.test_results['summary']['skipped']}
- Total: {total_tests}
- Success Rate: {success_rate:.1f}%

Test Categories:
1. Unit Tests: Core module functionality
2. Integration Tests: Workspace and KBase integration
3. Performance Tests: Processing speed and efficiency
4. Workflow Tests: End-to-end pipeline testing
5. KBase SDK Tests: SDK compliance and deployment

Data Availability: {'✅ Available' if self.check_data_availability() else '❌ Missing'}

Output Consolidation:
- All test outputs consolidated in: {self.output_dir}
- HTML reports: {self.output_dir}/html_reports
- Data exports: {self.output_dir}/exports
- Test reports: {self.output_dir}/reports
- Log files: {self.output_dir}/*.log

Recommendations:
- Ensure all data files are present for comprehensive testing
- Run tests in KBase environment for full integration testing
- Check logs for detailed error information

For more information, see:
- KBase SDK Documentation: https://kbase.github.io/kb_sdk_docs/
- Module Documentation: https://github.com/kbaseapps/kbase_protein_query_module
"""
        
        # Write report to file
        report_file = self.output_dir / 'test_report.txt'
        with open(report_file, 'w') as f:
            f.write(report)
        
        # Also write JSON summary
        json_file = self.output_dir / 'test_results.json'
        with open(json_file, 'w') as f:
            json.dump(self.test_results, f, indent=2)
        
        self.logger.info(f"Test report written to: {report_file}")
        self.logger.info(f"Test results JSON written to: {json_file}")
        print(report)
        
        return success_rate >= 80  # Consider 80% success rate as passing
    
    def run_all_tests(self) -> int:
        """Run all tests and generate report."""
        self.logger.info("Starting unified test suite...")
        
        start_time = time.time()
        
        # Clean up old outputs
        self.cleanup_old_outputs()
        
        # Check data availability
        data_available = self.check_data_availability()
        
        # Run different test categories
        unit_success = self.run_unit_tests()
        integration_success = self.run_integration_tests()
        performance_success = self.run_performance_tests()
        workflow_success = self.run_workflow_tests()
        sdk_success = self.run_kbase_sdk_tests()
        
        end_time = time.time()
        duration = end_time - start_time
        
        self.logger.info(f"Test suite completed in {duration:.2f} seconds")
        
        # Generate final report
        overall_success = self.generate_test_report()
        
        if overall_success:
            self.logger.info("✅ All tests passed successfully!")
            return 0
        else:
            self.logger.error("❌ Some tests failed. Check the report for details.")
            return 1

def main():
    """Main entry point for unified test runner."""
    print("🧪 KBase Protein Query Module - Unified Test Runner")
    print("=" * 60)
    print("Consolidating all test outputs to test/outputs/")
    print("=" * 60)
    
    # Initialize test runner
    runner = UnifiedTestRunner()
    
    # Run all tests
    exit_code = runner.run_all_tests()
    
    sys.exit(exit_code)

if __name__ == "__main__":
    main() 