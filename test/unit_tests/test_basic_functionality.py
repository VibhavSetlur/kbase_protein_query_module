"""
Basic functionality tests for KBase Protein Query Module
"""

import unittest
import sys
import os

# Add the lib directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'lib'))

class TestBasicFunctionality(unittest.TestCase):
    """Test basic module functionality"""
    
    def test_module_import(self):
        """Test that the module can be imported"""
        try:
            from kbase_protein_query_module.kbase_protein_query_moduleImpl import kbase_protein_query_module
            self.assertTrue(True, "Module imported successfully")
        except ImportError as e:
            self.fail(f"Failed to import module: {e}")
    
    def test_module_initialization(self):
        """Test that the module can be initialized"""
        try:
            from kbase_protein_query_module.kbase_protein_query_moduleImpl import kbase_protein_query_module
            config = {'scratch': '/tmp'}
            module = kbase_protein_query_module(config)
            self.assertIsNotNone(module)
            self.assertEqual(module.VERSION, "2.0.0")
        except Exception as e:
            self.fail(f"Failed to initialize module: {e}")
    
    def test_status_method(self):
        """Test the status method"""
        try:
            from kbase_protein_query_module.kbase_protein_query_moduleImpl import kbase_protein_query_module
            config = {'scratch': '/tmp'}
            module = kbase_protein_query_module(config)
            
            ctx = {}
            status_result = module.status(ctx)
            
            self.assertIsInstance(status_result, list)
            self.assertEqual(len(status_result), 1)
            self.assertEqual(status_result[0]['state'], 'OK')
            self.assertEqual(status_result[0]['version'], '2.0.0')
        except Exception as e:
            self.fail(f"Status method failed: {e}")

if __name__ == '__main__':
    unittest.main()
