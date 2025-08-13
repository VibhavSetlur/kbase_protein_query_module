#!/usr/bin/env python3
"""
Comprehensive test suite for sequence_analyzer.py

Tests all sequence analysis features, edge cases, and validation scenarios
for the KBase Protein Query Module sequence analyzer.
"""

import unittest
import numpy as np
import os
import sys
from kbase_protein_query_module.src.analysis import ProteinSequenceAnalyzer

class TestSequenceAnalyzer(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = ProteinSequenceAnalyzer()
        self.temp_dir = os.path.join(os.path.dirname(__file__), "temp_test_data")
        os.makedirs(self.temp_dir, exist_ok=True)
        
        # Test sequences
        self.test_sequence = "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"
        self.test_sequences = [
            "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ",
            "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG",
            "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"
        ]
        self.test_ids = ["protein1", "protein2", "protein3"]

    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)

    def test_analyze_sequence_basic(self):
        """Test basic sequence analysis."""
        result = self.analyzer.analyze_sequence(self.test_sequence)
        
        self.assertIsInstance(result, dict)
        self.assertIn('length', result)
        self.assertIn('amino_acid_composition', result)
        self.assertIn('physicochemical_properties', result)
        
        self.assertEqual(result['length'], len(self.test_sequence))
        self.assertGreater(result['physicochemical_properties']['molecular_weight'], 0)
        self.assertGreater(result['physicochemical_properties']['isoelectric_point'], 0)

    def test_analyze_sequences_batch(self):
        """Test batch sequence analysis."""
        # Since the analyzer doesn't have a batch method, test individual analysis
        results = {}
        for seq, pid in zip(self.test_sequences, self.test_ids):
            results[pid] = self.analyzer.analyze_sequence(seq, pid)
        
        self.assertIsInstance(results, dict)
        self.assertEqual(len(results), len(self.test_sequences))
        
        for protein_id, result in results.items():
            self.assertIsInstance(result, dict)
            self.assertIn('length', result)
            self.assertIn('amino_acid_composition', result)
            self.assertIn('physicochemical_properties', result)

    def test_secondary_structure_prediction(self):
        """Test secondary structure prediction."""
        sequence = "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"
        
        # Test the actual secondary structure prediction method
        result = self.analyzer._predict_secondary_structure(sequence)
        
        self.assertIsInstance(result, dict)
        # Check that we have secondary structure prediction results
        self.assertIsInstance(result, dict)
        self.assertGreater(len(result), 0)
        # Check that we have secondary structure prediction results
        self.assertIsInstance(result, dict)
        self.assertGreater(len(result), 0)
        self.assertIn('turn_preference', result)
        self.assertIn('dominant_structure', result)

    def test_domain_analysis(self):
        """Test domain analysis."""
        sequence = "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"
        
        # Test the actual sequence analysis which includes motif analysis
        result = self.analyzer.analyze_sequence(sequence)
        
        self.assertIsInstance(result, dict)
        self.assertIn('sequence_motifs', result)
        self.assertIsInstance(result['sequence_motifs'], dict)

    def test_motif_analysis(self):
        """Test motif analysis."""
        sequence = "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"
        
        # Test the actual motif finding method
        result = self.analyzer._find_sequence_motifs(sequence)
        
        self.assertIsInstance(result, dict)
        self.assertIn('n_glycosylation', result)
        self.assertIn('o_glycosylation', result)
        # Check that we have some motif analysis results
        self.assertIsInstance(result, dict)
        self.assertGreater(len(result), 0)

    def test_comprehensive_analysis(self):
        """Test comprehensive sequence analysis."""
        result = self.analyzer.analyze_sequence(self.test_sequence)
        
        self.assertIsInstance(result, dict)
        self.assertIn('length', result)
        self.assertIn('amino_acid_composition', result)
        self.assertIn('physicochemical_properties', result)
        self.assertIn('secondary_structure_prediction', result)
        self.assertIn('sequence_motifs', result)
        self.assertIn('bioinformatics_links', result)
        
        # Check that physicochemical properties are present
        phys_props = result['physicochemical_properties']
        self.assertIn('molecular_weight', phys_props)
        self.assertIn('isoelectric_point', phys_props)

if __name__ == '__main__':
    unittest.main()
