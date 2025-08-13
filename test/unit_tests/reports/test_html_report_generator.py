#!/usr/bin/env python3
"""
Test script for HTML Report Generator with Sequence Analyzer Integration
"""

import unittest
import numpy as np
import os
import sys
import tempfile
from kbase_protein_query_module.src.reports.html.generator import HTMLReportGenerator

# Add the lib directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))

def test_html_report_generator():
    """Test the HTML report generator with sequence analysis."""
    
    try:
        from kbase_protein_query_module.src.reports.html.generator import HTMLReportGenerator
        from kbase_protein_query_module.src.analysis.sequence_analyzer import ProteinSequenceAnalyzer
        
        print("✅ Successfully imported HTML report generator and sequence analyzer")
        
        # Create a temporary directory for testing
        with tempfile.TemporaryDirectory() as temp_dir:
            print(f"📁 Created temporary directory: {temp_dir}")
            
            # Initialize the HTML report generator
            report_generator = HTMLReportGenerator(output_dir=temp_dir)
            
            # Create sample pipeline results
            sample_sequence = "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"
            sample_protein_id = "TEST_PROTEIN_001"
            
            pipeline_results = {
                'protein_existence': {
                    'protein_id': sample_protein_id,
                    'exists': True,
                    'family_id': 'family_0',
                    'metadata': {'source': 'test'}
                },
                'embedding_generation': {
                    'embedding_dim': 1280,
                    'embedding_norm': 1.0,
                    'sequence_length': len(sample_sequence),
                    'sequence': sample_sequence,
                    'model_name': 'esm2_t6_8M_UR50D'
                },
                'family_assignment': {
                    'family_id': 'family_0',
                    'confidence': 0.85,
                    'eigenprotein_id': 'eigen_001'
                },
                'similarity_search': {
                    'matches': [
                        {'protein_id': 'PROT_001', 'similarity': 0.95},
                        {'protein_id': 'PROT_002', 'similarity': 0.92},
                        {'protein_id': 'PROT_003', 'similarity': 0.89}
                    ],
                    'similarity_stats': {
                        'max': 0.95,
                        'min': 0.89,
                        'mean': 0.92
                    },
                    'family_id': 'family_0',
                    'top_n': 3
                },
                'network_building': {
                    'network_nodes': 4,
                    'network_edges': 6,
                    'network_density': 1.0
                }
            }
            
            print("📊 Created sample pipeline results")
            
            # Generate comprehensive HTML report
            try:
                report_result = report_generator.generate_comprehensive_report(
                    pipeline_results, sample_protein_id, sample_sequence
                )
                
                print("✅ Successfully generated comprehensive HTML report")
                print(f"📄 HTML report path: {report_result['html_path']}")
                print(f"🕸️ Network visualization: {'Available' if report_result['network_viz_path'] else 'Not available'}")
                print(f"📊 Dashboard data: {'Available' if report_result['dashboard_data'] else 'Not available'}")
                
                # Check if HTML file was created
                if os.path.exists(report_result['html_path']):
                    print("✅ HTML file was created successfully")
                    
                    # Read and check HTML content
                    with open(report_result['html_path'], 'r', encoding='utf-8') as f:
                        html_content = f.read()
                    
                    # Check for key elements in the HTML
                    required_elements = [
                        'Comprehensive Protein Network Analysis Report',
                        'Sequence Analysis',
                        'Network Visualization',
                        'Statistics',
                        'Bioinformatics'
                    ]
                    
                    missing_elements = []
                    for element in required_elements:
                        if element not in html_content:
                            missing_elements.append(element)
                    
                    if missing_elements:
                        print(f"⚠️ Missing HTML elements: {missing_elements}")
                    else:
                        print("✅ All required HTML elements found")
                        
                else:
                    print("❌ HTML file was not created")
                    raise AssertionError("HTML file was not created")
                    
            except Exception as e:
                print(f"❌ Error generating HTML report: {e}")
                import traceback
                traceback.print_exc()
                raise
            
            # Test sequence analyzer directly
            try:
                sequence_analyzer = ProteinSequenceAnalyzer()
                analysis = sequence_analyzer.analyze_sequence(sample_sequence, sample_protein_id)
                
                print("✅ Sequence analysis completed successfully")
                print(f"📊 Analysis includes: {list(analysis.keys())}")
                
                # Check for key analysis components
                required_analysis = [
                    'amino_acid_composition',
                    'physicochemical_properties',
                    'secondary_structure_prediction',
                    'sequence_motifs',
                    'bioinformatics_links'
                ]
                
                missing_analysis = []
                for component in required_analysis:
                    if component not in analysis:
                        missing_analysis.append(component)
                
                if missing_analysis:
                    print(f"⚠️ Missing analysis components: {missing_analysis}")
                    raise AssertionError(f"Missing analysis components: {missing_analysis}")
                else:
                    print("✅ All required analysis components found")
                    
            except Exception as e:
                print(f"❌ Error in sequence analysis: {e}")
                import traceback
                traceback.print_exc()
                raise
            
            print("\n🎉 All tests passed! HTML report generator is working correctly.")
            
    except ImportError as e:
        print(f"❌ Import error: {e}")
        raise
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        raise

if __name__ == "__main__":
    print("🧪 Testing HTML Report Generator with Sequence Analyzer Integration")
    print("=" * 70)
    
    test_html_report_generator() 