#!/usr/bin/env python3
"""
Test script to verify network analysis output generation with proper dependencies.
"""

import sys
import os
import tempfile
import json

# Add lib to path
sys.path.insert(0, '/scratch/vsetlur/kbase_protein_query_module/lib')

from kbase_protein_query_module.src.analysis.network_analysis.network_analysis import NetworkAnalysis
from kbase_protein_query_module.src.output.analysis.network_analysis.output import NetworkAnalysisOutput

def test_network_analysis():
    """Test network analysis with mock data."""
    print("Testing Network Analysis with mock data...")
    
    # Create temporary output directory
    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"Using temp directory: {temp_dir}")
        
        # Create mock analysis data with the correct structure
        import numpy as np
        import pandas as pd
        
        mock_data = {
            'proteins': [
                {'protein_id': 'P12345', 'sequence': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'},
                {'protein_id': 'P67890', 'sequence': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'},
                {'protein_id': 'P11111', 'sequence': 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'}
            ],
            'similarity_results': [
                {'protein1': 'P12345', 'protein2': 'P67890', 'similarity': 0.8},
                {'protein1': 'P12345', 'protein2': 'P11111', 'similarity': 0.6},
                {'protein1': 'P67890', 'protein2': 'P11111', 'similarity': 0.7}
            ],
            'embeddings': np.random.rand(3, 320),  # Mock embeddings
            'protein_ids': ['P12345', 'P67890', 'P11111'],
            'metadata_df': pd.DataFrame({
                'Entry': ['P12345', 'P67890', 'P11111'],
                'protein_id': ['P12345', 'P67890', 'P11111'],
                'organism': ['E.coli', 'E.coli', 'E.coli'],
                'length': [100, 100, 100]
            }),
            'query_embedding': np.random.rand(320),
            'query_protein_id': 'P12345',
            'output_dir': temp_dir
        }
        
        # Test network analysis
        try:
            network_analysis = NetworkAnalysis()
            result = network_analysis.run(mock_data)
            
            if result and result.success:
                print("✓ Network analysis completed successfully")
                analysis_data = result.data.get('network_analysis', {})
                print(f"  - Analysis data keys: {list(analysis_data.keys())}")
                
                # Test output generation
                output_manager = NetworkAnalysisOutput()
                output_result = output_manager.generate_outputs(analysis_data, 'network_analysis')
                
                if output_result and output_result.success:
                    print("✓ Network analysis outputs generated successfully")
                    print(f"  - Output files: {output_result.output_files}")
                    
                    # Check if files were actually created
                    for file_path in output_result.output_files:
                        full_path = os.path.join(temp_dir, file_path)
                        if os.path.exists(full_path):
                            file_size = os.path.getsize(full_path)
                            print(f"  - {file_path}: {file_size} bytes")
                        else:
                            print(f"  - {file_path}: FILE NOT FOUND")
                    
                    return True
                else:
                    print("✗ Output generation failed")
                    if output_result:
                        print(f"  - Output result: {output_result}")
                        if hasattr(output_result, 'error_message'):
                            print(f"  - Error: {output_result.error_message}")
                    return False
            else:
                print("✗ Network analysis failed")
                if result:
                    print(f"  Error: {result.error_message}")
                return False
                
        except Exception as e:
            print(f"✗ Exception during network analysis: {e}")
            import traceback
            traceback.print_exc()
            return False

def test_analysis_config():
    """Test that analysis configuration is working."""
    print("Testing analysis configuration...")
    
    try:
        from kbase_protein_query_module.src.analysis.config import get_enabled_analyses, is_analysis_enabled
        
        enabled = get_enabled_analyses()
        print(f"✓ Enabled analyses: {list(enabled.keys())}")
        
        network_enabled = is_analysis_enabled('network_analysis')
        print(f"✓ Network analysis enabled: {network_enabled}")
        
        return True
        
    except Exception as e:
        print(f"✗ Analysis config test failed: {e}")
        return False

if __name__ == "__main__":
    print("=== Network Analysis Test ===")
    
    # Test configuration first
    config_ok = test_analysis_config()
    print()
    
    if config_ok:
        # Test actual network analysis
        analysis_ok = test_network_analysis()
        print()
        
        if analysis_ok:
            print("🎉 All tests passed! Network analysis is working correctly.")
        else:
            print("❌ Network analysis tests failed.")
    else:
        print("❌ Configuration tests failed.")
