#!/usr/bin/env python3
"""
Clean test script for validating pipeline output generation.
This script tests the core pipeline functionality and generates outputs in test/outputs.
"""

import sys
import os
import time
import json
from datetime import datetime

# Add lib to path for module imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))

def test_pipeline_output_generation():
    """
    Test pipeline output generation with clean, organized results.
    Generates timestamped output directories with comprehensive analysis results.
    """
    from kbase_protein_query_module.kbase_protein_query_moduleImpl import kbase_protein_query_module
    
    # Initialize implementation
    impl = kbase_protein_query_module({})
    ctx = {'provenance': [{'ws_name': 'test_workspace'}]}
    
    # Test configuration
    test_cases = [
        {
            'name': 'single_protein',
            'proteins': ['MKVLLTLLCVVVLVLPGCELASTKPVKNAEVLHIQTAKNSVSGVQTPTDVVVSSDAGVTFYQYKPFVVDVTAKLVDCNKFGPDPVKAVDKTDVVFHVDVDVNVNVDVDVNVNVDVDVNVN'],
            'stages': ['sequence_analysis', 'family_assignment', 'similarity_search']
        },
        {
            'name': 'multi_protein',
            'proteins': [
                'MKVLLTLLCVVVLVLPGCELASTKPVKNAEVLHIQTAKNSVSGVQTPTDVVVSSDAGVTFYQYKPFVVDVTAKLVDCNKFGPDPVKAVDKTDVVFHVDVDVNVNVDVDVNVNVDVDVNVN',
                'MKLVLSLLLVGFLGAIILAVVVIVMVVVGVVVGVVVGVVVGVVVGVVVGVVVGVVVGVVVGVVVGVVVGVVVGVVVGVVVGVVVGVVVGVVVGVVVGVVVGVVVGVVVGVVVGVVVGVVV'
            ],
            'stages': ['sequence_analysis', 'family_assignment', 'similarity_search', 'multi_protein_analysis']
        }
    ]
    
    print("Testing Pipeline Output Generation")
    print("=" * 50)
    
    for test_case in test_cases:
        print(f"\nTesting: {test_case['name']}")
        print("-" * 30)
        
        # Create test parameters
        params = {
            'workspace_name': 'test_workspace',
            'input_proteins': test_case['proteins'],
            'analysis_stages': test_case['stages'],
            'output_format': 'all',
            'output_report_name': f"{test_case['name']}_analysis_report",
            'output_data_name': f"{test_case['name']}_analysis_data"
        }
        
        try:
            # Run analysis
            result = impl.run_protein_query_analysis(ctx, params)
            
            # Validate result
            if result and len(result) > 0:
                output = result[0]
                print(f"✅ Analysis completed successfully")
                print(f"   Report: {output.get('report_name', 'N/A')}")
                print(f"   Stages: {len(output.get('stages_completed', []))}")
                print(f"   Proteins: {output.get('protein_count', 0)}")
                
                # Check for output directory
                if 'output_directory' in output and output['output_directory']:
                    print(f"   Output Dir: {output['output_directory']}")
                    
                    # List files in output directory
                    if os.path.exists(output['output_directory']):
                        files = os.listdir(output['output_directory'])
                        print(f"   Generated Files: {len(files)}")
                        for file in sorted(files)[:5]:  # Show first 5 files
                            print(f"     - {file}")
                else:
                    print("   Warning: No output directory generated")
            else:
                print("❌ Analysis failed - no result returned")
                
        except Exception as e:
            print(f"❌ Analysis failed: {e}")
            import traceback
            traceback.print_exc()
    
    # Check test/outputs directory
    print(f"\n📁 Final Output Directory Check")
    print("-" * 30)
    outputs_dir = "test/outputs"
    if os.path.exists(outputs_dir):
        runs = [d for d in os.listdir(outputs_dir) if d.startswith('pipeline_run_')]
        print(f"Total pipeline runs: {len(runs)}")
        for run in sorted(runs)[-3:]:  # Show last 3 runs
            run_path = os.path.join(outputs_dir, run)
            files = os.listdir(run_path) if os.path.exists(run_path) else []
            print(f"  {run}: {len(files)} files")
    else:
        print("No outputs directory found")

if __name__ == "__main__":
    test_pipeline_output_generation()
