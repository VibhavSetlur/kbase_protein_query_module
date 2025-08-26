#!/usr/bin/env python3
"""
Clean and Run Script for KBase Protein Query Module Tests

This script cleans up old output directories and runs the unified test runner.
"""

import os
import shutil
from pathlib import Path

def cleanup_old_outputs():
    """Clean up old output directories."""
    print("🧹 Cleaning up old output directories...")
    
    # Directories to clean up
    cleanup_dirs = [
        'exports',
        'html_reports', 
        'test_outputs',
        'reports',
        'output',
        'test_outputs'
    ]
    
    for dir_name in cleanup_dirs:
        dir_path = Path(dir_name)
        if dir_path.exists() and dir_path.is_dir():
            try:
                shutil.rmtree(dir_path)
                print(f"✅ Cleaned up: {dir_path}")
            except Exception as e:
                print(f"⚠️  Could not clean up {dir_path}: {e}")
    
    print("✅ Cleanup completed!")

def main():
    """Main entry point."""
    print("🧪 KBase Protein Query Module - Clean and Run")
    print("=" * 50)
    
    # Clean up old outputs
    cleanup_old_outputs()
    
    print("\n🚀 Running unified test suite...")
    print("=" * 50)
    
    # Run the unified test runner
    import run_tests
    run_tests.main()

if __name__ == "__main__":
    main()
