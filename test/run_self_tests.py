#!/usr/bin/env python3
"""
Self-test runner for all module components.

This script runs the main() functions from all Python files that have self-tests.
It's independent of KBase SDK tests and can be run directly.

Usage:
    python test/run_self_tests.py
    # Or from project root:
    python -m test.run_self_tests
"""

import sys
import os
import subprocess
from pathlib import Path
from typing import List, Tuple

# Add lib to path
project_root = Path(__file__).parent.parent
lib_path = project_root / "lib"
sys.path.insert(0, str(lib_path))

# Color output for better visibility
class Colors:
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    END = '\033[0m'
    BOLD = '\033[1m'

def print_header(text: str):
    """Print a formatted header."""
    print(f"\n{Colors.BOLD}{Colors.BLUE}{'='*60}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.BLUE}{text}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.BLUE}{'='*60}{Colors.END}\n")

def print_success(text: str):
    """Print success message."""
    print(f"{Colors.GREEN}✅ {text}{Colors.END}")

def print_error(text: str):
    """Print error message."""
    print(f"{Colors.RED}❌ {text}{Colors.END}")

def print_warning(text: str):
    """Print warning message."""
    print(f"{Colors.YELLOW}⚠️  {text}{Colors.END}")

def run_test(module_path: str, test_name: str) -> Tuple[bool, str]:
    """
    Run a test by executing the main() function in a module.
    
    Args:
        module_path: Python module path (e.g., "kbase_protein_query_module.src.util.embeddings.generator")
        test_name: Human-readable test name
        
    Returns:
        Tuple of (success: bool, output: str)
    """
    try:
        # Import and run the main function
        module_parts = module_path.split('.')
        module = __import__(module_path, fromlist=[module_parts[-1]])
        
        if hasattr(module, 'main'):
            print(f"Running {test_name}...")
            result = module.main()
            
            # main() should return 0 for success, non-zero for failure
            if result == 0:
                return True, "PASSED"
            else:
                return False, f"FAILED with exit code {result}"
        else:
            return False, "No main() function found"
    except ImportError as e:
        return False, f"Import error: {e}"
    except Exception as e:
        return False, f"Error: {e}"

def run_test_via_subprocess(script_path: str, test_name: str) -> Tuple[bool, str]:
    """
    Run a test by executing a Python script directly.
    
    Args:
        script_path: Path to Python script
        test_name: Human-readable test name
        
    Returns:
        Tuple of (success: bool, output: str)
    """
    try:
        script_full_path = project_root / script_path
        if not script_full_path.exists():
            return False, f"Script not found: {script_path}"
        
        print(f"Running {test_name}...")
        result = subprocess.run(
            [sys.executable, str(script_full_path)],
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        
        if result.returncode == 0:
            return True, "PASSED"
        else:
            error_msg = result.stderr.strip() or result.stdout.strip()
            return False, f"FAILED: {error_msg[:200]}"
    except subprocess.TimeoutExpired:
        return False, "TIMEOUT (exceeded 5 minutes)"
    except Exception as e:
        return False, f"Error: {e}"

def main():
    """Run all self-tests."""
    print_header("Self-Test Runner")
    
    # List of tests to run
    # Format: (test_name, module_path_or_script_path, use_subprocess)
    tests = [
        # Core components
        ("Analysis Config", "kbase_protein_query_module.src.analysis.config", False),
        ("Analysis Manager", "kbase_protein_query_module.src.analysis.analysis_manager", False),
        ("Input Manager", "kbase_protein_query_module.src.input.input_manager", False),
        ("Output Manager", "kbase_protein_query_module.src.output.output_manager", False),
        ("Workflow Orchestrator", "kbase_protein_query_module.src.core.workflow_orchestrator", False),
        
        # Utility components
        ("Embedding Generator", "kbase_protein_query_module.src.util.embeddings.generator", False),
        ("Protein Storage", "kbase_protein_query_module.src.util.storage.storage", False),
        ("Similarity Search", "kbase_protein_query_module.src.util.storage.similarity_search", False),
        ("UniProt API", "kbase_protein_query_module.src.util.uniprot.api", False),
        
        # Analysis components (may require data files)
        ("Network Analysis", "kbase_protein_query_module.src.analysis.network_analysis.network_analysis", False),
    ]
    
    results = []
    passed = 0
    failed = 0
    skipped = 0
    
    for test_name, test_path, use_subprocess in tests:
        try:
            if use_subprocess:
                success, message = run_test_via_subprocess(test_path, test_name)
            else:
                success, message = run_test(test_path, test_name)
            
            if success:
                print_success(f"{test_name}: {message}")
                passed += 1
            else:
                print_error(f"{test_name}: {message}")
                failed += 1
            
            results.append((test_name, success, message))
        except KeyboardInterrupt:
            print_warning(f"Test interrupted by user")
            break
        except Exception as e:
            print_error(f"{test_name}: Unexpected error - {e}")
            failed += 1
            results.append((test_name, False, f"Unexpected error: {e}"))
    
    # Print summary
    print_header("Test Summary")
    print(f"Total tests: {len(results)}")
    print_success(f"Passed: {passed}")
    if failed > 0:
        print_error(f"Failed: {failed}")
    if skipped > 0:
        print_warning(f"Skipped: {skipped}")
    
    # Print failed tests
    if failed > 0:
        print_header("Failed Tests")
        for test_name, success, message in results:
            if not success:
                print_error(f"{test_name}: {message}")
    
    # Return exit code
    return 0 if failed == 0 else 1

if __name__ == "__main__":
    sys.exit(main())

