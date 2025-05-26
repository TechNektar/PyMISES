#!/usr/bin/env python3
"""
PyMISES Test Runner

This script runs the PyMISES test suite. It can run all tests or specific test modules.

Usage:
    python run_tests.py             # Run all tests
    python run_tests.py core        # Run tests in the core module
    python run_tests.py grid euler  # Run grid and euler tests
    python run_tests.py -l          # List available test modules
    python run_tests.py -v          # Run tests in verbose mode
    python run_tests.py -i          # Run integration tests only
    python run_tests.py -u          # Run unit tests only
    python run_tests.py -h          # Show help

Examples:
    python run_tests.py                       # Run all tests
    python run_tests.py core                  # Run all core module tests
    python run_tests.py boundary_conditions   # Run all boundary condition tests
    python run_tests.py grid newton           # Run grid and newton tests
    python run_tests.py -i                   # Run integration tests only
    python run_tests.py -v grid              # Run grid tests with verbose output
"""

import sys
import os
import unittest
import argparse
import importlib
import fnmatch
import time
import traceback
import logging
import datetime

# Configure logging
def setup_logger():
    """Set up a logger for error reporting."""
    # Generate timestamp for log file name
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = f"test_run_{timestamp}.log"
    
    # Create or ensure logs directory exists
    logs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "logs")
    os.makedirs(logs_dir, exist_ok=True)
    log_path = os.path.join(logs_dir, log_file)
    
    # Configure the logger
    logger = logging.getLogger("pymises_tests")
    logger.setLevel(logging.DEBUG)
    
    # Create file handler
    file_handler = logging.FileHandler(log_path, mode="w")
    file_handler.setLevel(logging.DEBUG)
    
    # Create formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    file_handler.setFormatter(formatter)
    
    # Add handlers to logger
    logger.addHandler(file_handler)
    
    # Write header
    logger.info("PyMISES Test Runner Log")
    logger.info("=" * 50)
    logger.info(f"Test run started at: {datetime.datetime.now()}")
    logger.info(f"Python version: {sys.version}")
    logger.info(f"Working directory: {os.getcwd()}")
    logger.info("=" * 50)
    
    print(f"\nLog file: {log_path}\n")
    
    # Return logger and log path
    return logger, log_path

def run_unit_tests(module_names=None, verbosity=1, logger=None):
    """Run unit tests for specific modules or all if none specified."""
    # Define the test directory path
    project_root = os.path.dirname(os.path.abspath(__file__))
    test_dir = os.path.join(project_root, 'tests')
    
    # Ensure the current working directory is the project root
    # This helps unittest discovery find modules correctly
    os.chdir(project_root)
    
    # Get all test directories (exclude special directories)
    def is_test_dir(d):
        if not os.path.isdir(os.path.join(test_dir, d)):
            return False
        if d.startswith('__'):  # Skip __pycache__, __init__, etc.
            return False
        if d == 'integration':   # Integration tests are handled separately
            return False
        return True
    
    test_dirs = [d for d in os.listdir(test_dir) if is_test_dir(d)]
    
    # If specific modules are requested, filter to only those
    filtered_dirs = test_dirs
    if module_names:
        filtered_dirs = []
        for name in module_names:
            # Check for exact directory match
            if name in test_dirs:
                filtered_dirs.append(name)
            else:
                # Check for module name match (e.g., 'grid' matching 'test_grid.py')
                for d in test_dirs:
                    module_files = [f[5:-3] for f in os.listdir(os.path.join(test_dir, d)) 
                                   if f.startswith('test_') and f.endswith('.py')]
                    if name in module_files:
                        filtered_dirs.append(d)
                        break
    
    # Create test suite
    test_suite = unittest.TestSuite()
    
    # Discover and add tests from each directory
    for dir_name in filtered_dirs:
        dir_path = os.path.join(test_dir, dir_name)
        
        # If specific modules are requested, filter tests to only those
        pattern = 'test_*.py'
        if module_names:
            patterns = []
            for name in module_names:
                patterns.append(f'test_{name}.py')
            
            # If no specific test matches, use all tests from the directory
            if not patterns:
                patterns = [pattern]
            
            # Create a test suite for each pattern
            for p in patterns:
                try:
                    # Check if the directory exists before discovering tests
                    if not os.path.isdir(dir_path):
                        print(f"Warning: Directory {dir_path} not found, skipping")
                        continue
                        
                    # Check if the directory contains any test files matching the pattern
                    test_files = [f for f in os.listdir(dir_path) if f.startswith('test_') and f.endswith('.py')]
                    if not any(p == pattern or fnmatch.fnmatch(f, p) for f in test_files for p in patterns):
                        print(f"No test files matching {patterns} found in {dir_path}, skipping")
                        continue
                    
                    suite = unittest.defaultTestLoader.discover(dir_path, pattern=p, top_level_dir=project_root)
                    test_suite.addTest(suite)
                except Exception as e:
                    error_msg = f"Error discovering tests in {dir_path} with pattern {p}: {e}"
                    print(error_msg)
                    if logger:
                        logger.error(error_msg)
                        logger.error(f"Current working directory: {os.getcwd()}")
                        logger.error(f"Project root: {project_root}")
                        logger.error("Traceback: " + ''.join(traceback.format_exc()))
                    else:
                        print(f"Current working directory: {os.getcwd()}")
                        print(f"Project root: {project_root}")
                        print("Traceback:")
                        traceback.print_exc()
                    continue  # Continue with other patterns instead of raising
        else:
            # Load all tests from the directory
            try:
                # Check if the directory exists before discovering tests
                if not os.path.isdir(dir_path):
                    print(f"Warning: Directory {dir_path} not found, skipping")
                    continue
                    
                # Check if the directory contains any test files matching the pattern
                test_files = [f for f in os.listdir(dir_path) if f.startswith('test_') and f.endswith('.py')]
                if not any(fnmatch.fnmatch(f, pattern) for f in test_files):
                    print(f"No test files matching {pattern} found in {dir_path}, skipping")
                    continue
                
                suite = unittest.defaultTestLoader.discover(dir_path, pattern=pattern, top_level_dir=project_root)
                test_suite.addTest(suite)
            except Exception as e:
                error_msg = f"Error discovering tests in {dir_path} with pattern {pattern}: {e}"
                print(error_msg)
                if logger:
                    logger.error(error_msg)
                    logger.error(f"Current working directory: {os.getcwd()}")
                    logger.error(f"Project root: {project_root}")
                    logger.error("Traceback: " + ''.join(traceback.format_exc()))
                else:
                    print(f"Current working directory: {os.getcwd()}")
                    print(f"Project root: {project_root}")
                    print("Traceback:")
                    traceback.print_exc()
                continue  # Continue with other directories
    
    # Run tests
    test_runner = unittest.TextTestRunner(verbosity=verbosity)
    return test_runner.run(test_suite)

def run_integration_tests(module_names=None, verbosity=1, logger=None):
    """Run integration tests for specific modules or all if none specified."""
    # Define the integration test directory path
    project_root = os.path.dirname(os.path.abspath(__file__))
    integration_dir = os.path.join(project_root, 'tests', 'integration')
    
    # Ensure the current working directory is the project root
    # This helps unittest discovery find modules correctly
    os.chdir(project_root)
    
    # Check if the integration directory exists
    if not os.path.isdir(integration_dir):
        print(f"Warning: Integration test directory '{integration_dir}' not found")
        return unittest.TestResult()
    
    # Create test suite
    test_suite = unittest.TestSuite()
    
    # Discover and add tests from the integration directory
    if module_names:
        for name in module_names:
            pattern = f'test_*{name}*.py'
            try:
                # Check if there are any test files matching the pattern
                test_files = [f for f in os.listdir(integration_dir) if f.startswith('test_') and f.endswith('.py')]
                if not any(fnmatch.fnmatch(f, pattern) for f in test_files):
                    print(f"No integration test files matching {pattern} found, skipping")
                    continue
                    
                suite = unittest.defaultTestLoader.discover(integration_dir, pattern=pattern, top_level_dir=project_root)
                test_suite.addTest(suite)
            except Exception as e:
                error_msg = f"Error discovering tests in {integration_dir} with pattern {pattern}: {e}"
                print(error_msg)
                if logger:
                    logger.error(error_msg)
                    logger.error(f"Current working directory: {os.getcwd()}")
                    logger.error(f"Project root: {project_root}")
                    logger.error("Traceback: " + ''.join(traceback.format_exc()))
                else:
                    print(f"Current working directory: {os.getcwd()}")
                    print(f"Project root: {project_root}")
                    print("Traceback:")
                    traceback.print_exc()
                continue  # Continue with other modules instead of raising
    else:
        # Load all integration tests
        pattern = 'test_*.py'
        try:
            # Check if there are any test files matching the pattern
            test_files = [f for f in os.listdir(integration_dir) if f.startswith('test_') and f.endswith('.py')]
            if not test_files:
                print("No integration test files found")
            else:
                suite = unittest.defaultTestLoader.discover(integration_dir, pattern=pattern, top_level_dir=project_root)
                test_suite.addTest(suite)
        except Exception as e:
            error_msg = f"Error discovering integration tests: {e}"
            print(error_msg)
            if logger:
                logger.error(error_msg)
                logger.error(f"Current working directory: {os.getcwd()}")
                logger.error(f"Project root: {project_root}")
                logger.error("Traceback: " + ''.join(traceback.format_exc()))
            else:
                print(f"Current working directory: {os.getcwd()}")
                print(f"Project root: {project_root}")
                print("Traceback:")
                traceback.print_exc()
    
    # Run tests
    test_runner = unittest.TextTestRunner(verbosity=verbosity)
    return test_runner.run(test_suite)

def list_available_modules():
    """List available test modules."""
    print("Available test modules:")
    
    # Get test directory
    test_dir = os.path.join(os.path.dirname(__file__), 'tests')
    
    # Get all test directories
    test_dirs = [d for d in os.listdir(test_dir) 
                 if os.path.isdir(os.path.join(test_dir, d))]
    
    for dir_name in sorted(test_dirs):
        print(f"\n{dir_name.upper()}:")
        dir_path = os.path.join(test_dir, dir_name)
        
        # Get all test files in the directory
        test_files = [f for f in os.listdir(dir_path) 
                     if f.startswith('test_') and f.endswith('.py')]
        
        # Extract module names
        module_names = [f[5:-3] for f in test_files]
        
        for module in sorted(module_names):
            print(f"  - {module}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PyMISES Test Runner')
    parser.add_argument('modules', nargs='*', help='Specific modules to test')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    parser.add_argument('-l', '--list', action='store_true', help='List available test modules')
    parser.add_argument('-i', '--integration', action='store_true', help='Run integration tests only')
    parser.add_argument('-u', '--unit', action='store_true', help='Run unit tests only')
    
    args = parser.parse_args()
    
    if args.list:
        list_available_modules()
        sys.exit(0)
    
    verbosity = 2 if args.verbose else 1
    start_time = time.time()
    
    # Set up logger
    logger, log_path = setup_logger()
    logger.info(f"Command-line arguments: {args}")
    
    # Determine which tests to run
    run_units = not args.integration or args.unit
    run_integrations = not args.unit or args.integration
    
    # Run the tests
    unit_result = None
    integration_result = None
    
    if run_units:
        logger.info("Starting unit tests")
        print("\n========== RUNNING UNIT TESTS ==========\n")
        unit_result = run_unit_tests(args.modules, verbosity, logger)
        logger.info(f"Unit tests completed: {unit_result.testsRun} run, "
                  f"{len(unit_result.failures)} failures, "
                  f"{len(unit_result.errors)} errors")
    
    if run_integrations:
        logger.info("Starting integration tests")
        print("\n========== RUNNING INTEGRATION TESTS ==========\n")
        integration_result = run_integration_tests(args.modules, verbosity, logger)
        logger.info(f"Integration tests completed: {integration_result.testsRun} run, "
                  f"{len(integration_result.failures)} failures, "
                  f"{len(integration_result.errors)} errors")
    
    # Calculate total time
    total_time = time.time() - start_time
    
    # Print summary
    print("\n========== TEST SUMMARY ==========\n")
    logger.info("========== TEST SUMMARY ==========")
    
    unit_run = 0
    unit_failures = 0
    unit_errors = 0
    
    integration_run = 0
    integration_failures = 0
    integration_errors = 0
    
    if unit_result:
        unit_run = unit_result.testsRun
        unit_failures = len(unit_result.failures)
        unit_errors = len(unit_result.errors)
        
        print(f"Unit Tests: {unit_run} tests, {unit_failures} failures, {unit_errors} errors")
        
        # Log unit test failures
        if unit_failures > 0:
            logger.error("Unit Test Failures:")
            for i, (test, traceback) in enumerate(unit_result.failures):
                logger.error(f"Failure {i+1}: {test}")
                logger.error(f"Traceback: {traceback}\n")
        
        # Log unit test errors
        if unit_errors > 0:
            logger.error("Unit Test Errors:")
            for i, (test, traceback) in enumerate(unit_result.errors):
                logger.error(f"Error {i+1}: {test}")
                logger.error(f"Traceback: {traceback}\n")
    
    if integration_result:
        integration_run = integration_result.testsRun
        integration_failures = len(integration_result.failures)
        integration_errors = len(integration_result.errors)
        
        print(f"Integration Tests: {integration_run} tests, {integration_failures} failures, {integration_errors} errors")
        
        # Log integration test failures
        if integration_failures > 0:
            logger.error("Integration Test Failures:")
            for i, (test, traceback) in enumerate(integration_result.failures):
                logger.error(f"Failure {i+1}: {test}")
                logger.error(f"Traceback: {traceback}\n")
        
        # Log integration test errors
        if integration_errors > 0:
            logger.error("Integration Test Errors:")
            for i, (test, traceback) in enumerate(integration_result.errors):
                logger.error(f"Error {i+1}: {test}")
                logger.error(f"Traceback: {traceback}\n")
    
    total_run = unit_run + integration_run
    total_failures = unit_failures + integration_failures
    total_errors = unit_errors + integration_errors
    
    summary_msg = f"\nTotal: {total_run} tests, {total_failures} failures, {total_errors} errors"
    time_msg = f"Time: {total_time:.2f} seconds"
    
    print(summary_msg)
    print(time_msg)
    
    logger.info(summary_msg)
    logger.info(time_msg)
    
    # Final log message
    result_status = "FAILED" if total_failures > 0 or total_errors > 0 else "PASSED"
    logger.info(f"Test run {result_status}")
    logger.info("End of test run")
    
    # Print reminder about log file if there were errors
    if total_errors > 0 or total_failures > 0:
        print(f"\nDetailed error information saved to: {log_path}")
    
    # Set exit code based on test results
    if total_failures > 0 or total_errors > 0:
        sys.exit(1)
    else:
        sys.exit(0)
