"""
PyMISES test package.

This package contains tests for the PyMISES implementation of the Multiple-blade
Interacting Streamtube Euler Solver originally developed by Mark Drela at MIT.
"""

# Make sure the pymises package is importable from tests
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
