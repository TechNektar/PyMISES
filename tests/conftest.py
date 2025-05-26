"""
PyTest configuration and fixtures for PyMISES tests.

This module provides pytest fixtures that can be used across different test files
to avoid code duplication and ensure consistent test setup.
"""

import pytest
import numpy as np
import os
import sys

from pymises.core.geometry import AirfoilGeometry, BladeGeometry, CascadeGeometry
from pymises.core.grid import StreamlineGrid, GridGenerator
from pymises.core.euler import EulerSolver
from pymises.core.boundary_layer import BoundaryLayerSolver, BoundaryLayerFactory
from pymises.core.newton import NewtonSolver

from tests.helpers import create_naca_airfoil, create_simple_grid, create_circular_grid

# Set to True to save test artifacts (grids, solutions, etc.) to disk
SAVE_TEST_ARTIFACTS = False

@pytest.fixture
def naca0012():
    """Fixture providing a NACA 0012 airfoil geometry."""
    x, y = create_naca_airfoil('0012', n_points=101)
    return AirfoilGeometry(x, y)

@pytest.fixture
def naca4412():
    """Fixture providing a NACA 4412 airfoil geometry."""
    x, y = create_naca_airfoil('4412', n_points=101)
    return AirfoilGeometry(x, y)

@pytest.fixture
def simple_cascade(naca0012):
    """Fixture providing a simple cascade geometry."""
    return CascadeGeometry(
        blade=naca0012,
        stagger_angle=30.0,
        pitch=1.0,
        chord=1.0
    )

@pytest.fixture
def simple_grid():
    """Fixture providing a simple rectangular grid."""
    x, y = create_simple_grid(ni=11, nj=11)
    return StreamlineGrid(x, y)

@pytest.fixture
def circular_grid():
    """Fixture providing a circular/O-grid."""
    x, y = create_circular_grid(ni=36, nj=10)
    return StreamlineGrid(x, y)

@pytest.fixture
def airfoil_o_grid(naca0012):
    """Fixture providing an O-grid around a NACA 0012 airfoil."""
    grid_gen = GridGenerator(naca0012, {
        'ni': 31,
        'nj': 15,
        'far_field_distance': 10.0
    })
    return grid_gen.generate_grid(grid_type='o-grid')

@pytest.fixture
def airfoil_c_grid(naca0012):
    """Fixture providing a C-grid around a NACA 0012 airfoil."""
    grid_gen = GridGenerator(naca0012, {
        'ni': 31,
        'nj': 15,
        'far_field_distance': 10.0,
        'wake_length': 5.0
    })
    return grid_gen.generate_grid(grid_type='c-grid')

@pytest.fixture
def cascade_grid(simple_cascade):
    """Fixture providing a grid for a cascade."""
    grid_gen = GridGenerator(simple_cascade, {
        'ni': 31,
        'nj': 15
    })
    return grid_gen.generate_grid(grid_type='cascade')

@pytest.fixture
def euler_solver(airfoil_o_grid):
    """Fixture providing an Euler solver with a basic grid."""
    return EulerSolver(airfoil_o_grid)

@pytest.fixture
def boundary_layer_solver():
    """Fixture providing a simple boundary layer solver for a flat plate."""
    n_points = 101
    x = np.linspace(0.001, 1.0, n_points)
    edge_velocity = np.ones_like(x)
    reynolds = 1.0e6
    
    return BoundaryLayerSolver(
        x=x,
        edge_velocity=edge_velocity,
        reynolds_number=reynolds
    )

@pytest.fixture
def newton_solver():
    """Fixture providing a Newton solver for simple test problems."""
    # Define a simple nonlinear system: F(x) = 0
    # F1(x) = x1^2 + x2^2 - 4 = 0 (circle with radius 2)
    # F2(x) = x1 - x2 = 0 (line through origin)
    
    def residual_function(x):
        F = np.zeros(2)
        F[0] = x[0]**2 + x[1]**2 - 4
        F[1] = x[0] - x[1]
        return F
    
    def jacobian_function(x):
        J = np.zeros((2, 2))
        J[0, 0] = 2 * x[0]
        J[0, 1] = 2 * x[1]
        J[1, 0] = 1
        J[1, 1] = -1
        return J
    
    # Initial guess
    x0 = np.array([1.0, 0.0])
    
    return NewtonSolver(
        residual_function=residual_function,
        jacobian_function=jacobian_function,
        solution=x0
    )

@pytest.fixture
def save_test_results(request):
    """Fixture for saving test results (grids, solutions, etc.) to disk."""
    if not SAVE_TEST_ARTIFACTS:
        return lambda name, data: None
    
    # Get the test function name
    test_name = request.function.__name__
    
    # Create output directory if it doesn't exist
    output_dir = os.path.join(os.path.dirname(__file__), 'output')
    os.makedirs(output_dir, exist_ok=True)
    
    def _save_data(name, data):
        """Save data to disk."""
        # Create filename
        filename = f"{test_name}_{name}"
        
        # Save based on data type
        if isinstance(data, np.ndarray):
            np.save(os.path.join(output_dir, filename), data)
        elif isinstance(data, dict):
            np.savez(os.path.join(output_dir, filename), **data)
        else:
            # Just save as text
            with open(os.path.join(output_dir, f"{filename}.txt"), 'w') as f:
                f.write(str(data))
    
    return _save_data
