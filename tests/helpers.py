"""
Helper functions and utilities for PyMISES tests.

This module provides common utilities and helper functions used across
different test modules to reduce code duplication and simplify tests.
"""

import numpy as np
import os
from typing import Tuple, Dict, List, Optional, Union

def create_naca_airfoil(naca_code: str = '0012', n_points: int = 101) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create coordinates for a NACA 4-digit airfoil.
    
    Args:
        naca_code: The 4-digit NACA code (e.g., '0012').
        n_points: Number of points to generate.
        
    Returns:
        Tuple containing arrays of x and y coordinates.
    """
    # Parse NACA digits
    try:
        if len(naca_code) != 4:
            raise ValueError(f"NACA code must be 4 digits, got {naca_code}")
        
        m = int(naca_code[0]) / 100.0
        p = int(naca_code[1]) / 10.0
        t = int(naca_code[2:]) / 100.0
    except ValueError:
        raise ValueError(f"Invalid NACA code: {naca_code}")
    
    # Create x coordinates with clustering near leading and trailing edges
    beta = np.linspace(0, np.pi, n_points)
    x = 0.5 * (1 - np.cos(beta))  # Cosine clustering
    
    # Calculate thickness distribution
    y_t = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)
    
    # Calculate camber line and slope
    y_c = np.zeros_like(x)
    dyc_dx = np.zeros_like(x)
    
    if m > 0 and p > 0:  # Cambered airfoil
        # Front part of the camber line (x < p)
        mask_front = x <= p
        y_c[mask_front] = m * (x[mask_front] / p**2) * (2 * p - x[mask_front])
        dyc_dx[mask_front] = m * (2 * p - 2 * x[mask_front]) / p**2
        
        # Rear part of the camber line (x ≥ p)
        mask_rear = x > p
        y_c[mask_rear] = m * ((1 - x[mask_rear]) / (1 - p)**2) * (1 + x[mask_rear] - 2 * p)
        dyc_dx[mask_rear] = m * (2 * p - 2 * x[mask_rear]) / (1 - p)**2
    
    # Calculate the coordinates of the upper and lower surface
    theta = np.arctan(dyc_dx)
    
    x_upper = x - y_t * np.sin(theta)
    y_upper = y_c + y_t * np.cos(theta)
    
    x_lower = x + y_t * np.sin(theta)
    y_lower = y_c - y_t * np.cos(theta)
    
    # Combine the upper and lower surface points
    # Start from the trailing edge, go around the leading edge, and back to the trailing edge
    x_airfoil = np.concatenate([np.flipud(x_upper), x_lower[1:]])
    y_airfoil = np.concatenate([np.flipud(y_upper), y_lower[1:]])
    
    return x_airfoil, y_airfoil

def create_flat_plate(length: float = 1.0, n_points: int = 101) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create coordinates for a flat plate (useful for validation tests).
    
    Args:
        length: Length of the flat plate.
        n_points: Number of points to generate.
        
    Returns:
        Tuple containing arrays of x and y coordinates.
    """
    # Create x coordinates with clustering near leading and trailing edges
    beta = np.linspace(0, np.pi, n_points)
    x = length * 0.5 * (1 - np.cos(beta))  # Cosine clustering
    
    # Flat plate has y=0
    y = np.zeros_like(x)
    
    return x, y

def create_simple_grid(ni: int = 11, nj: int = 11, 
                      x_min: float = 0.0, x_max: float = 1.0,
                      y_min: float = 0.0, y_max: float = 1.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create a simple rectangular grid (useful for simple test cases).
    
    Args:
        ni: Number of points in i-direction.
        nj: Number of points in j-direction.
        x_min, x_max: x-coordinate range.
        y_min, y_max: y-coordinate range.
        
    Returns:
        Tuple containing 2D arrays of x and y coordinates with shape (ni, nj).
    """
    # Create 1D coordinate arrays
    x_1d = np.linspace(x_min, x_max, ni)
    y_1d = np.linspace(y_min, y_max, nj)
    
    # Create 2D grid
    x, y = np.meshgrid(x_1d, y_1d, indexing='ij')
    
    return x, y

def create_circular_grid(ni: int = 36, nj: int = 10, 
                        r_inner: float = 1.0, r_outer: float = 5.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create a circular/O-grid (useful for testing O-grid functionality).
    
    Args:
        ni: Number of points around the circle.
        nj: Number of points in radial direction.
        r_inner: Inner radius.
        r_outer: Outer radius.
        
    Returns:
        Tuple containing 2D arrays of x and y coordinates with shape (ni, nj).
    """
    # Create 1D arrays for angle and radius
    theta = np.linspace(0, 2*np.pi, ni, endpoint=False)
    r = np.linspace(r_inner, r_outer, nj)
    
    # Initialize 2D arrays
    x = np.zeros((ni, nj))
    y = np.zeros((ni, nj))
    
    # Fill the grid
    for i in range(ni):
        for j in range(nj):
            x[i, j] = r[j] * np.cos(theta[i])
            y[i, j] = r[j] * np.sin(theta[i])
    
    return x, y

def get_test_data_path(filename: str) -> str:
    """
    Get the full path to a test data file.
    
    Args:
        filename: Name of the test data file.
        
    Returns:
        Full path to the test data file.
    """
    test_data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(test_data_dir, filename)

def verify_solution_consistency(solution: Dict[str, np.ndarray]) -> bool:
    """
    Verify that a solution dictionary contains consistent and physically valid data.
    
    Args:
        solution: Solution dictionary with flow variables.
        
    Returns:
        True if the solution is consistent and physically valid.
        
    Raises:
        ValueError: If the solution is inconsistent or physically invalid.
    """
    # Check that required fields are present
    required_fields = ['density', 'pressure', 'velocity_x', 'velocity_y']
    for field in required_fields:
        if field not in solution:
            raise ValueError(f"Solution is missing required field: {field}")
    
    # Check for positive density and pressure
    if np.any(solution['density'] <= 0):
        raise ValueError("Solution contains non-positive density values")
    
    if np.any(solution['pressure'] <= 0):
        raise ValueError("Solution contains non-positive pressure values")
    
    # Check for finite values
    for field, values in solution.items():
        if np.any(~np.isfinite(values)):
            raise ValueError(f"Solution contains non-finite values in {field}")
    
    return True

def calculate_error_norms(numerical: np.ndarray, exact: np.ndarray) -> Dict[str, float]:
    """
    Calculate error norms between numerical and exact solutions.
    
    Args:
        numerical: Numerical solution array.
        exact: Exact solution array.
        
    Returns:
        Dictionary containing L1, L2, and Linf error norms.
    """
    # Compute absolute error
    error = numerical - exact
    
    # Calculate error norms
    l1_norm = np.mean(np.abs(error))
    l2_norm = np.sqrt(np.mean(error**2))
    linf_norm = np.max(np.abs(error))
    
    return {
        'L1': l1_norm,
        'L2': l2_norm,
        'Linf': linf_norm
    }
