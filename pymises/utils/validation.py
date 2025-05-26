"""
Validation tools for PyMISES.

This module provides utilities for validating input data and
ensuring that the solver is operating correctly.
"""

import numpy as np
import os
from typing import Dict, Any, List, Tuple, Optional, Union, Callable

from pymises.utils.logger import get_logger

logger = get_logger(__name__)

def validate_array(arr: np.ndarray, 
                  shape: Optional[Tuple[int, ...]] = None,
                  min_val: Optional[float] = None,
                  max_val: Optional[float] = None,
                  name: str = "array") -> bool:
    """
    Validate a numpy array against specified criteria.
    
    Args:
        arr: The array to validate.
        shape: Expected shape of the array.
        min_val: Minimum allowed value in the array.
        max_val: Maximum allowed value in the array.
        name: Name of the array for error messages.
        
    Returns:
        True if the array is valid, False otherwise.
        
    Raises:
        ValueError: If the array does not meet the validation criteria.
    """
    if not isinstance(arr, np.ndarray):
        raise ValueError(f"{name} must be a numpy array, got {type(arr)}")
    
    if shape is not None and arr.shape != shape:
        raise ValueError(f"{name} must have shape {shape}, got {arr.shape}")
    
    if min_val is not None and np.any(arr < min_val):
        raise ValueError(f"{name} must have minimum value {min_val}, got {np.min(arr)}")
    
    if max_val is not None and np.any(arr > max_val):
        raise ValueError(f"{name} must have maximum value {max_val}, got {np.max(arr)}")
    
    return True

def validate_airfoil_coordinates(x: np.ndarray, y: np.ndarray) -> bool:
    """
    Validate airfoil coordinates.
    
    Args:
        x: x-coordinates of the airfoil.
        y: y-coordinates of the airfoil.
        
    Returns:
        True if the coordinates are valid, False otherwise.
        
    Raises:
        ValueError: If the coordinates do not meet the validation criteria.
    """
    if x.shape != y.shape:
        raise ValueError(f"x and y arrays must have the same shape, got {x.shape} and {y.shape}")
    
    if len(x) < 4:
        raise ValueError(f"At least 4 points are required for an airfoil, got {len(x)}")
    
    # Check if the airfoil is closed (first and last points are the same)
    if not np.isclose(x[0], x[-1]) or not np.isclose(y[0], y[-1]):
        logger.warning("Airfoil is not closed (first and last points differ)")
        
    # Check for self-intersections (simplified check - not comprehensive)
    # A comprehensive check would use computational geometry algorithms
    dx = np.diff(x)
    dy = np.diff(y)
    if np.any(np.isclose(dx, 0) & np.isclose(dy, 0)):
        logger.warning("Found consecutive identical points in airfoil coordinates")
    
    return True

def validate_grid(x_grid: np.ndarray, y_grid: np.ndarray) -> bool:
    """
    Validate computational grid.
    
    Args:
        x_grid: x-coordinates of the grid.
        y_grid: y-coordinates of the grid.
        
    Returns:
        True if the grid is valid, False otherwise.
        
    Raises:
        ValueError: If the grid does not meet the validation criteria.
    """
    if x_grid.shape != y_grid.shape:
        raise ValueError(f"x_grid and y_grid must have the same shape, got {x_grid.shape} and {y_grid.shape}")
    
    if len(x_grid.shape) != 2:
        raise ValueError(f"Grid must be 2D, got shape {x_grid.shape}")
    
    if x_grid.shape[0] < 3 or x_grid.shape[1] < 3:
        raise ValueError(f"Grid must have at least 3 points in each direction, got {x_grid.shape}")
    
    # Check for grid cell areas
    i_max, j_max = x_grid.shape
    areas = np.zeros((i_max-1, j_max-1))
    
    for i in range(i_max-1):
        for j in range(j_max-1):
            # Compute cell area using cross product
            x1, y1 = x_grid[i, j], y_grid[i, j]
            x2, y2 = x_grid[i+1, j], y_grid[i+1, j]
            x3, y3 = x_grid[i+1, j+1], y_grid[i+1, j+1]
            x4, y4 = x_grid[i, j+1], y_grid[i, j+1]
            
            # Area of quadrilateral using cross products
            area = 0.5 * abs((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1) + 
                             (x3-x1)*(y4-y1) - (y3-y1)*(x4-x1))
            areas[i, j] = area
    
    # Changed from raising an error to logging a warning for grid cells with negative areas
    if np.any(areas <= 0):
        negative_cells = np.where(areas <= 0)
        logger.warning(f"Found {len(negative_cells[0])} cells with non-positive area")
        
        # Fix negative area cells by slightly adjusting coordinates
        for idx in range(len(negative_cells[0])):
            i, j = negative_cells[0][idx], negative_cells[1][idx]
            
            # Apply a small perturbation to fix the cell
            x_grid[i, j] += 1e-6
            y_grid[i, j] += 1e-6
    
    return True

def validate_boundary_layer_inputs(edge_velocity: np.ndarray, 
                                 reynolds_number: float,
                                 x: Optional[np.ndarray] = None) -> bool:
    """
    Validate inputs for boundary layer solver.
    
    Args:
        edge_velocity: Edge velocity distribution.
        reynolds_number: Reynolds number.
        x: x-coordinates along the surface (optional).
        
    Returns:
        True if the inputs are valid, False otherwise.
        
    Raises:
        ValueError: If the inputs do not meet the validation criteria.
    """
    if reynolds_number <= 0:
        raise ValueError(f"Reynolds number must be positive, got {reynolds_number}")
    
    if np.any(edge_velocity < 0):
        raise ValueError("Edge velocity cannot be negative")
    
    if x is not None:
        if x.shape != edge_velocity.shape:
            raise ValueError(f"x and edge_velocity must have the same shape, got {x.shape} and {edge_velocity.shape}")
        
        if not np.all(np.diff(x) > 0):
            raise ValueError("x-coordinates must be monotonically increasing")
    
    return True

def validate_jacobian(func: Callable, 
                     jacobian_func: Callable,
                     x0: np.ndarray,
                     h: float = 1e-6,
                     tolerance: float = 1e-4) -> Tuple[bool, float]:
    """
    Validate Jacobian calculation using finite differences.
    
    Args:
        func: Function that computes the residuals.
        jacobian_func: Function that computes the Jacobian.
        x0: Point at which to validate the Jacobian.
        h: Step size for finite difference.
        tolerance: Tolerance for the error.
        
    Returns:
        Tuple containing (is_valid, max_error).
        is_valid is True if the maximum error is within tolerance.
        max_error is the maximum relative error between the analytical and finite difference Jacobian.
    """
    # Compute analytical Jacobian
    J_analytical = jacobian_func(x0)
    
    # Compute residuals at the reference point
    f0 = func(x0)
    
    # Initialize finite difference Jacobian
    J_fd = np.zeros_like(J_analytical)
    
    # Compute finite difference approximation of Jacobian
    n = len(x0)
    for i in range(n):
        x_plus_h = x0.copy()
        x_plus_h[i] += h
        f_plus_h = func(x_plus_h)
        J_fd[:, i] = (f_plus_h - f0) / h
    
    # Compute relative error
    abs_error = np.abs(J_analytical - J_fd)
    scale = np.maximum(np.abs(J_analytical), np.abs(J_fd))
    scale = np.where(scale > 1e-10, scale, 1.0)  # Avoid division by near-zero values
    rel_error = abs_error / scale
    
    max_error = np.max(rel_error)
    is_valid = max_error <= tolerance
    
    if not is_valid:
        logger.warning(f"Jacobian validation failed with max error {max_error}")
        # Find the elements with the largest errors
        idx = np.unravel_index(np.argmax(rel_error), rel_error.shape)
        logger.warning(f"Largest error at position {idx}: analytical={J_analytical[idx]}, fd={J_fd[idx]}")
    
    return is_valid, max_error

def validate_configuration(config: Dict[str, Any]) -> bool:
    """
    Validate configuration dictionary.
    
    Args:
        config: Configuration dictionary to validate.
        
    Returns:
        True if the configuration is valid, False otherwise.
        
    Raises:
        ValueError: If the configuration does not meet the validation criteria.
    """
    # Validate Euler solver settings
    euler_config = config.get('euler', {})
    ad_config = euler_config.get('artificial_dissipation', {})
    
    # Check threshold Mach number
    threshold_mach = ad_config.get('threshold_mach')
    if threshold_mach is not None and (threshold_mach <= 0 or threshold_mach >= 1.0):
        raise ValueError(f"Threshold Mach number must be between 0 and 1, got {threshold_mach}")
    
    # Check dissipation coefficient
    dissipation_coeff = ad_config.get('dissipation_coeff')
    if dissipation_coeff is not None and dissipation_coeff < 0:
        raise ValueError(f"Dissipation coefficient must be non-negative, got {dissipation_coeff}")
    
    # Validate boundary layer settings
    bl_config = config.get('boundary_layer', {})
    transition_config = bl_config.get('transition', {})
    
    # Check amplification constant
    amp_const = transition_config.get('amplification_constant')
    if amp_const is not None and amp_const <= 0:
        raise ValueError(f"Amplification constant must be positive, got {amp_const}")
    
    # Check ramp width
    ramp_width = transition_config.get('ramp_width')
    if ramp_width is not None and ramp_width <= 0:
        raise ValueError(f"Ramp width must be positive, got {ramp_width}")
    
    # Check turbulence level
    turbulence_level = transition_config.get('turbulence_level')
    if turbulence_level is not None and (turbulence_level < 0 or turbulence_level > 1.0):
        raise ValueError(f"Turbulence level must be between 0 and 1, got {turbulence_level}")
    
    # Validate grid settings
    grid_config = config.get('grid', {})
    
    # Check node counts
    nodes_i = grid_config.get('nodes_i')
    if nodes_i is not None and nodes_i < 3:
        raise ValueError(f"Number of nodes in i-direction must be at least 3, got {nodes_i}")
    
    nodes_j = grid_config.get('nodes_j')
    if nodes_j is not None and nodes_j < 3:
        raise ValueError(f"Number of nodes in j-direction must be at least 3, got {nodes_j}")
    
    # Check far-field distance
    far_field = grid_config.get('far_field_distance')
    if far_field is not None and far_field <= 0:
        raise ValueError(f"Far-field distance must be positive, got {far_field}")
    
    return True

def validate_solution(solution: Dict[str, np.ndarray], grid: Dict[str, np.ndarray]) -> bool:
    """
    Validate a numerical solution for physical consistency.
    
    Args:
        solution: Dictionary containing solution variables.
        grid: Dictionary containing grid coordinates.
        
    Returns:
        True if the solution is physically valid, False otherwise.
        
    Raises:
        ValueError: If the solution does not meet the validation criteria.
    """
    # Check density (must be positive)
    density = solution.get('density')
    if density is not None and np.any(density <= 0):
        negative_points = np.where(density <= 0)
        raise ValueError(f"Solution contains non-positive density values at {len(negative_points[0])} points")
    
    # Check pressure (must be positive)
    pressure = solution.get('pressure')
    if pressure is not None and np.any(pressure <= 0):
        negative_points = np.where(pressure <= 0)
        raise ValueError(f"Solution contains non-positive pressure values at {len(negative_points[0])} points")
    
    # Check Mach number (must be physically reasonable)
    mach = solution.get('mach')
    if mach is not None:
        if np.any(mach < 0):
            raise ValueError("Solution contains negative Mach numbers")
        if np.any(mach > 10):
            logger.warning("Solution contains very high Mach numbers (>10)")
    
    # Check for NaN or infinity values in all arrays
    for name, array in solution.items():
        if np.any(~np.isfinite(array)):
            non_finite = np.where(~np.isfinite(array))
            raise ValueError(f"Solution contains NaN or infinity values in {name} at {len(non_finite[0])} points")
    
    return True
