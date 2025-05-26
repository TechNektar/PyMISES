"""
Jacobian calculation for the Euler solver.

This module provides functions for computing the Jacobian matrix of the Euler equations
using finite differences or analytical methods.
"""

import numpy as np
import scipy.sparse as sp
import logging
from typing import Dict, Any, Callable, Optional, Tuple, List, Union

logger = logging.getLogger(__name__)

def compute_euler_jacobian_fd(residual_function: Callable, solution_vector: np.ndarray,
                            step_size: float = 1e-6) -> sp.csr_matrix:
    """
    Compute the Jacobian matrix of the Euler equations using finite differences.

    Args:
        residual_function: Function that computes the residual vector.
        solution_vector: Current solution vector.
        step_size: Step size for finite differences.

    Returns:
        Sparse Jacobian matrix.
    """
    n = len(solution_vector)

    # Compute base residual
    base_residual = residual_function(solution_vector)
    m = len(base_residual)

    # Initialize sparse matrix storage
    rows = []
    cols = []
    data = []

    # Estimate sparsity pattern based on typical Euler discretization
    # For each point, we expect dependencies on the point itself and its neighbors
    for j in range(n):
        # Perturb solution vector
        perturbed_solution = solution_vector.copy()
        perturbed_solution[j] += step_size

        # Compute perturbed residual
        try:
            perturbed_residual = residual_function(perturbed_solution)

            # Calculate derivatives
            derivatives = (perturbed_residual - base_residual) / step_size

            # Store non-zero derivatives
            for i in range(m):
                if abs(derivatives[i]) > 1e-14:  # Filter out very small values
                    rows.append(i)
                    cols.append(j)
                    data.append(derivatives[i])
        except Exception as e:
            logger.warning(f"Error computing Jacobian column {j}: {str(e)}")
            # Add diagonal entry to maintain matrix structure
            rows.append(j % m)
            cols.append(j)
            data.append(1.0)

    # Create sparse matrix
    jacobian = sp.csr_matrix((data, (rows, cols)), shape=(m, n))

    return jacobian

def compute_euler_jacobian_block(euler_solver, solution_vector: np.ndarray) -> sp.csr_matrix:
    """
    Compute the Jacobian matrix of the Euler equations using a block structure.

    This function computes the Jacobian matrix using a block structure based on
    the physical structure of the Euler equations. This is more efficient than
    a full finite difference approach.

    Args:
        euler_solver: EulerSolver instance.
        solution_vector: Current solution vector.

    Returns:
        Sparse Jacobian matrix.
    """
    # Get grid dimensions
    ni, nj = euler_solver.grid.ni, euler_solver.grid.nj
    n_points = ni * nj
    n_equations = euler_solver.n_equations
    n_vars = n_equations * n_points

    # Initialize sparse matrix
    jacobian = sp.lil_matrix((n_vars, n_vars))

    # Update solution from vector
    euler_solver.set_solution_from_vector(solution_vector)

    # Get solution variables
    density = euler_solver.solution['density'].flatten()
    velocity_x = euler_solver.solution['velocity_x'].flatten()
    velocity_y = euler_solver.solution['velocity_y'].flatten()
    pressure = euler_solver.solution['pressure'].flatten()

    # Calculate grid metrics (approximate cell sizes)
    x = euler_solver.grid.x
    y = euler_solver.grid.y

    # Compute approximate cell sizes
    dx = np.zeros_like(x)
    dy = np.zeros_like(y)

    # For interior points, use central differences
    dx[:, 1:-1] = 0.5 * (x[:, 2:] - x[:, :-2])
    dy[:, 1:-1] = 0.5 * (y[:, 2:] - y[:, :-2])

    # For boundary points, use one-sided differences
    dx[:, 0] = x[:, 1] - x[:, 0]
    dx[:, -1] = x[:, -1] - x[:, -2]
    dy[:, 0] = y[:, 1] - y[:, 0]
    dy[:, -1] = y[:, -1] - y[:, -2]

    # Ensure positive values and avoid division by zero
    dx = np.abs(dx) + 1e-10
    dy = np.abs(dy) + 1e-10

    # Flatten for use in calculations
    dx = dx.flatten()
    dy = dy.flatten()

    # Specific heat ratio
    gamma = euler_solver.thermo.gamma

    # Compute sound speed
    sound_speed = np.sqrt(gamma * pressure / density)

    # Compute local Mach number
    velocity_magnitude = np.sqrt(velocity_x**2 + velocity_y**2)
    mach_number = velocity_magnitude / sound_speed

    # Compute local time step (CFL condition)
    dt = euler_solver.config.get('cfl', 0.8) * np.minimum(dx, dy) / (sound_speed + velocity_magnitude)

    # Compute flux Jacobians for each point
    for i in range(ni):
        for j in range(nj):
            idx = j * ni + i

            # Local variables
            rho = density[idx]
            u = velocity_x[idx]
            v = velocity_y[idx]
            p = pressure[idx]
            c = sound_speed[idx]

            # Compute flux Jacobians
            # A = dF/dU (x-direction flux Jacobian)
            A = np.zeros((n_equations, n_equations))
            A[0, 0] = 0
            A[0, 1] = 1
            A[0, 2] = 0
            A[1, 0] = -u**2 + 0.5 * (gamma - 1) * (u**2 + v**2)
            A[1, 1] = (3 - gamma) * u
            A[1, 2] = (1 - gamma) * v
            A[2, 0] = -u * v
            A[2, 1] = v
            A[2, 2] = u

            # B = dG/dU (y-direction flux Jacobian)
            B = np.zeros((n_equations, n_equations))
            B[0, 0] = 0
            B[0, 1] = 0
            B[0, 2] = 1
            B[1, 0] = -u * v
            B[1, 1] = v
            B[1, 2] = u
            B[2, 0] = -v**2 + 0.5 * (gamma - 1) * (u**2 + v**2)
            B[2, 1] = (1 - gamma) * u
            B[2, 2] = (3 - gamma) * v

            # Compute local Jacobian contribution
            local_jacobian = np.eye(n_equations) + dt[idx] * (A / dx[idx] + B / dy[idx])

            # Add to global Jacobian
            for eq1 in range(n_equations):
                for eq2 in range(n_equations):
                    row = eq1 * n_points + idx
                    col = eq2 * n_points + idx
                    jacobian[row, col] = local_jacobian[eq1, eq2]

    # Apply boundary conditions to Jacobian
    for bc in euler_solver.boundary_conditions:
        if hasattr(bc, 'apply_jacobian'):
            bc.apply_jacobian(jacobian, euler_solver.solution)

    return jacobian.tocsr()

def compute_euler_jacobian(euler_solver, solution_vector: np.ndarray,
                         method: str = 'block') -> sp.csr_matrix:
    """
    Compute the Jacobian matrix of the Euler equations.

    Args:
        euler_solver: EulerSolver instance.
        solution_vector: Current solution vector.
        method: Method to use for Jacobian computation ('fd' or 'block').

    Returns:
        Sparse Jacobian matrix.
    """
    if method == 'fd':
        return compute_euler_jacobian_fd(euler_solver.compute_residuals, solution_vector)
    elif method == 'block':
        return compute_euler_jacobian_block(euler_solver, solution_vector)
    else:
        logger.warning(f"Unknown Jacobian method: {method}, using block method")
        return compute_euler_jacobian_block(euler_solver, solution_vector)
