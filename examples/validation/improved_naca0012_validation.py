"""
PyMISES - Improved NACA 0012 Validation Example

This example demonstrates the use of the improved PyMISES solver for validating
the NACA 0012 airfoil against experimental data.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import logging
import time

# Add the PyMISES directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from pymises.core.geometry import AirfoilGeometry
from pymises.core.grid import GridGenerator
from pymises.core.euler import EulerSolver
from pymises.core.newton import NewtonSolver
from pymises.boundary_conditions.wall import InviscidWallBC
from pymises.boundary_conditions.farfield import VortexFarfieldBC

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def analyze_naca0012(alpha, mach=0.3, grid_size=(51, 31)):
    """
    Analyze NACA 0012 airfoil at specified conditions.
    
    Parameters
    ----------
    alpha : float
        Angle of attack in degrees
    mach : float
        Freestream Mach number
    grid_size : tuple
        Grid dimensions (ni, nj)
        
    Returns
    -------
    dict
        Dictionary containing solution data
    """
    print(f"Analyzing NACA 0012 at α={alpha}°, M={mach}")
    
    # Convert angle of attack to radians
    alpha_rad = np.radians(alpha)
    
    # 1. Create airfoil geometry
    airfoil = AirfoilGeometry.create_naca('0012', n_points=101)
    
    # 2. Generate computational grid
    ni, nj = grid_size
    grid_gen = GridGenerator(airfoil, {
        'ni': ni,  # Number of points in streamwise direction
        'nj': nj,  # Number of points in normal direction
        'far_field_distance': 10.0,  # Far field boundary distance in chord lengths
        'le_clustering': 0.2,  # Leading edge clustering factor
        'te_clustering': 0.3,  # Trailing edge clustering factor
        'wall_clustering': 0.3,  # Wall clustering factor
        'clustering_method': 'tanh'  # Clustering method ('tanh', 'exp', 'sine')
    })
    
    # Generate O-grid around the airfoil
    grid = grid_gen.generate_grid(grid_type='o-grid')
    
    # 3. Set up boundary conditions
    
    # Wall boundary condition on airfoil surface
    # For O-grid, the airfoil surface is the first j-line (j=0)
    airfoil_indices = [i for i in range(grid.ni)]
    wall_bc = InviscidWallBC(airfoil_indices, normal_direction="inner")
    
    # Far-field boundary condition
    # For O-grid, the far-field boundary is the last j-line (j=nj-1)
    farfield_indices = [i + (grid.nj-1) * grid.ni for i in range(grid.ni)]
    
    # Freestream conditions (standard atmosphere at sea level)
    p_inf = 101325.0  # Pa
    T_inf = 288.15    # K
    
    # Initially set circulation to zero (will be updated)
    circulation = 0.0
    farfield_bc = VortexFarfieldBC(
        farfield_indices,
        mach_inf=mach, 
        alpha=alpha_rad,
        p0=p_inf * (1 + 0.2*mach**2)**3.5,  # Total pressure
        T0=T_inf * (1 + 0.2*mach**2),       # Total temperature
        circulation=circulation,
        airfoil_x=0.25,  # Quarter-chord position
        airfoil_y=0.0
    )
    
    # 4. Initialize Euler solver with improved settings
    euler_solver = EulerSolver(grid, {
        'cfl': 0.5,  # Conservative CFL number
        'jacobian_method': 'block',  # Use block Jacobian method
        'artificial_viscosity': 0.1  # Add artificial viscosity for stability
    })
    euler_solver.add_boundary_condition(wall_bc)
    euler_solver.add_boundary_condition(farfield_bc)
    
    # Initialize the flow field
    euler_solver.initialize(
        mach=mach,
        alpha=alpha_rad,
        p0=p_inf * (1 + 0.2*mach**2)**3.5,  # Total pressure
        T0=T_inf * (1 + 0.2*mach**2)        # Total temperature
    )
    
    # 5. Set up the Newton solver with improved settings
    newton_solver = NewtonSolver(
        config={
            'max_iterations': 30,
            'tolerance': 1e-4,  # Relaxed tolerance
            'relaxation_strategy': 'adaptive',
            'verbose': 1
        },
        residual_function=euler_solver.compute_residuals,
        jacobian_function=euler_solver.compute_jacobian,
        solution=euler_solver.get_solution_vector()
    )
    
    # 6. Run the Newton solver with conservative settings
    print("Running Euler solution...")
    start_time = time.time()
    
    solution, convergence_history = newton_solver.solve(
        max_iter=30,
        tolerance=1e-4,
        relaxation=0.3,  # Start with small relaxation
        adaptive_relaxation=True  # Use adaptive relaxation
    )
    
    end_time = time.time()
    print(f"Solution completed in {end_time - start_time:.2f} seconds")
    
    # Update the Euler solver with the converged solution
    euler_solver.set_solution_from_vector(solution)
    
    # 7. Compute forces and extract solution
    forces = euler_solver.compute_forces()
    solution_dict = euler_solver.get_solution()
    
    # Add additional data to the solution
    solution_dict.update({
        'cl': forces['cl'],
        'cd': forces['cd'],
        'cm': forces['cm'],
        'convergence_history': convergence_history,
        'grid': grid,
        'airfoil': airfoil,
        'mach': mach,
        'alpha': alpha,
        'p_inf': p_inf,
        'T_inf': T_inf
    })
    
    # 8. Extract surface pressure for visualization
    x_coords = grid.x.flatten()[airfoil_indices]
    y_coords = grid.y.flatten()[airfoil_indices]
    p = solution_dict['pressure'].flatten()[airfoil_indices]
    
    # Calculate pressure coefficient
    q_inf = 0.5 * 1.4 * p_inf * mach**2
    cp = (p - p_inf) / q_inf
    
    # Separate upper and lower surface
    upper_surface = np.where(y_coords >= 0)[0]
    lower_surface = np.where(y_coords < 0)[0]
    
    x_upper = x_coords[upper_surface]
    cp_upper = cp[upper_surface]
    x_lower = x_coords[lower_surface]
    cp_lower = cp[lower_surface]
    
    # Sort by x-coordinate
    upper_sort = np.argsort(x_upper)
    lower_sort = np.argsort(x_lower)
    
    x_upper = x_upper[upper_sort]
    cp_upper = cp_upper[upper_sort]
    x_lower = x_lower[lower_sort]
    cp_lower = cp_lower[lower_sort]
    
    solution_dict.update({
        'x_upper': x_upper,
        'cp_upper': cp_upper,
        'x_lower': x_lower,
        'cp_lower': cp_lower
    })
    
    # 9. Create visualization plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot pressure distribution
    ax1.scatter(x_upper, -cp_upper, color='blue', marker='o', label='Upper Surface')
    ax1.scatter(x_lower, -cp_lower, color='red', marker='o', label='Lower Surface')
    ax1.set_xlabel('x/c')
    ax1.set_ylabel('-Cp')
    ax1.set_title(f'NACA 0012 Pressure Distribution, α={alpha}°, M={mach}')
    ax1.grid(True)
    ax1.legend()
    
    # Plot convergence history
    ax2.semilogy(range(1, len(convergence_history) + 1), convergence_history, 'b-')
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('Residual Norm')
    ax2.set_title('Convergence History')
    ax2.grid(True)
    
    fig.tight_layout()
    fig.savefig(f'naca0012_a{alpha}_m{mach}.png', dpi=300)
    
    # Print summary of results
    print("\nAnalysis Results:")
    print(f"Lift coefficient (CL): {solution_dict['cl']:.4f}")
    print(f"Drag coefficient (CD): {solution_dict['cd']:.6f}")
    print(f"Moment coefficient (CM): {solution_dict['cm']:.4f}")
    
    return solution_dict, fig

def run_validation_sweep(alpha_range=None, mach=0.3, grid_size=(51, 31)):
    """
    Run validation sweep for NACA 0012 across a range of angles of attack.
    
    Parameters
    ----------
    alpha_range : list or None
        List of angles of attack to analyze
        If None, defaults to [-4, -2, 0, 2, 4, 6, 8]
    mach : float
        Freestream Mach number
    grid_size : tuple
        Grid dimensions (ni, nj)
        
    Returns
    -------
    dict
        Dictionary containing all results
    """
    print(f"Running NACA 0012 validation sweep at M={mach}")
    
    # Default alpha range if not specified
    if alpha_range is None:
        alpha_range = [-4, -2, 0, 2, 4, 6, 8]
    
    # Initialize results dictionary
    results = {
        'alpha': [],
        'cl': [],
        'cd': [],
        'cm': [],
        'solutions': []
    }
    
    # Run analysis for each angle of attack
    for alpha in alpha_range:
        print(f"\n{'='*60}")
        print(f"Analyzing α = {alpha}°")
        print(f"{'='*60}")
        
        # Run analysis
        solution, _ = analyze_naca0012(
            alpha=alpha,
            mach=mach,
            grid_size=grid_size
        )
        
        # Store results
        results['alpha'].append(alpha)
        results['cl'].append(solution['cl'])
        results['cd'].append(solution['cd'])
        results['cm'].append(solution['cm'])
        results['solutions'].append(solution)
    
    # Create summary plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot lift curve
    ax1.plot(results['alpha'], results['cl'], 'bo-')
    ax1.set_xlabel('Angle of Attack (deg)')
    ax1.set_ylabel('Lift Coefficient (CL)')
    ax1.set_title('NACA 0012 Lift Curve')
    ax1.grid(True)
    
    # Plot drag polar
    ax2.plot(results['cd'], results['cl'], 'ro-')
    ax2.set_xlabel('Drag Coefficient (CD)')
    ax2.set_ylabel('Lift Coefficient (CL)')
    ax2.set_title('NACA 0012 Drag Polar')
    ax2.grid(True)
    
    fig.tight_layout()
    fig.savefig(f'naca0012_validation_m{mach}.png', dpi=300)
    
    # Print summary of results
    print("\nValidation Summary:")
    print(f"{'Alpha':^10} | {'CL':^10} | {'CD':^10} | {'CM':^10}")
    print(f"{'-'*10} | {'-'*10} | {'-'*10} | {'-'*10}")
    
    for i, alpha in enumerate(results['alpha']):
        print(f"{alpha:^10.1f} | {results['cl'][i]:^10.4f} | {results['cd'][i]:^10.6f} | {results['cm'][i]:^10.4f}")
    
    return results, fig

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="PyMISES Improved NACA 0012 Validation")
    parser.add_argument('--alpha', type=float, default=None, help="Single angle of attack to analyze")
    parser.add_argument('--mach', type=float, default=0.3, help="Mach number")
    parser.add_argument('--ni', type=int, default=51, help="Number of points in streamwise direction")
    parser.add_argument('--nj', type=int, default=31, help="Number of points in normal direction")
    parser.add_argument('--sweep', action='store_true', help="Run a sweep of angles of attack")
    
    args = parser.parse_args()
    
    grid_size = (args.ni, args.nj)
    
    if args.sweep:
        # Run validation sweep
        results, fig = run_validation_sweep(
            mach=args.mach,
            grid_size=grid_size
        )
    elif args.alpha is not None:
        # Run single analysis
        solution, fig = analyze_naca0012(
            alpha=args.alpha,
            mach=args.mach,
            grid_size=grid_size
        )
    else:
        # Default to alpha=4 if no options specified
        solution, fig = analyze_naca0012(
            alpha=4.0,
            mach=args.mach,
            grid_size=grid_size
        )
    
    # Show plots
    plt.show()
