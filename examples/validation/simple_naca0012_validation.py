"""
PyMISES - Simple NACA 0012 Validation Example

This is a simplified validation example that focuses on just the inviscid solution
with a very small grid and relaxed convergence criteria.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import logging

# Add the PyMISES directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from pymises.core.geometry import AirfoilGeometry
from pymises.core.grid import GridGenerator
from pymises.core.euler import EulerSolver
from pymises.core.newton import NewtonSolver
from pymises.boundary_conditions.wall import InviscidWallBC
from pymises.boundary_conditions.farfield import VortexFarfieldBC
from pymises.physics.thermo import Thermodynamics

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_simple_validation():
    """
    Run a simple validation case for NACA 0012 airfoil.
    """
    print("Running simple NACA 0012 validation")
    
    # Flow conditions
    mach = 0.3        # Low subsonic flow
    alpha = 2.0       # Angle of attack in degrees
    
    # Convert angle of attack to radians
    alpha_rad = np.radians(alpha)
    
    # 1. Create airfoil geometry
    print("Creating NACA 0012 airfoil geometry...")
    airfoil = AirfoilGeometry.create_naca('0012', n_points=101)
    
    # 2. Generate computational grid
    print("Generating computational grid...")
    grid_gen = GridGenerator(airfoil, {
        'ni': 51,  # Increased number of points in streamwise direction
        'nj': 31,  # Increased number of points in normal direction
        'far_field_distance': 15.0,  # Increased far field boundary distance in chord lengths
        'le_clustering': 0.1,  # Increased leading edge clustering factor
        'te_clustering': 0.2,  # Trailing edge clustering factor
        'wall_clustering': 0.1,  # Increased wall clustering factor
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
    
    # Create a thermodynamic model (needed for proper initialization)
    thermo = Thermodynamics(
        gamma=1.4,         # Ratio of specific heats
        R=287.0            # Gas constant for air in J/(kg·K)
    )
    
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
    
    # 4. Initialize Euler solver
    print("Setting up Euler solver...")
    euler_solver = EulerSolver(grid, config={
        'stability_factor': 0.3,      # Added stability factor
        'reference_length': 1.0,      # Chord length
        'thermo': thermo,             # Added thermodynamic model
        'dissipation_coefficient': 0.1  # Lower dissipation for better stability
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
    
    # 5. Set up the Newton solver for the Euler equations
    print("Running inviscid (Euler) solution...")
    newton_solver = NewtonSolver(
        config={
            'max_iterations': 50,
            'tolerance': 1e-4,
            'relaxation_strategy': 'adaptive',
            'verbose': 1
        },
        residual_function=euler_solver.compute_residuals,
        jacobian_function=euler_solver.compute_jacobian,
        solution=euler_solver.get_solution_vector()
    )
    
    # Run the Newton solver with improved settings
    inviscid_solution, euler_convergence = newton_solver.solve(
        max_iter=30,           # More iterations for better convergence
        tolerance=1e-4,        # Still relaxed but reasonable tolerance
        relaxation=0.1,        # Start with lower relaxation factor for stability
        adaptive_relaxation=True  # Enable adaptive relaxation
    )
    
    # Update the Euler solver with the converged solution
    euler_solver.set_solution_from_vector(inviscid_solution)
    
    # Get the solution dictionary
    solution = euler_solver.get_solution()
    
    # Calculate forces
    forces = euler_solver.compute_forces()
    solution.update({
        'cl': forces['cl'],
        'cd': forces['cd'],
        'cm': forces['cm'],
        'convergence_history': euler_convergence,
        'grid': grid,
        'airfoil': airfoil,
        'mach': mach,
        'alpha': alpha,
        'p_inf': p_inf,
        'T_inf': T_inf
    })
    
    # Print summary of results
    print("\nSimple Validation Results:")
    print(f"Lift coefficient (CL): {solution['cl']:.4f}")
    print(f"Drag coefficient (CD): {solution['cd']:.6f}")
    print(f"Moment coefficient (CM): {solution['cm']:.4f}")
    print(f"Final residual: {euler_convergence[-1]:.2e}")
    
    # Create a simple pressure plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Extract surface pressure
    x_coords = grid.x.flatten()[airfoil_indices]
    y_coords = grid.y.flatten()[airfoil_indices]
    p = solution['pressure'].flatten()[airfoil_indices]
    
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
    
    # Plot pressure distribution
    ax.scatter(x_upper, -cp_upper, color='blue', marker='o', label='Upper Surface')
    ax.scatter(x_lower, -cp_lower, color='red', marker='o', label='Lower Surface')
    
    ax.set_xlabel('x/c')
    ax.set_ylabel('-Cp')
    ax.set_title(f'NACA 0012 Pressure Distribution, M={mach}, α={alpha}°')
    ax.grid(True)
    ax.legend()
    
    # Add a plot of convergence history
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    ax2.semilogy(euler_convergence, 'b-')
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('Residual (log scale)')
    ax2.set_title('Convergence History')
    ax2.grid(True)
    
    # Save the figures
    fig.savefig('naca0012_simple_validation.png', dpi=300)
    fig2.savefig('naca0012_simple_convergence.png', dpi=300)
    
    # Show the plots
    plt.show()
    
    return solution

if __name__ == "__main__":
    # Run the simple validation
    solution = run_simple_validation()
