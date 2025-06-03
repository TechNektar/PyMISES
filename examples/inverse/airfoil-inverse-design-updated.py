"""
PyMISES - Airfoil Inverse Design Example

This example demonstrates the use of PyMISES for the inverse design
of an airfoil based on a target pressure distribution.
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
from pymises.boundary_conditions.inverse import PressureSpecificationBC, GeometricConstraintBC, MixedInverseBC
from pymises.postprocessing.visualize import plot_pressure, plot_geometry_comparison

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_airfoil_inverse_design(base_airfoil='naca0012', target_pressure_file=None, 
                             mach=0.5, alpha=2.0, mixed_design=True):
    """
    Run inverse design for an airfoil based on a target pressure distribution.
    
    Parameters
    ----------
    base_airfoil : str
        Name of the baseline airfoil
    target_pressure_file : str, optional
        File containing target pressure distribution
        If None, a synthetic target will be created
    mach : float
        Design Mach number
    alpha : float
        Design angle of attack in degrees
    mixed_design : bool
        If True, use mixed inverse design with fixed leading edge
        If False, use full inverse design
        
    Returns
    -------
    dict
        Solution dictionary containing all results
    """
    print(f"Running inverse design based on {base_airfoil}")
    print(f"Design condition: M={mach}, α={alpha}°")
    print(f"Design type: {'Mixed' if mixed_design else 'Full'} inverse")
    
    # Convert angle of attack to radians
    alpha_rad = np.radians(alpha)
    
    # 1. Create baseline airfoil geometry
    if base_airfoil.lower().startswith('naca'):
        # Use NACA airfoil generator
        airfoil = AirfoilGeometry.create_naca(base_airfoil, n_points=101)
    else:
        # Load coordinates from file
        airfoil = AirfoilGeometry.import_from_file(f'{base_airfoil}.dat')
    
    # 2. Generate computational grid
    print("Generating computational grid...")
    grid_gen = GridGenerator(airfoil, {
        'ni': 101,  # Number of points in streamwise direction
        'nj': 41,   # Number of points in normal direction
        'far_field_distance': 15.0,  # Far field boundary distance in chord lengths
        'le_clustering': 0.2,  # Leading edge clustering factor
        'te_clustering': 0.3,  # Trailing edge clustering factor
        'wall_clustering': 0.3,  # Wall clustering factor
        'clustering_method': 'tanh'  # Clustering method ('tanh', 'exp', 'sine')
    })
    
    # Generate O-grid around the airfoil
    grid = grid_gen.generate_grid(grid_type='o-grid')
    
    # 3. Analyze the baseline airfoil to get initial flow solution
    print("Analyzing baseline airfoil...")
    
    # Wall boundary condition on airfoil surface
    # For O-grid, the airfoil surface is the first j-line (j=0)
    airfoil_indices = [i for i in range(grid.ni)]
    wall_bc = InviscidWallBC(airfoil_indices, normal_direction='inner')
    
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
    
    # Initialize Euler solver
    euler_solver = EulerSolver(grid)
    euler_solver.add_boundary_condition(wall_bc)
    euler_solver.add_boundary_condition(farfield_bc)
    
    # Initialize the flow field
    euler_solver.initialize(
        mach=mach,
        alpha=alpha_rad,
        p0=p_inf * (1 + 0.2*mach**2)**3.5,  # Total pressure
        T0=T_inf * (1 + 0.2*mach**2)        # Total temperature
    )
    
    # Set up the Newton solver for the Euler equations
    newton_solver = NewtonSolver(
        residual_function=euler_solver.compute_residuals,
        jacobian_function=euler_solver.compute_jacobian,
        solution=euler_solver.get_solution_vector()
    )
    
    # Run the Newton solver
    baseline_solution, convergence_history = newton_solver.solve(
        max_iter=20, 
        tolerance=1e-6,
        relaxation=0.7
    )
    
    # Update the Euler solver with the converged solution
    euler_solver.set_solution_from_vector(baseline_solution)
    
    # Get the baseline solution dictionary
    baseline_solution_dict = euler_solver.get_solution()
    
    # 4. Create target pressure distribution
    print("Setting up target pressure distribution...")
    
    if target_pressure_file is not None:
        # Load target pressure from file
        target_data = np.loadtxt(target_pressure_file, delimiter=',')
        target_x = target_data[:, 0]
        target_cp = target_data[:, 1]
        
        # Map to solver grid points
        target_pressure = map_pressure_to_grid(target_x, target_cp, grid, airfoil_indices, p_inf, mach)
    else:
        # Create synthetic target by modifying baseline pressure
        # For example, accelerate flow on upper surface, decelerate on lower surface
        
        # Get baseline pressure
        baseline_pressure = baseline_solution_dict['pressure']
        
        # Create synthetic target
        target_pressure = np.copy(baseline_pressure)
        
        # Get geometric information
        x_coords = grid.x.flatten()[airfoil_indices]
        y_coords = grid.y.flatten()[airfoil_indices]
        
        # Find upper and lower surface points
        upper_surface = y_coords > 0
        lower_surface = y_coords < 0
        
        # Modify pressure on upper surface (decrease pressure)
        for i, idx in enumerate(airfoil_indices):
            if upper_surface[i] and 0.1 < x_coords[i] < 0.8:
                # Decrease pressure on upper surface by 10%
                target_pressure.flatten()[idx] *= 0.9
            elif lower_surface[i] and 0.1 < x_coords[i] < 0.8:
                # Increase pressure on lower surface by 10%
                target_pressure.flatten()[idx] *= 1.1
    
    # 5. Set up inverse design boundary conditions
    print("Setting up inverse design problem...")
    
    if mixed_design:
        # Mixed inverse design with fixed leading and trailing edges
        
        # Define segments for design
        x_coords = grid.x.flatten()[airfoil_indices]
        
        # Identify leading edge point (minimum x)
        le_idx = np.argmin(x_coords)
        
        # Fixed leading edge region (first 10% of chord)
        fixed_indices = [idx for i, idx in enumerate(airfoil_indices) 
                         if x_coords[i] < 0.1 or x_coords[i] > 0.9]
        
        # Pressure-specified region (10-90% of chord)
        design_indices = [idx for i, idx in enumerate(airfoil_indices) 
                          if 0.1 <= x_coords[i] <= 0.9]
        
        # Create target pressure values for design region
        design_pressures = target_pressure.flatten()[design_indices]
        
        # Define design segments
        design_segments = [
            {
                'type': 'fixed',
                'indices': fixed_indices
            },
            {
                'type': 'pressure',
                'indices': design_indices,
                'values': design_pressures
            }
        ]
        
        # Create mixed inverse boundary condition
        inverse_bc = MixedInverseBC(design_segments=design_segments, normal_direction='inner')
    else:
        # Full inverse design with leading edge closure constraint
        
        # Identify leading edge point (minimum x)
        x_coords = grid.x.flatten()[airfoil_indices]
        le_idx = np.argmin(x_coords)
        
        # Create pressure specification BC
        inverse_bc = PressureSpecificationBC(
            airfoil_indices, 
            target_pressure.flatten()[airfoil_indices],
            normal_direction='inner'
        )
        
        # Add leading edge geometry constraint
        le_constraint = GeometricConstraintBC(
            [airfoil_indices[le_idx]], 
            constraint_type='fixed',
            normal_direction='inner'
        )
        
        # Add trailing edge thickness constraint
        te_constraint = GeometricConstraintBC(
            [airfoil_indices[0], airfoil_indices[-1]], 
            constraint_type='thickness',
            constraint_value=0.0,  # sharp trailing edge
            normal_direction='inner'
        )
        
        # Add constraints to solver
        euler_solver.add_boundary_condition(le_constraint)
        euler_solver.add_boundary_condition(te_constraint)
    
    # 6. Run inverse design
    print("Running inverse design optimization...")
    
    # Remove the wall boundary condition
    euler_solver.remove_boundary_condition(wall_bc)
    
    # Add inverse design boundary condition
    euler_solver.add_boundary_condition(inverse_bc)
    
    # Set up Newton solver for the inverse design
    inverse_newton = NewtonSolver(
        residual_function=euler_solver.compute_residuals,
        jacobian_function=euler_solver.compute_jacobian,
        solution=euler_solver.get_solution_vector()
    )
    
    # Enable geometry mode in solver
    euler_solver.enable_geometry_mode()
    
    # Run the inverse design
    inverse_solution, inverse_convergence = inverse_newton.solve(
        max_iter=50, 
        tolerance=1e-5,
        relaxation=0.5
    )
    
    # Update the Euler solver with the converged solution
    euler_solver.set_solution_from_vector(inverse_solution)
    
    # 7. Post-process and visualize results
    print("Post-processing inverse design results...")
    
    # Extract final geometry and solution
    final_geometry = euler_solver.get_geometry()
    final_solution_dict = euler_solver.get_solution()
    
    # Compute aerodynamic coefficients
    forces = euler_solver.compute_forces()
    final_solution_dict.update({
        'cl': forces['cl'],
        'cd': forces['cd'],
        'cm': forces['cm'],
        'convergence_history': inverse_convergence,
        'grid': grid,
        'airfoil': final_geometry
    })
    
    # Create result plots
    fig1 = plot_pressure(final_solution_dict, final_geometry)
    fig1.savefig(f'{base_airfoil}_inverse_pressure.png', dpi=300)
    
    fig2 = plot_geometry_comparison(airfoil, final_geometry)
    fig2.savefig(f'{base_airfoil}_inverse_geometry.png', dpi=300)
    
    # Plot pressure comparison
    fig3, ax3 = plt.subplots(figsize=(10, 6))
    
    # Get surface coordinates for plotting
    x_coords = grid.x.flatten()[airfoil_indices]
    
    # Sort by x-coordinate for proper plotting
    sort_idx = np.argsort(x_coords)
    x_sorted = x_coords[sort_idx]
    
    # Get baseline, target, and final pressures
    baseline_p = baseline_solution_dict['pressure'].flatten()[airfoil_indices][sort_idx]
    target_p = target_pressure.flatten()[airfoil_indices][sort_idx]
    final_p = final_solution_dict['pressure'].flatten()[airfoil_indices][sort_idx]
    
    # Convert to pressure coefficient for better visualization
    p_inf = 101325.0  # Pa
    q_inf = 0.5 * 1.4 * p_inf * mach**2
    
    baseline_cp = (baseline_p - p_inf) / q_inf
    target_cp = (target_p - p_inf) / q_inf
    final_cp = (final_p - p_inf) / q_inf
    
    # Plot pressure coefficients (negative for conventional orientation)
    ax3.plot(x_sorted, -baseline_cp, 'b-', label='Baseline')
    ax3.plot(x_sorted, -target_cp, 'r--', label='Target')
    ax3.plot(x_sorted, -final_cp, 'g-', label='Final')
    ax3.set_title('Pressure Distribution Comparison')
    ax3.set_xlabel('x/c')
    ax3.set_ylabel('-Cp')
    ax3.legend()
    ax3.grid(True)
    ax3.invert_yaxis()  # Invert y-axis for conventional Cp plot
    
    fig3.savefig(f'{base_airfoil}_inverse_pressure_comparison.png', dpi=300)
    
    # Print summary of results
    print("\nInverse Design Results:")
    print(f"Baseline CL: {baseline_solution_dict['cl']:.4f}")
    print(f"Final CL: {final_solution_dict['cl']:.4f}")
    print(f"Lift change: {((final_solution_dict['cl'] - baseline_solution_dict['cl']) / baseline_solution_dict['cl'] * 100):.1f}%")
    
    print(f"Baseline CD: {baseline_solution_dict['cd']:.6f}")
    print(f"Final CD: {final_solution_dict['cd']:.6f}")
    print(f"Drag change: {((final_solution_dict['cd'] - baseline_solution_dict['cd']) / baseline_solution_dict['cd'] * 100):.1f}%")
    
    # Save final geometry to file
    final_geometry.save_to_file(f'{base_airfoil}_inverse_result.dat')
    print(f"Final geometry saved to {base_airfoil}_inverse_result.dat")
    
    return {
        'baseline_solution': baseline_solution_dict,
        'target_pressure': target_pressure,
        'final_solution': final_solution_dict,
        'baseline_geometry': airfoil,
        'final_geometry': final_geometry
    }

def map_pressure_to_grid(target_x, target_cp, grid, airfoil_indices, p_inf, mach):
    """
    Map target pressure coefficients to grid points.
    
    Parameters
    ----------
    target_x : np.ndarray
        X-coordinates of target pressure distribution
    target_cp : np.ndarray
        Pressure coefficients at target points
    grid : Grid
        Computational grid
    airfoil_indices : List[int]
        Indices of grid points on the airfoil surface
    p_inf : float
        Freestream static pressure
    mach : float
        Freestream Mach number
        
    Returns
    -------
    np.ndarray
        Target pressure values at grid points
    """
    # Get coordinates of grid points on airfoil
    grid_x = grid.x.flatten()[airfoil_indices]
    grid_y = grid.y.flatten()[airfoil_indices]
    
    # Create target pressure array (same size as grid pressures)
    target_pressure = np.full_like(grid.x, p_inf)
    
    # Convert target Cp to pressure
    # p = p_inf + 0.5 * Cp * rho_inf * V_inf^2
    # For compressible flow: p = p_inf + 0.5 * gamma * p_inf * M^2 * Cp
    target_pressure_values = p_inf * (1 + 0.7 * mach**2 * target_cp / 2)
    
    # Interpolate target pressure to grid points
    from scipy.interpolate import interp1d
    
    # Create interpolation function
    # Handle potential issues with non-monotonic x values
    if np.all(np.diff(target_x) > 0):
        # X is monotonically increasing
        interp_func = interp1d(target_x, target_pressure_values, kind='linear', 
                              bounds_error=False, fill_value='extrapolate')
    else:
        # Sort the target data by x
        sort_idx = np.argsort(target_x)
        target_x_sorted = target_x[sort_idx]
        target_p_sorted = target_pressure_values[sort_idx]
        
        # Create interpolation function
        interp_func = interp1d(target_x_sorted, target_p_sorted, kind='linear', 
                              bounds_error=False, fill_value='extrapolate')
    
    # Apply interpolation to each airfoil point
    for i, idx in enumerate(airfoil_indices):
        # Get flattened index
        flat_idx = idx
        
        # Interpolate pressure at this x-coordinate
        target_pressure.flatten()[flat_idx] = interp_func(grid_x[i])
    
    return target_pressure

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="PyMISES Airfoil Inverse Design Example")
    parser.add_argument('--airfoil', type=str, default='naca0012', help="Baseline airfoil name")
    parser.add_argument('--target', type=str, default=None, help="File with target pressure distribution")
    parser.add_argument('--mach', type=float, default=0.5, help="Design Mach number")
    parser.add_argument('--alpha', type=float, default=2.0, help="Design angle of attack (degrees)")
    parser.add_argument('--full', action='store_true', help="Use full inverse design (default is mixed)")
    
    args = parser.parse_args()
    
    # Run the inverse design
    result = run_airfoil_inverse_design(
        base_airfoil=args.airfoil,
        target_pressure_file=args.target,
        mach=args.mach,
        alpha=args.alpha,
        mixed_design=not args.full
    )
    
    # Show plots
    plt.show()
