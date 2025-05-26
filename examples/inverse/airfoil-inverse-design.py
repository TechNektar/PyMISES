"""
PyMISES - Airfoil Inverse Design Example

This example demonstrates the use of PyMISES for the inverse design
of an airfoil based on a target pressure distribution.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add the PyMISES directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from pymises.core.geometry import AirfoilGeometry
from pymises.core.grid import GridGenerator
from pymises.core.euler import EulerSolver
from pymises.core.newton import NewtonSolver
from pymises.boundary_conditions.wall import InviscidWallBC
from pymises.boundary_conditions.farfield import VortexFarfieldBC
from pymises.boundary_conditions.inverse import PressureSpecificationBC, GeometricConstraintBC, MixedInverseBC
from pymises.postprocessing.visualize import plot_pressure, plot_geometry_comparison

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
        airfoil = AirfoilGeometry.load_from_file(f'{base_airfoil}.dat')
    
    # 2. Generate computational grid
    print("Generating computational grid...")
    grid_gen = GridGenerator(airfoil, grid_type='o-grid')
    grid = grid_gen.generate_grid(n_normal=41, n_wake=21, far_field_radius=15.0)
    
    # 3. Analyze the baseline airfoil to get initial flow solution
    print("Analyzing baseline airfoil...")
    
    # Wall boundary condition on airfoil surface
    airfoil_indices = grid.get_boundary_indices('airfoil')
    wall_bc = InviscidWallBC(airfoil_indices, normal_direction='inner')
    
    # Far-field boundary condition
    farfield_indices = grid.get_boundary_indices('farfield')
    
    # Freestream conditions (standard atmosphere at sea level)
    p_inf = 101325.0  # Pa
    T_inf = 288.15    # K
    
    # Initially set circulation to zero (will be updated)
    circulation = 0.0
    farfield_bc = VortexFarfieldBC(
        farfield_indices, mach, alpha_rad, p_inf, T_inf, circulation,
        airfoil_x=0.25, airfoil_y=0.0  # Quarter-chord position
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
        solution=euler_solver.get_solution()
    )
    
    # Run the Newton solver
    baseline_solution, _ = newton_solver.solve(
        max_iter=20, 
        tolerance=1e-6,
        relaxation=0.7
    )
    
    # Update the Euler solver with the converged solution
    euler_solver.set_solution(baseline_solution)
    
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
        baseline_pressure = baseline_solution['pressure']
        
        # Create synthetic target
        target_pressure = np.copy(baseline_pressure)
        
        # Get geometric information
        x_coords = grid.x[airfoil_indices]
        y_coords = grid.y[airfoil_indices]
        
        # Find upper and lower surface points
        upper_surface = y_coords > 0
        lower_surface = y_coords < 0
        
        # Modify pressure on upper surface (decrease pressure)
        for i, idx in enumerate(airfoil_indices):
            if upper_surface[i] and 0.1 < x_coords[i] < 0.8:
                # Decrease pressure on upper surface by 10%
                target_pressure[idx] *= 0.9
            elif lower_surface[i] and 0.1 < x_coords[i] < 0.8:
                # Increase pressure on lower surface by 10%
                target_pressure[idx] *= 1.1
    
    # 5. Set up inverse design boundary conditions
    print("Setting up inverse design problem...")
    
    if mixed_design:
        # Mixed inverse design with fixed leading and trailing edges
        
        # Define segments for design
        x_coords = grid.x[airfoil_indices]
        
        # Identify leading edge point (minimum x)
        le_idx = np.argmin(x_coords)
        
        # Fixed leading edge region (first 10% of chord)
        fixed_indices = [idx for i, idx in enumerate(airfoil_indices) 
                         if x_coords[i] < 0.1 or x_coords[i] > 0.9]
        
        # Pressure-specified region (10-90% of chord)
        design_indices = [idx for i, idx in enumerate(airfoil_indices) 
                          if 0.1 <= x_coords[i] <= 0.9]
        
        # Create target pressure values for design region
        design_pressures = target_pressure[design_indices]
        
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
        inverse_bc = MixedInverseBC(design_segments)
    else:
        # Full inverse design with leading edge closure constraint
        
        # Identify leading edge point (minimum x)
        x_coords = grid.x[airfoil_indices]
        le_idx = np.argmin(x_coords)
        
        # Create pressure specification BC
        inverse_bc = PressureSpecificationBC(airfoil_indices, target_pressure[airfoil_indices])
        
        # Add leading edge geometry constraint
        le_constraint = GeometricConstraintBC([airfoil_indices[le_idx]], 'fixed')
        
        # Add trailing edge thickness constraint
        te_constraint = GeometricConstraintBC([airfoil_indices[0], airfoil_indices[-1]], 
                                             'thickness',
                                             constraint_value=0.0)  # sharp trailing edge
    
    # 6. Run inverse design
    print("Running inverse design optimization...")
    
    # Add inverse design boundary condition
    euler_solver.add_boundary_condition(inverse_bc)
    
    # Set up Newton solver for the inverse design
    inverse_newton = NewtonSolver(
        residual_function=euler_solver.compute_residuals,
        jacobian_function=euler_solver.compute_jacobian,
        solution=euler_solver.get_solution()
    )
    
    # Enable geometry mode in solver
    euler_solver.enable_geometry_mode()
    
    # Run the inverse design
    inverse_solution, inverse_convergence = inverse_newton.solve(
        max_iter=50, 
        tolerance=1e-5,
        relaxation=0.5
    )
    
    # 7. Post-process and visualize results
    print("Post-processing inverse design results...")
    
    # Extract final geometry and solution
    final_geometry = euler_solver.get_geometry()
    final_solution = inverse_solution
    
    # Compute aerodynamic coefficients
    forces = euler_solver.compute_forces()
    final_solution.update({
        'cl': forces['cl'],
        'cd': forces['cd'],
        'cm': forces['cm'],
        'convergence_history': inverse_convergence
    })
    
    # Create result plots
    fig1 = plot_pressure(final_solution, final_geometry)
    fig1.savefig(f'{base_airfoil}_inverse_pressure.png', dpi=300)
    
    fig2 = plot_geometry_comparison(airfoil, final_geometry)
    fig2.savefig(f'{base_airfoil}_inverse_geometry.png', dpi=300)
    
    # Plot pressure comparison
    fig3, ax3 = plt.subplots(figsize=(10, 6))
    ax3.plot(x_coords, baseline_solution['pressure'][airfoil_indices], 'b-', label='Baseline')
    ax3.plot(x_coords, target_pressure[airfoil_indices], 'r--', label='Target')
    ax3.plot(x_coords, final_solution['pressure'][airfoil_indices], 'g-', label='Final')
    ax3.set_title('Pressure Distribution Comparison')
    ax3.set_xlabel('x/c')
    ax3.set_ylabel('Pressure')
    ax3.legend()
    ax3.grid(True)
    fig3.savefig(f'{base_airfoil}_inverse_pressure_comparison.png', dpi=300)
    
    # Print summary of results
    print("\nInverse Design Results:")
    print(f"Baseline CL: {baseline_solution['cl']:.4f}")
    print(f"Final CL: {final_solution['cl']:.4f}")
    print(f"Lift change: {((final_solution['cl'] - baseline_solution['cl']) / baseline_solution['cl'] * 100):.1f}%")
    
    print(f"Baseline CD: {baseline_solution['cd']:.6f}")
    print(f"Final CD: {final_solution['cd']:.6f}")
    print(f"Drag change: {((final_solution['cd'] - baseline_solution['cd']) / baseline_solution['cd'] * 100):.1f}%")
    
    # Save final geometry to file
    final_geometry.save_to_file(f'{base_airfoil}_inverse_result.dat')
    print(f"Final geometry saved to {base_airfoil}_inverse_result.dat")
    
    return {
        'baseline_solution': baseline_solution,
        'target_pressure': target_pressure,
        'final_solution': final_solution,
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
    grid_x = grid.x[airfoil_indices]
    grid_y = grid.y[airfoil_indices]
    
    # Create target pressure array (same size as grid pressures)
    target_pressure = np.full_like(grid.get_variable('pressure'), p_inf)
    
    # Convert target Cp to pressure
    # p = p_inf + 0.5 * Cp * rho_inf * V_inf^2
    target_pressure_values = p_inf * (1 + 0.7 * mach**2 * target_cp / 2)
    
    # Interpolate target pressure to grid points
    for i, idx in enumerate(airfoil_indices):
        # Find nearest point in target distribution
        # This is a simple nearest-point interpolation
        # A more accurate implementation would use proper interpolation
        closest_idx = np.argmin(np.abs(target_x - grid_x[i]))
        target_pressure[idx] = target_pressure_values[closest_idx]
    
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