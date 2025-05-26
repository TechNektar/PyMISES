"""
PyMISES - RAE 2822 Validation Example (Case 6)

This example validates the PyMISES solver against the well-known RAE 2822 airfoil
Case 6 test case, which features transonic flow with a shock on the upper surface.

Reference data from:
1. Cook, P. H., McDonald, M. A., and Firmin, M. C. P., "Aerofoil RAE 2822 - Pressure
   Distributions, and Boundary Layer and Wake Measurements," AGARD Advisory
   Report No. 138, 1979.
2. NPARC Alliance Validation Archive, https://www.grc.nasa.gov/www/wind/valid/raetaf/raetaf.html
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd
from pathlib import Path

# Add the PyMISES directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from pymises.core.geometry import AirfoilGeometry
from pymises.core.grid import GridGenerator
from pymises.core.euler import EulerSolver
from pymises.core.boundary_layer import BoundaryLayerFactory
from pymises.core.coupling import CoupledSolver
from pymises.core.newton import NewtonSolver
from pymises.boundary_conditions.wall import InviscidWallBC, ViscousWallBC
from pymises.boundary_conditions.farfield import VortexFarfieldBC
from pymises.postprocessing.visualize import (
    plot_pressure,
    plot_boundary_layer,
    plot_mach_contours,
    plot_convergence_history
)
from pymises.postprocessing.export import (
    export_solution_to_vtk,
    export_pressure_distribution,
    export_performance_report
)

def load_experimental_data(data_dir):
    """
    Load experimental data for RAE 2822 Case 6.

    Parameters
    ----------
    data_dir : str
        Directory containing experimental data files

    Returns
    -------
    dict
        Dictionary containing experimental data
    """
    data_path = Path(data_dir)

    # Load surface pressure data
    pressure_file = data_path / 'rae2822_case6_cp.csv'
    if pressure_file.exists():
        pressure_data = pd.read_csv(pressure_file)
    else:
        # Create dummy data if file doesn't exist
        print(f"Warning: Experimental pressure data file not found: {pressure_file}")
        pressure_data = pd.DataFrame({
            'x/c': np.linspace(0, 1, 50),
            'cp_upper': np.zeros(50),
            'cp_lower': np.zeros(50)
        })

    # Load boundary layer data
    bl_file = data_path / 'rae2822_case6_bl.csv'
    if bl_file.exists():
        bl_data = pd.read_csv(bl_file)
    else:
        # Create dummy data if file doesn't exist
        print(f"Warning: Experimental boundary layer data file not found: {bl_file}")
        bl_data = pd.DataFrame({
            'x/c': np.linspace(0, 1, 20),
            'delta_star_upper': np.zeros(20),
            'theta_upper': np.zeros(20),
            'H_upper': np.zeros(20),
            'cf_upper': np.zeros(20)
        })

    return {
        'pressure': pressure_data,
        'boundary_layer': bl_data
    }

def run_rae2822_case6_validation(export_dir=None, compare_with_exp=True):
    """
    Run validation for RAE 2822 Case 6.

    Parameters
    ----------
    export_dir : str, optional
        Directory to export results to, if None, no export is performed
    compare_with_exp : bool, optional
        If True, compare with experimental data

    Returns
    -------
    dict
        Dictionary containing solution data
    """
    print("Running RAE 2822 Case 6 Validation")

    # Case 6 flow conditions
    mach = 0.725       # Freestream Mach number
    alpha = 2.31       # Angle of attack in degrees
    reynolds = 6.5e6   # Reynolds number

    # Convert angle of attack to radians
    alpha_rad = np.radians(alpha)

    # 1. Create airfoil geometry
    print("Loading RAE 2822 airfoil geometry...")

    # Check if geometry file exists
    geometry_file = Path(__file__).parent / 'data' / 'rae2822.dat'
    if geometry_file.exists():
        airfoil = AirfoilGeometry.import_from_file(str(geometry_file))
    else:
        # Create simplified NACA 0012 if RAE 2822 file not found
        print(f"Warning: RAE 2822 geometry file not found: {geometry_file}")
        print("Using NACA 0012 as a substitute")
        airfoil = AirfoilGeometry.create_naca('0012', n_points=201)

    # 2. Generate computational grid
    print("Generating computational grid...")
    grid_gen = GridGenerator(airfoil, {
        'ni': 51,  # Number of points in streamwise direction
        'nj': 31,   # Number of points in normal direction
        'far_field_distance': 20.0,  # Far field boundary distance in chord lengths
        'le_clustering': 0.15,  # Leading edge clustering factor
        'te_clustering': 0.25,  # Trailing edge clustering factor
        'wall_clustering': 0.2,  # Wall clustering factor
        'clustering_method': 'tanh'  # Clustering method ('tanh', 'exp', 'sine')
    })

    # Generate O-grid around the airfoil (C-grid not fully implemented)
    grid = grid_gen.generate_grid(grid_type='o-grid')

    # 3. Set up boundary conditions

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
        farfield_indices, mach, alpha_rad, p_inf, T_inf, circulation,
        airfoil_x=0.25, airfoil_y=0.0  # Quarter-chord position
    )

    # 4. Initialize Euler solver
    print("Setting up Euler solver...")
    euler_solver = EulerSolver(grid, {
        'cfl': 0.5,  # Conservative CFL number
        'jacobian_method': 'fd',  # Use finite difference Jacobian method
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

    # 5. Set up the Newton solver for the Euler equations
    print("Running inviscid (Euler) solution...")
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

    # Run the Newton solver
    inviscid_solution, euler_convergence = newton_solver.solve(
        max_iter=30,
        tolerance=1e-4,
        relaxation=0.5,
        adaptive_relaxation=True
    )

    # Update the Euler solver with the converged solution
    euler_solver.set_solution_from_vector(inviscid_solution)

    # Get the solution dictionary
    inviscid_solution_dict = euler_solver.get_solution()

    # Calculate lift and update circulation
    forces = euler_solver.compute_forces()
    lift = forces['cl']
    circulation = lift * 1.0  # Simple relationship based on Kutta-Joukowski

    # Update far-field boundary condition with calculated circulation
    farfield_bc.circulation = circulation

    # Run another iteration of the Euler solver with the updated circulation
    inviscid_solution, _ = newton_solver.solve(
        max_iter=10,
        tolerance=1e-4,
        relaxation=0.5,
        adaptive_relaxation=True
    )

    # Update the Euler solver with the converged solution
    euler_solver.set_solution_from_vector(inviscid_solution)

    # Get the updated solution dictionary
    inviscid_solution_dict = euler_solver.get_solution()

    # 6. Set up boundary layer solver
    print("Setting up viscous (boundary layer) solution...")

    # Create boundary layer solver factory
    bl_factory = BoundaryLayerFactory(reynolds, transition_model='modified_ags')

    # Create viscous wall boundary condition (replaces inviscid wall BC)
    def displacement_thickness_provider(idx):
        # This will be updated by the coupled solver
        return 0.0

    viscous_wall_bc = ViscousWallBC(
        airfoil_indices,
        normal_direction='inner',
        displacement_thickness_provider=displacement_thickness_provider
    )

    # Replace inviscid wall BC with viscous wall BC
    euler_solver.remove_boundary_condition(wall_bc)
    euler_solver.add_boundary_condition(viscous_wall_bc)

    # Create coupled solver
    coupled_solver = CoupledSolver(euler_solver, bl_factory)

    # Initialize boundary layer solver with inviscid solution
    coupled_solver.initialize(inviscid_solution_dict)

    # Set up Newton solver for the coupled system
    coupled_newton = NewtonSolver(
        config={
            'max_iterations': 30,
            'tolerance': 1e-4,  # Relaxed tolerance
            'relaxation_strategy': 'adaptive',
            'verbose': 1
        },
        residual_function=coupled_solver.compute_residuals,
        jacobian_function=coupled_solver.compute_jacobian,
        solution=coupled_solver.get_solution_vector()
    )

    # Run the coupled solution
    print("Running viscous-inviscid interaction...")
    viscous_solution, viscous_convergence = coupled_newton.solve(
        max_iter=30,
        tolerance=1e-4,
        relaxation=0.5,
        adaptive_relaxation=True
    )

    # Update the coupled solver with the converged solution
    coupled_solver.set_solution_from_vector(viscous_solution)

    # Get final solution
    final_solution = coupled_solver.get_solution()

    # Extract boundary layer properties for validation
    bl_properties = coupled_solver.get_boundary_layer_properties()
    final_solution.update(bl_properties)

    # 7. Post-process and visualize results
    print("Post-processing results...")

    # Compute final aerodynamic coefficients
    forces = euler_solver.compute_forces()
    final_solution.update({
        'cl': forces['cl'],
        'cd': forces['cd'],
        'cm': forces['cm'],
        'euler_convergence': euler_convergence,
        'viscous_convergence': viscous_convergence
    })

    # 8. Compare with experimental data if requested
    if compare_with_exp:
        print("Comparing with experimental data...")

        # Load experimental data
        data_dir = Path(__file__).parent / 'data'
        exp_data = load_experimental_data(data_dir)

        # Extract numerical results for comparison
        # Surface pressure
        surface_indices = airfoil.get_surface_indices()
        x = final_solution['grid_x'][surface_indices]
        p = final_solution['pressure'][surface_indices]

        # Calculate pressure coefficient
        cp = (p - p_inf) / (0.5 * p_inf * 1.4 * mach**2)

        # Separate upper and lower surface
        upper_surface = np.where(final_solution['grid_y'][surface_indices] >= 0)[0]
        lower_surface = np.where(final_solution['grid_y'][surface_indices] < 0)[0]

        x_upper = x[upper_surface]
        cp_upper = cp[upper_surface]
        x_lower = x[lower_surface]
        cp_lower = cp[lower_surface]

        # Sort by x-coordinate
        upper_sort = np.argsort(x_upper)
        lower_sort = np.argsort(x_lower)

        x_upper = x_upper[upper_sort]
        cp_upper = cp_upper[upper_sort]
        x_lower = x_lower[lower_sort]
        cp_lower = cp_lower[lower_sort]

        # Create comparison plots
        fig_cp, ax_cp = plt.subplots(figsize=(10, 6))

        # Plot numerical results
        ax_cp.scatter(x_upper, -cp_upper, color='blue', marker='o', label='PyMISES (Upper)')
        ax_cp.scatter(x_lower, -cp_lower, color='red', marker='o', label='PyMISES (Lower)')

        # Plot experimental data
        exp_x = exp_data['pressure']['x/c']
        exp_cp_upper = exp_data['pressure']['cp_upper']
        exp_cp_lower = exp_data['pressure']['cp_lower']

        ax_cp.scatter(exp_x, -exp_cp_upper, color='blue', marker='x', s=30, label='Experiment (Upper)')
        ax_cp.scatter(exp_x, -exp_cp_lower, color='red', marker='x', s=30, label='Experiment (Lower)')

        ax_cp.set_xlabel('x/c')
        ax_cp.set_ylabel('-Cp')
        ax_cp.set_title('RAE 2822 Case 6 - Pressure Coefficient Comparison')
        ax_cp.grid(True)
        ax_cp.legend()

        # Add text box with flow conditions
        textstr = f"Mach = {mach:.3f}\nReynolds = {reynolds:.1e}\nα = {alpha:.2f}°"
        props = dict(boxstyle='round', facecolor='white', alpha=0.8)
        ax_cp.text(0.05, 0.95, textstr, transform=ax_cp.transAxes, fontsize=10,
                  verticalalignment='top', bbox=props)

        # Boundary layer comparison
        if 'delta_star' in final_solution and 'theta' in final_solution:
            fig_bl, ax_bl = plt.subplots(figsize=(10, 6))

            # Plot numerical results
            delta_star_upper = final_solution['delta_star'][upper_surface][upper_sort]
            theta_upper = final_solution['theta'][upper_surface][upper_sort]
            H_upper = delta_star_upper / theta_upper

            ax_bl.plot(x_upper, H_upper, 'b-', label='PyMISES H (Upper)')

            # Plot experimental data if available
            if 'boundary_layer' in exp_data:
                exp_x_bl = exp_data['boundary_layer']['x/c']
                exp_H_upper = exp_data['boundary_layer']['H_upper']

                ax_bl.scatter(exp_x_bl, exp_H_upper, color='blue', marker='x', s=30, label='Experiment H (Upper)')

            ax_bl.set_xlabel('x/c')
            ax_bl.set_ylabel('Shape Factor (H)')
            ax_bl.set_title('RAE 2822 Case 6 - Boundary Layer Shape Factor')
            ax_bl.grid(True)
            ax_bl.legend()

        # Save figures if export directory is specified
        if export_dir is not None:
            export_path = Path(export_dir)
            export_path.mkdir(parents=True, exist_ok=True)

            fig_cp.savefig(export_path / 'rae2822_case6_cp_comparison.png', dpi=300)

            if 'delta_star' in final_solution and 'theta' in final_solution:
                fig_bl.savefig(export_path / 'rae2822_case6_bl_comparison.png', dpi=300)

    # 9. Export results if requested
    if export_dir is not None:
        export_path = Path(export_dir)
        export_path.mkdir(parents=True, exist_ok=True)

        # Export solutions
        print(f"Exporting results to {export_path}...")

        # Export pressure distribution
        export_pressure_distribution(final_solution, airfoil, export_path / 'rae2822_case6_pressure.csv')

        # Export solution to VTK for visualization
        export_solution_to_vtk(final_solution, export_path / 'rae2822_case6_solution.vtk')

        # Export boundary layer data
        if 'delta_star' in final_solution and 'theta' in final_solution:
            bl_data = {
                'x': final_solution['grid_x'][surface_indices],
                'y': final_solution['grid_y'][surface_indices],
                'delta_star': final_solution['delta_star'],
                'theta': final_solution['theta'],
                'cf': final_solution.get('cf', np.zeros_like(final_solution['delta_star']))
            }

            # Export as CSV
            pd.DataFrame(bl_data).to_csv(export_path / 'rae2822_case6_bl_data.csv', index=False)

        # Export convergence history
        pd.DataFrame({
            'iteration': np.arange(1, len(euler_convergence) + 1),
            'euler_residual': euler_convergence,
            'viscous_residual': viscous_convergence[:len(euler_convergence)] if len(viscous_convergence) >= len(euler_convergence) else np.pad(viscous_convergence, (0, len(euler_convergence) - len(viscous_convergence)), 'constant', constant_values=np.nan)
        }).to_csv(export_path / 'rae2822_case6_convergence.csv', index=False)

        # Export performance report
        export_performance_report(final_solution, export_path / 'rae2822_case6_report.txt')

    # 10. Print summary of results
    print("\nRAE 2822 Case 6 Results:")
    print(f"Lift coefficient (CL): {final_solution['cl']:.4f}")
    print(f"Drag coefficient (CD): {final_solution['cd']:.6f}")
    print(f"Moment coefficient (CM): {final_solution['cm']:.4f}")

    if 'transition_location' in final_solution:
        trans_loc = final_solution['transition_location']
        print(f"Transition location (x/c):")
        print(f"  Upper surface: {trans_loc.get('upper', 'N/A'):.4f}")
        print(f"  Lower surface: {trans_loc.get('lower', 'N/A'):.4f}")

    if 'separation_location' in final_solution:
        sep_loc = final_solution['separation_location']
        print(f"Separation location (x/c):")
        if sep_loc.get('upper', -1) > 0:
            print(f"  Upper surface: {sep_loc['upper']:.4f}")
        else:
            print(f"  Upper surface: None")

        if sep_loc.get('lower', -1) > 0:
            print(f"  Lower surface: {sep_loc['lower']:.4f}")
        else:
            print(f"  Lower surface: None")

    return final_solution

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="PyMISES RAE 2822 Case 6 Validation")
    parser.add_argument('--export-dir', type=str, default='results',
                        help="Directory to export results to")
    parser.add_argument('--no-compare', action='store_true',
                        help="Skip comparison with experimental data")

    args = parser.parse_args()

    # Run validation
    solution = run_rae2822_case6_validation(
        export_dir=args.export_dir,
        compare_with_exp=not args.no_compare
    )

    # Show plots
    plt.show()