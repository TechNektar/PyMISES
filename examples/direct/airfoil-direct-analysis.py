"""
PyMISES - Airfoil Direct Analysis Example

This example demonstrates the use of PyMISES for the direct analysis
of a NACA 0012 airfoil in subsonic flow.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

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
from pymises.postprocessing.visualize import plot_pressure, plot_grid, plot_boundary_layer

def run_airfoil_analysis(airfoil_name='naca0012', mach=0.5, alpha=2.0, reynolds=5e6):
    """
    Run direct analysis for an airfoil.

    Parameters
    ----------
    airfoil_name : str
        Name of the airfoil to analyze ('naca0012', 'naca4412', etc.)
    mach : float
        Freestream Mach number
    alpha : float
        Angle of attack in degrees
    reynolds : float
        Reynolds number based on chord

    Returns
    -------
    dict
        Solution dictionary containing all results
    """
    print(f"Running analysis for {airfoil_name} at M={mach}, α={alpha}°, Re={reynolds:.2e}")

    # Convert angle of attack to radians
    alpha_rad = np.radians(alpha)

    # 1. Create airfoil geometry
    if airfoil_name.lower().startswith('naca'):
        # Use NACA airfoil generator
        airfoil = AirfoilGeometry.create_naca(airfoil_name, n_points=101)
    else:
        # Load coordinates from file
        airfoil = AirfoilGeometry.import_from_file(f'{airfoil_name}.dat')

    # 2. Generate computational grid
    print("Generating computational grid...")
    grid_gen = GridGenerator(airfoil, {
        'nj': 41,  # n_normal
        'far_field_distance': 15.0  # far_field_radius
    })
    grid = grid_gen.generate_grid(grid_type='o-grid')

    # 3. Set up boundary conditions

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

    # 4. Initialize Euler solver
    print("Setting up Euler solver...")
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

    # 5. Set up the Newton solver for the Euler equations
    print("Running inviscid (Euler) solution...")
    newton_solver = NewtonSolver(
        residual_function=euler_solver.compute_residuals,
        jacobian_function=euler_solver.compute_jacobian,
        solution=euler_solver.get_solution()
    )

    # Run the Newton solver
    inviscid_solution, convergence_history = newton_solver.solve(
        max_iter=20,
        tolerance=1e-6,
        relaxation=0.7
    )

    # Update the Euler solver with the converged solution
    euler_solver.set_solution(inviscid_solution)

    # Calculate lift and update circulation
    lift, drag, moment = euler_solver.compute_forces()
    circulation = lift * 1.0  # Simplified relation, proper calculation would use Kutta-Joukowski

    # Update far-field boundary condition with calculated circulation
    farfield_bc.circulation = circulation

    # Run another iteration of the Euler solver with the updated circulation
    inviscid_solution, _ = newton_solver.solve(
        max_iter=10,
        tolerance=1e-6,
        relaxation=0.8
    )

    # 6. Set up boundary layer solver if Reynolds number is specified
    if reynolds > 0:
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
        coupled_solver.initialize(inviscid_solution)

        # Set up Newton solver for the coupled system
        coupled_newton = NewtonSolver(
            residual_function=coupled_solver.compute_residuals,
            jacobian_function=coupled_solver.compute_jacobian,
            solution=coupled_solver.get_solution()
        )

        # Run the coupled solution
        viscous_solution, viscous_convergence = coupled_newton.solve(
            max_iter=30,
            tolerance=1e-5,
            relaxation=0.6
        )

        # Get final solution
        final_solution = viscous_solution

        # Extract boundary layer properties for visualization
        bl_properties = coupled_solver.get_boundary_layer_properties()
        final_solution.update(bl_properties)
    else:
        # If inviscid only, the final solution is the inviscid solution
        final_solution = inviscid_solution

    # 7. Post-process and visualize results
    print("Post-processing results...")

    # Compute aerodynamic coefficients
    forces = euler_solver.compute_forces()
    final_solution.update({
        'cl': forces['cl'],
        'cd': forces['cd'],
        'cm': forces['cm'],
        'convergence_history': convergence_history
    })

    # Create result plots
    fig1 = plot_pressure(final_solution, airfoil)
    fig1.savefig(f'{airfoil_name}_M{mach}_a{alpha}_pressure.png', dpi=300)

    fig2 = plot_grid(grid)
    fig2.savefig(f'{airfoil_name}_M{mach}_a{alpha}_grid.png', dpi=300)

    if reynolds > 0:
        fig3 = plot_boundary_layer(final_solution, airfoil)
        fig3.savefig(f'{airfoil_name}_M{mach}_a{alpha}_boundary_layer.png', dpi=300)

    # Print summary of results
    print("\nAnalysis Results:")
    print(f"Lift coefficient (CL): {final_solution['cl']:.4f}")
    print(f"Drag coefficient (CD): {final_solution['cd']:.6f}")
    print(f"Moment coefficient (CM): {final_solution['cm']:.4f}")

    if reynolds > 0 and 'transition_location' in final_solution:
        print(f"Transition location (x/c):")
        print(f"  Upper surface: {final_solution['transition_location']['upper']:.4f}")
        print(f"  Lower surface: {final_solution['transition_location']['lower']:.4f}")

    if reynolds > 0 and 'separation_location' in final_solution:
        print(f"Separation location (x/c):")
        print(f"  Upper surface: {final_solution['separation_location']['upper']:.4f}")
        print(f"  Lower surface: {final_solution['separation_location']['lower']:.4f}")

    return final_solution

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="PyMISES Airfoil Direct Analysis Example")
    parser.add_argument('--airfoil', type=str, default='naca0012', help="Airfoil name")
    parser.add_argument('--mach', type=float, default=0.5, help="Mach number")
    parser.add_argument('--alpha', type=float, default=2.0, help="Angle of attack (degrees)")
    parser.add_argument('--reynolds', type=float, default=5e6, help="Reynolds number")
    parser.add_argument('--inviscid', action='store_true', help="Run inviscid analysis only")

    args = parser.parse_args()

    # If inviscid flag is set, set Reynolds to 0
    if args.inviscid:
        args.reynolds = 0

    # Run the analysis
    solution = run_airfoil_analysis(
        airfoil_name=args.airfoil,
        mach=args.mach,
        alpha=args.alpha,
        reynolds=args.reynolds
    )

    # Show plots
    plt.show()