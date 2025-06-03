"""
PyMISES - Airfoil Validation Example

This example demonstrates the validation of PyMISES against experimental data
for a NACA 0012 airfoil at various angles of attack.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import logging
import pandas as pd

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

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def analyze_airfoil(airfoil_name='naca0012', mach=0.5, alpha=2.0, reynolds=5e6):
    """
    Analyze an airfoil using PyMISES.
    
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
    print(f"Analyzing {airfoil_name} at M={mach}, α={alpha}°, Re={reynolds:.2e}")
    
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
        solution=euler_solver.get_solution_vector()
    )
    
    # Run the Newton solver
    inviscid_solution, convergence_history = newton_solver.solve(
        max_iter=20,
        tolerance=1e-6,
        relaxation=0.7
    )
    
    # Update the Euler solver with the converged solution
    euler_solver.set_solution_from_vector(inviscid_solution)
    
    # Calculate lift and update circulation
    forces = euler_solver.compute_forces()
    lift = forces['cl']
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
            solution=coupled_solver.get_solution_vector()
        )
        
        # Run the coupled solution
        viscous_solution, viscous_convergence = coupled_newton.solve(
            max_iter=30,
            tolerance=1e-5,
            relaxation=0.6
        )
        
        # Update the coupled solver with the converged solution
        coupled_solver.set_solution_from_vector(viscous_solution)
        
        # Get final solution
        final_solution = coupled_solver.get_solution()
        
        # Extract boundary layer properties for visualization
        bl_properties = coupled_solver.get_boundary_layer_properties()
        final_solution.update(bl_properties)
    else:
        # If inviscid only, the final solution is the inviscid solution
        final_solution = euler_solver.get_solution()
    
    # 7. Post-process and visualize results
    print("Post-processing results...")
    
    # Compute aerodynamic coefficients
    forces = euler_solver.compute_forces()
    final_solution.update({
        'cl': forces['cl'],
        'cd': forces['cd'],
        'cm': forces['cm'],
        'convergence_history': convergence_history,
        'grid': grid,
        'airfoil': airfoil,
        'mach': mach,
        'alpha': alpha,
        'reynolds': reynolds,
        'p_inf': p_inf,
        'T_inf': T_inf
    })
    
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

def load_experimental_data(filename):
    """
    Load experimental data from a CSV file.
    
    Parameters
    ----------
    filename : str
        Path to the CSV file containing experimental data
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing experimental data
    """
    try:
        data = pd.read_csv(filename)
        return data
    except Exception as e:
        logger.error(f"Error loading experimental data: {str(e)}")
        return None

def run_validation_sweep(airfoil_name='naca0012', mach=0.5, alpha_range=None, reynolds=5e6, 
                       exp_data_file=None):
    """
    Run a validation sweep across multiple angles of attack.
    
    Parameters
    ----------
    airfoil_name : str
        Name of the airfoil to analyze
    mach : float
        Freestream Mach number
    alpha_range : list or None
        List of angles of attack to analyze
        If None, defaults to [-4, -2, 0, 2, 4, 6, 8, 10]
    reynolds : float
        Reynolds number based on chord
    exp_data_file : str or None
        Path to CSV file containing experimental data
        If None, no experimental comparison is performed
        
    Returns
    -------
    dict
        Dictionary containing all results
    """
    # Default alpha range if not specified
    if alpha_range is None:
        alpha_range = [-4, -2, 0, 2, 4, 6, 8, 10]
    
    print(f"Running validation sweep for {airfoil_name} at M={mach}, Re={reynolds:.2e}")
    print(f"Angles of attack: {alpha_range}")
    
    # Load experimental data if provided
    exp_data = None
    if exp_data_file is not None:
        exp_data = load_experimental_data(exp_data_file)
        if exp_data is not None:
            print(f"Loaded experimental data from {exp_data_file}")
    
    # Initialize results dictionary
    results = {
        'airfoil': airfoil_name,
        'mach': mach,
        'reynolds': reynolds,
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
        solution = analyze_airfoil(
            airfoil_name=airfoil_name,
            mach=mach,
            alpha=alpha,
            reynolds=reynolds
        )
        
        # Store results
        results['alpha'].append(alpha)
        results['cl'].append(solution['cl'])
        results['cd'].append(solution['cd'])
        results['cm'].append(solution['cm'])
        results['solutions'].append(solution)
        
        # Create pressure distribution plot
        fig = plot_pressure(solution, solution['airfoil'])
        fig.savefig(f'{airfoil_name}_M{mach}_a{alpha}_pressure.png', dpi=300)
        plt.close(fig)
        
        # Create boundary layer plot if viscous analysis was performed
        if reynolds > 0 and 'displacement_thickness' in solution:
            fig = plot_boundary_layer(solution, solution['airfoil'])
            fig.savefig(f'{airfoil_name}_M{mach}_a{alpha}_boundary_layer.png', dpi=300)
            plt.close(fig)
    
    # Create validation plots
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    
    # Plot lift curve
    ax1.plot(results['alpha'], results['cl'], 'bo-', label='PyMISES')
    ax1.set_xlabel('Angle of Attack (deg)')
    ax1.set_ylabel('Lift Coefficient (CL)')
    ax1.grid(True)
    
    # Plot drag polar
    ax2.plot(results['cd'], results['cl'], 'bo-', label='PyMISES')
    ax2.set_xlabel('Drag Coefficient (CD)')
    ax2.set_ylabel('Lift Coefficient (CL)')
    ax2.grid(True)
    
    # Plot moment curve
    ax3.plot(results['alpha'], results['cm'], 'bo-', label='PyMISES')
    ax3.set_xlabel('Angle of Attack (deg)')
    ax3.set_ylabel('Moment Coefficient (CM)')
    ax3.grid(True)
    
    # Add experimental data if available
    if exp_data is not None:
        if 'alpha' in exp_data.columns and 'cl' in exp_data.columns:
            ax1.plot(exp_data['alpha'], exp_data['cl'], 'ro-', label='Experiment')
        
        if 'cd' in exp_data.columns and 'cl' in exp_data.columns:
            ax2.plot(exp_data['cd'], exp_data['cl'], 'ro-', label='Experiment')
        
        if 'alpha' in exp_data.columns and 'cm' in exp_data.columns:
            ax3.plot(exp_data['alpha'], exp_data['cm'], 'ro-', label='Experiment')
    
    # Add legends
    ax1.legend()
    ax2.legend()
    ax3.legend()
    
    # Add title
    fig.suptitle(f'{airfoil_name} Validation, M={mach}, Re={reynolds:.2e}')
    fig.tight_layout()
    
    # Save validation plot
    fig.savefig(f'{airfoil_name}_M{mach}_validation.png', dpi=300)
    
    # Create a summary table
    print("\nValidation Summary:")
    print(f"{'Alpha':^10} | {'CL':^10} | {'CD':^10} | {'CM':^10}")
    print(f"{'-'*10} | {'-'*10} | {'-'*10} | {'-'*10}")
    
    for i, alpha in enumerate(results['alpha']):
        print(f"{alpha:^10.1f} | {results['cl'][i]:^10.4f} | {results['cd'][i]:^10.6f} | {results['cm'][i]:^10.4f}")
    
    return results, fig

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="PyMISES Airfoil Validation Example")
    parser.add_argument('--airfoil', type=str, default='naca0012', help="Airfoil name")
    parser.add_argument('--mach', type=float, default=0.5, help="Mach number")
    parser.add_argument('--reynolds', type=float, default=5e6, help="Reynolds number")
    parser.add_argument('--alpha_min', type=float, default=-4.0, help="Minimum angle of attack")
    parser.add_argument('--alpha_max', type=float, default=10.0, help="Maximum angle of attack")
    parser.add_argument('--alpha_step', type=float, default=2.0, help="Angle of attack step size")
    parser.add_argument('--exp_data', type=str, default=None, help="Experimental data file (CSV)")
    parser.add_argument('--inviscid', action='store_true', help="Run inviscid analysis only")
    
    args = parser.parse_args()
    
    # If inviscid flag is set, set Reynolds to 0
    if args.inviscid:
        args.reynolds = 0
    
    # Create alpha range
    alpha_range = np.arange(args.alpha_min, args.alpha_max + 0.1, args.alpha_step).tolist()
    
    # Run validation sweep
    results, fig = run_validation_sweep(
        airfoil_name=args.airfoil,
        mach=args.mach,
        alpha_range=alpha_range,
        reynolds=args.reynolds,
        exp_data_file=args.exp_data
    )
    
    # Show plots
    plt.show()
