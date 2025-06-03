"""
PyMISES - Transonic Airfoil Validation Example

This example validates the PyMISES solver against experimental data for transonic flow
over the NACA 0012 and RAE 2822 airfoils, focusing on shock capturing capabilities.

Reference data from:
1. Harris, C. D., "Two-Dimensional Aerodynamic Characteristics of the NACA 0012 Airfoil
   in the Langley 8-Foot Transonic Pressure Tunnel," NASA TM-81927, 1981.
2. Cook, P. H., McDonald, M. A., and Firmin, M. C. P., "Aerofoil RAE 2822 - Pressure
   Distributions, and Boundary Layer and Wake Measurements," AGARD Advisory
   Report No. 138, 1979.
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

def load_experimental_data(data_dir, case):
    """
    Load experimental data for transonic airfoil cases.
    
    Parameters
    ----------
    data_dir : str
        Directory containing experimental data files
    case : str
        Case identifier ('naca0012_m080', 'rae2822_case9', etc.)
        
    Returns
    -------
    dict
        Dictionary containing experimental data
    """
    data_path = Path(data_dir)
    
    # Load surface pressure data
    pressure_file = data_path / f'{case}_cp.csv'
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
    
    return {
        'pressure': pressure_data
    }

def analyze_transonic_airfoil(airfoil_name, mach, alpha, reynolds, grid_density='medium'):
    """
    Analyze a transonic airfoil case.
    
    Parameters
    ----------
    airfoil_name : str
        Name of the airfoil ('naca0012', 'rae2822', etc.)
    mach : float
        Freestream Mach number
    alpha : float
        Angle of attack in degrees
    reynolds : float
        Reynolds number based on chord
    grid_density : str
        Grid density level ('coarse', 'medium', 'fine')
        
    Returns
    -------
    dict
        Dictionary containing solution data
    """
    print(f"Analyzing {airfoil_name} at M={mach}, α={alpha}°, Re={reynolds:.1e}")
    
    # Convert angle of attack to radians
    alpha_rad = np.radians(alpha)
    
    # 1. Create airfoil geometry
    if airfoil_name.lower().startswith('naca'):
        # Use NACA airfoil generator
        airfoil = AirfoilGeometry.create_naca(airfoil_name, n_points=201)
    else:
        # Load coordinates from file
        geometry_file = Path(__file__).parent / 'data' / f'{airfoil_name}.dat'
        if geometry_file.exists():
            airfoil = AirfoilGeometry.import_from_file(geometry_file)
        else:
            # Fallback to NACA 0012 if file not found
            print(f"Warning: Airfoil geometry file not found: {geometry_file}")
            print("Using NACA 0012 as a substitute")
            airfoil = AirfoilGeometry.create_naca('0012', n_points=201)
    
    # 2. Generate computational grid
    print("Generating computational grid...")
    
    # Set grid parameters based on density level
    if grid_density == 'coarse':
        ni, nj = 151, 41
        far_field_distance = 15.0
    elif grid_density == 'fine':
        ni, nj = 301, 101
        far_field_distance = 20.0
    else:  # medium (default)
        ni, nj = 201, 61
        far_field_distance = 15.0
    
    grid_gen = GridGenerator(airfoil, {
        'ni': ni,  # Number of points in streamwise direction
        'nj': nj,  # Number of points in normal direction
        'far_field_distance': far_field_distance,  # Far field boundary distance in chord lengths
        'le_clustering': 0.15,  # Leading edge clustering factor
        'te_clustering': 0.25,  # Trailing edge clustering factor
        'wall_clustering': 0.2,  # Wall clustering factor
        'clustering_method': 'tanh'  # Clustering method ('tanh', 'exp', 'sine')
    })
    
    # Generate C-grid around the airfoil for better shock capturing
    grid = grid_gen.generate_grid(grid_type='c-grid')
    
    # 3. Set up boundary conditions
    
    # Wall boundary condition on airfoil surface
    # For C-grid, the airfoil surface is typically along j=0 from i=ni/3 to i=2*ni/3
    ni, nj = grid.ni, grid.nj
    airfoil_indices = [i for i in range(ni // 3, 2 * ni // 3)]
    wall_bc = InviscidWallBC(airfoil_indices, normal_direction="inner")
    
    # Far-field boundary condition
    # For C-grid, the far-field boundary is typically along i=0, j=nj-1, and i=ni-1
    farfield_indices = [i for i in range(ni // 3)] + \
                      [i + (nj-1) * ni for i in range(ni)] + \
                      [i for i in range(2 * ni // 3, ni)]
    
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
    inviscid_solution, euler_convergence = newton_solver.solve(
        max_iter=100,  # More iterations for transonic cases
        tolerance=1e-6,
        relaxation=0.5  # Lower relaxation for stability
    )
    
    # Update the Euler solver with the converged solution
    euler_solver.set_solution_from_vector(inviscid_solution)
    
    # Calculate lift and update circulation
    inviscid_forces = euler_solver.compute_forces()
    lift = inviscid_forces['cl']
    circulation = lift * 1.0  # Simple relationship based on Kutta-Joukowski
    
    # Update far-field boundary condition with calculated circulation
    farfield_bc.circulation = circulation
    
    # Run another iteration of the Euler solver with the updated circulation
    inviscid_solution, _ = newton_solver.solve(
        max_iter=50, 
        tolerance=1e-6,
        relaxation=0.6
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
        print("Running viscous-inviscid interaction...")
        viscous_solution, viscous_convergence = coupled_newton.solve(
            max_iter=50, 
            tolerance=1e-5,
            relaxation=0.5  # Lower relaxation for stability
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
        viscous_convergence = []
    
    # 7. Post-process results
    print("Post-processing results...")
    
    # Compute final aerodynamic coefficients
    forces = euler_solver.compute_forces()
    final_solution.update({
        'cl': forces['cl'],
        'cd': forces['cd'],
        'cm': forces['cm'],
        'euler_convergence': euler_convergence,
        'viscous_convergence': viscous_convergence,
        'grid': grid,
        'airfoil': airfoil,
        'mach': mach,
        'alpha': alpha,
        'reynolds': reynolds,
        'p_inf': p_inf,
        'T_inf': T_inf
    })
    
    # Extract surface pressure
    x_coords = grid.x.flatten()[airfoil_indices]
    y_coords = grid.y.flatten()[airfoil_indices]
    p = final_solution['pressure'].flatten()[airfoil_indices]
    
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
    
    final_solution.update({
        'x_upper': x_upper,
        'cp_upper': cp_upper,
        'x_lower': x_lower,
        'cp_lower': cp_lower
    })
    
    return final_solution

def run_transonic_validation_cases(cases=None, export_dir=None, compare_with_exp=True, grid_density='medium'):
    """
    Run validation for multiple transonic airfoil cases.
    
    Parameters
    ----------
    cases : list or None
        List of case dictionaries with airfoil, mach, alpha, reynolds
        If None, defaults to standard test cases
    export_dir : str, optional
        Directory to export results to, if None, no export is performed
    compare_with_exp : bool, optional
        If True, compare with experimental data
    grid_density : str
        Grid density level ('coarse', 'medium', 'fine')
        
    Returns
    -------
    dict
        Dictionary containing all results
    """
    print("Running Transonic Airfoil Validation Cases")
    
    # Default cases if not specified
    if cases is None:
        cases = [
            {
                'name': 'naca0012_m080',
                'airfoil': 'naca0012',
                'mach': 0.80,
                'alpha': 0.0,
                'reynolds': 9.0e6
            },
            {
                'name': 'naca0012_m085',
                'airfoil': 'naca0012',
                'mach': 0.85,
                'alpha': 1.0,
                'reynolds': 9.0e6
            },
            {
                'name': 'rae2822_case9',
                'airfoil': 'rae2822',
                'mach': 0.73,
                'alpha': 3.19,
                'reynolds': 6.5e6
            },
            {
                'name': 'rae2822_case10',
                'airfoil': 'rae2822',
                'mach': 0.75,
                'alpha': 2.81,
                'reynolds': 6.2e6
            }
        ]
    
    # Initialize results dictionary
    results = {
        'cases': cases,
        'solutions': []
    }
    
    # Run analysis for each case
    for case in cases:
        print(f"\n{'='*60}")
        print(f"Analyzing case: {case['name']}")
        print(f"{'='*60}")
        
        # Run analysis
        solution = analyze_transonic_airfoil(
            airfoil_name=case['airfoil'],
            mach=case['mach'],
            alpha=case['alpha'],
            reynolds=case['reynolds'],
            grid_density=grid_density
        )
        
        # Store results
        results['solutions'].append(solution)
        
        # Compare with experimental data if requested
        if compare_with_exp:
            print(f"Comparing with experimental data for {case['name']}...")
            
            # Load experimental data
            data_dir = Path(__file__).parent / 'data'
            exp_data = load_experimental_data(data_dir, case['name'])
            
            # Create comparison plot
            fig_cp, ax_cp = plt.subplots(figsize=(10, 6))
            
            # Plot numerical results
            ax_cp.scatter(solution['x_upper'], -solution['cp_upper'], color='blue', marker='o', s=20, label='PyMISES (Upper)')
            ax_cp.scatter(solution['x_lower'], -solution['cp_lower'], color='red', marker='o', s=20, label='PyMISES (Lower)')
            
            # Plot experimental data if available
            if 'pressure' in exp_data:
                if 'x/c' in exp_data['pressure'] and 'cp_upper' in exp_data['pressure'] and 'cp_lower' in exp_data['pressure']:
                    exp_x = exp_data['pressure']['x/c']
                    exp_cp_upper = exp_data['pressure']['cp_upper']
                    exp_cp_lower = exp_data['pressure']['cp_lower']
                    
                    ax_cp.scatter(exp_x, -exp_cp_upper, color='blue', marker='x', s=30, label='Experiment (Upper)')
                    ax_cp.scatter(exp_x, -exp_cp_lower, color='red', marker='x', s=30, label='Experiment (Lower)')
                elif 'x/c' in exp_data['pressure'] and 'cp' in exp_data['pressure']:
                    # Single column of Cp values
                    exp_x = exp_data['pressure']['x/c']
                    exp_cp = exp_data['pressure']['cp']
                    
                    ax_cp.scatter(exp_x, -exp_cp, color='black', marker='x', s=30, label='Experiment')
            
            ax_cp.set_xlabel('x/c')
            ax_cp.set_ylabel('-Cp')
            ax_cp.set_title(f"{case['name']} - Pressure Coefficient Comparison")
            ax_cp.grid(True)
            ax_cp.legend()
            
            # Add text box with flow conditions
            textstr = f"Mach = {case['mach']:.3f}\nReynolds = {case['reynolds']:.1e}\nα = {case['alpha']:.2f}°"
            props = dict(boxstyle='round', facecolor='white', alpha=0.8)
            ax_cp.text(0.05, 0.95, textstr, transform=ax_cp.transAxes, fontsize=10,
                      verticalalignment='top', bbox=props)
            
            # Create Mach contour plot
            fig_mach = plot_mach_contours(solution, solution['grid'])
            fig_mach.suptitle(f"{case['name']} - Mach Number Contours")
            
            # Save figures if export directory is specified
            if export_dir is not None:
                export_path = Path(export_dir)
                export_path.mkdir(parents=True, exist_ok=True)
                
                fig_cp.savefig(export_path / f"{case['name']}_cp_comparison.png", dpi=300)
                fig_mach.savefig(export_path / f"{case['name']}_mach_contours.png", dpi=300)
        
        # Create additional visualization plots
        fig_pressure = plot_pressure(solution, solution['airfoil'])
        fig_pressure.suptitle(f"{case['name']} - Surface Pressure Distribution")
        
        # Save additional figures if export directory is specified
        if export_dir is not None:
            export_path = Path(export_dir)
            export_path.mkdir(parents=True, exist_ok=True)
            
            fig_pressure.savefig(export_path / f"{case['name']}_pressure.png", dpi=300)
            
            # Export pressure distribution data
            pd.DataFrame({
                'x/c': np.concatenate([solution['x_upper'], solution['x_lower']]),
                'y/c': np.concatenate([
                    solution['airfoil'].y[solution['airfoil'].y >= 0],
                    solution['airfoil'].y[solution['airfoil'].y < 0]
                ]),
                'cp': np.concatenate([solution['cp_upper'], solution['cp_lower']])
            }).to_csv(export_path / f"{case['name']}_pressure.csv", index=False)
    
    # Print summary of results
    print("\nTransonic Validation Summary:")
    print(f"{'Case':^15} | {'Mach':^8} | {'Alpha':^8} | {'CL':^8} | {'CD':^8} | {'CM':^8}")
    print(f"{'-'*15} | {'-'*8} | {'-'*8} | {'-'*8} | {'-'*8} | {'-'*8}")
    
    for i, case in enumerate(cases):
        solution = results['solutions'][i]
        print(f"{case['name']:^15} | {case['mach']:^8.3f} | {case['alpha']:^8.2f} | {solution['cl']:^8.4f} | {solution['cd']:^8.6f} | {solution['cm']:^8.4f}")
    
    return results

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="PyMISES Transonic Airfoil Validation")
    parser.add_argument('--export-dir', type=str, default='results',
                        help="Directory to export results to")
    parser.add_argument('--no-compare', action='store_true',
                        help="Skip comparison with experimental data")
    parser.add_argument('--grid-density', type=str, default='medium', choices=['coarse', 'medium', 'fine'],
                        help="Grid density level")
    parser.add_argument('--case', type=str, default=None,
                        help="Run specific case only (naca0012_m080, rae2822_case9, etc.)")
    
    args = parser.parse_args()
    
    # Filter cases if a specific case is requested
    if args.case is not None:
        cases = [
            {
                'name': 'naca0012_m080',
                'airfoil': 'naca0012',
                'mach': 0.80,
                'alpha': 0.0,
                'reynolds': 9.0e6
            },
            {
                'name': 'naca0012_m085',
                'airfoil': 'naca0012',
                'mach': 0.85,
                'alpha': 1.0,
                'reynolds': 9.0e6
            },
            {
                'name': 'rae2822_case9',
                'airfoil': 'rae2822',
                'mach': 0.73,
                'alpha': 3.19,
                'reynolds': 6.5e6
            },
            {
                'name': 'rae2822_case10',
                'airfoil': 'rae2822',
                'mach': 0.75,
                'alpha': 2.81,
                'reynolds': 6.2e6
            }
        ]
        
        filtered_cases = [case for case in cases if case['name'] == args.case]
        if not filtered_cases:
            print(f"Error: Case '{args.case}' not found")
            sys.exit(1)
        
        cases = filtered_cases
    else:
        cases = None
    
    # Run validation
    results = run_transonic_validation_cases(
        cases=cases,
        export_dir=args.export_dir,
        compare_with_exp=not args.no_compare,
        grid_density=args.grid_density
    )
    
    # Show plots
    plt.show()
