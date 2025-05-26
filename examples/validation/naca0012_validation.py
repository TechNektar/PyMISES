"""
PyMISES - NACA 0012 Validation Example

This example validates the PyMISES solver against experimental data for the NACA 0012 airfoil
at various angles of attack and Reynolds numbers.

Reference data from:
1. Abbott, I. H., and Von Doenhoff, A. E., "Theory of Wing Sections," Dover Publications, 1959.
2. Gregory, N., and O'Reilly, C. L., "Low-Speed Aerodynamic Characteristics of NACA 0012 Airfoil
   Section, Including the Effects of Upper-Surface Roughness Simulating Hoar Frost," NPL Aero
   Report 1308, 1970.
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

def load_experimental_data(data_dir):
    """
    Load experimental data for NACA 0012.
    
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
    
    # Load lift-drag polar data
    polar_file = data_path / 'naca0012_polar.csv'
    if polar_file.exists():
        polar_data = pd.read_csv(polar_file)
    else:
        # Create dummy data if file doesn't exist
        print(f"Warning: Experimental polar data file not found: {polar_file}")
        # Create approximate data based on Abbott & Von Doenhoff
        alpha = np.linspace(-10, 10, 21)
        cl = 0.1 * alpha  # Approximate lift slope of 0.1 per degree
        cd = 0.008 + 0.0025 * (alpha / 10)**2  # Approximate drag polar
        
        polar_data = pd.DataFrame({
            'alpha': alpha,
            'cl': cl,
            'cd': cd
        })
    
    # Load pressure distribution data
    pressure_file = data_path / 'naca0012_pressure.csv'
    if pressure_file.exists():
        pressure_data = pd.read_csv(pressure_file)
    else:
        # Create dummy data if file doesn't exist
        print(f"Warning: Experimental pressure data file not found: {pressure_file}")
        pressure_data = pd.DataFrame({
            'x/c': np.linspace(0, 1, 50),
            'cp_alpha_0': np.zeros(50),
            'cp_alpha_4': np.zeros(50),
            'cp_alpha_8': np.zeros(50)
        })
    
    return {
        'polar': polar_data,
        'pressure': pressure_data
    }

def analyze_naca0012(alpha, mach=0.15, reynolds=3e6):
    """
    Analyze NACA 0012 airfoil at specified conditions.
    
    Parameters
    ----------
    alpha : float
        Angle of attack in degrees
    mach : float
        Freestream Mach number
    reynolds : float
        Reynolds number based on chord
        
    Returns
    -------
    dict
        Dictionary containing solution data
    """
    print(f"Analyzing NACA 0012 at α={alpha}°, M={mach}, Re={reynolds:.1e}")
    
    # Convert angle of attack to radians
    alpha_rad = np.radians(alpha)
    
    # 1. Create airfoil geometry
    airfoil = AirfoilGeometry.create_naca('0012', n_points=201)
    
    # 2. Generate computational grid
    grid_gen = GridGenerator(airfoil, {
        'ni': 201,  # Number of points in streamwise direction
        'nj': 61,   # Number of points in normal direction
        'far_field_distance': 15.0,  # Far field boundary distance in chord lengths
        'le_clustering': 0.15,  # Leading edge clustering factor
        'te_clustering': 0.25,  # Trailing edge clustering factor
        'wall_clustering': 0.2,  # Wall clustering factor
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
    newton_solver = NewtonSolver(
        residual_function=euler_solver.compute_residuals,
        jacobian_function=euler_solver.compute_jacobian,
        solution=euler_solver.get_solution_vector()
    )
    
    # Run the Newton solver
    inviscid_solution, euler_convergence = newton_solver.solve(
        max_iter=30, 
        tolerance=1e-6,
        relaxation=0.7
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
        max_iter=10, 
        tolerance=1e-6,
        relaxation=0.8
    )
    
    # 6. Set up boundary layer solver
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
    
    # Extract boundary layer properties
    bl_properties = coupled_solver.get_boundary_layer_properties()
    final_solution.update(bl_properties)
    
    # 7. Post-process results
    
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

def run_naca0012_validation(alpha_range=None, export_dir=None, compare_with_exp=True):
    """
    Run validation for NACA 0012 across a range of angles of attack.
    
    Parameters
    ----------
    alpha_range : list or None
        List of angles of attack to analyze
        If None, defaults to [-4, -2, 0, 2, 4, 6, 8, 10]
    export_dir : str, optional
        Directory to export results to, if None, no export is performed
    compare_with_exp : bool, optional
        If True, compare with experimental data
        
    Returns
    -------
    dict
        Dictionary containing all results
    """
    print("Running NACA 0012 Validation")
    
    # Default alpha range if not specified
    if alpha_range is None:
        alpha_range = [-4, -2, 0, 2, 4, 6, 8, 10]
    
    # Flow conditions
    mach = 0.15        # Low subsonic flow
    reynolds = 3.0e6   # Reynolds number from Abbott & Von Doenhoff
    
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
        solution = analyze_naca0012(
            alpha=alpha,
            mach=mach,
            reynolds=reynolds
        )
        
        # Store results
        results['alpha'].append(alpha)
        results['cl'].append(solution['cl'])
        results['cd'].append(solution['cd'])
        results['cm'].append(solution['cm'])
        results['solutions'].append(solution)
    
    # Compare with experimental data if requested
    if compare_with_exp:
        print("Comparing with experimental data...")
        
        # Load experimental data
        data_dir = Path(__file__).parent / 'data'
        exp_data = load_experimental_data(data_dir)
        
        # Create comparison plots
        fig_polar, (ax_cl, ax_cd, ax_cm) = plt.subplots(1, 3, figsize=(15, 5))
        
        # Plot lift curve
        ax_cl.plot(results['alpha'], results['cl'], 'bo-', label='PyMISES')
        ax_cl.set_xlabel('Angle of Attack (deg)')
        ax_cl.set_ylabel('Lift Coefficient (CL)')
        ax_cl.grid(True)
        
        # Plot drag polar
        ax_cd.plot(results['cd'], results['cl'], 'bo-', label='PyMISES')
        ax_cd.set_xlabel('Drag Coefficient (CD)')
        ax_cd.set_ylabel('Lift Coefficient (CL)')
        ax_cd.grid(True)
        
        # Plot moment curve
        ax_cm.plot(results['alpha'], results['cm'], 'bo-', label='PyMISES')
        ax_cm.set_xlabel('Angle of Attack (deg)')
        ax_cm.set_ylabel('Moment Coefficient (CM)')
        ax_cm.grid(True)
        
        # Add experimental data if available
        if 'polar' in exp_data:
            exp_alpha = exp_data['polar']['alpha']
            exp_cl = exp_data['polar']['cl']
            exp_cd = exp_data['polar']['cd']
            
            ax_cl.plot(exp_alpha, exp_cl, 'ro-', label='Experiment')
            ax_cd.plot(exp_cd, exp_cl, 'ro-', label='Experiment')
            
            if 'cm' in exp_data['polar'].columns:
                exp_cm = exp_data['polar']['cm']
                ax_cm.plot(exp_alpha, exp_cm, 'ro-', label='Experiment')
        
        # Add legends
        ax_cl.legend()
        ax_cd.legend()
        ax_cm.legend()
        
        # Add title
        fig_polar.suptitle(f'NACA 0012 Validation, M={mach}, Re={reynolds:.1e}')
        fig_polar.tight_layout()
        
        # Create pressure distribution comparison for specific angles
        key_angles = [0, 4, 8]
        key_indices = [results['alpha'].index(a) if a in results['alpha'] else None for a in key_angles]
        
        if any(idx is not None for idx in key_indices) and 'pressure' in exp_data:
            fig_cp, ax_cp = plt.subplots(len([idx for idx in key_indices if idx is not None]), 1, figsize=(10, 10), sharex=True)
            
            # Handle case with only one subplot
            if len([idx for idx in key_indices if idx is not None]) == 1:
                ax_cp = [ax_cp]
            
            plot_idx = 0
            for i, idx in enumerate(key_indices):
                if idx is None:
                    continue
                
                alpha = key_angles[i]
                solution = results['solutions'][idx]
                
                # Plot numerical results
                ax_cp[plot_idx].scatter(solution['x_upper'], -solution['cp_upper'], color='blue', marker='o', s=20, label='PyMISES (Upper)')
                ax_cp[plot_idx].scatter(solution['x_lower'], -solution['cp_lower'], color='red', marker='o', s=20, label='PyMISES (Lower)')
                
                # Plot experimental data if available
                exp_x = exp_data['pressure']['x/c']
                exp_cp_col = f'cp_alpha_{alpha}'
                
                if exp_cp_col in exp_data['pressure'].columns:
                    exp_cp = exp_data['pressure'][exp_cp_col]
                    ax_cp[plot_idx].scatter(exp_x, -exp_cp, color='black', marker='x', s=30, label='Experiment')
                
                ax_cp[plot_idx].set_ylabel('-Cp')
                ax_cp[plot_idx].set_title(f'α = {alpha}°')
                ax_cp[plot_idx].grid(True)
                ax_cp[plot_idx].legend()
                
                plot_idx += 1
            
            ax_cp[-1].set_xlabel('x/c')
            fig_cp.suptitle(f'NACA 0012 Pressure Distribution, M={mach}, Re={reynolds:.1e}')
            fig_cp.tight_layout()
        
        # Save figures if export directory is specified
        if export_dir is not None:
            export_path = Path(export_dir)
            export_path.mkdir(parents=True, exist_ok=True)
            
            fig_polar.savefig(export_path / 'naca0012_polar_comparison.png', dpi=300)
            
            if any(idx is not None for idx in key_indices) and 'pressure' in exp_data:
                fig_cp.savefig(export_path / 'naca0012_pressure_comparison.png', dpi=300)
    
    # Export results if requested
    if export_dir is not None:
        export_path = Path(export_dir)
        export_path.mkdir(parents=True, exist_ok=True)
        
        # Export solutions
        print(f"Exporting results to {export_path}...")
        
        # Export polar data
        pd.DataFrame({
            'alpha': results['alpha'],
            'cl': results['cl'],
            'cd': results['cd'],
            'cm': results['cm']
        }).to_csv(export_path / 'naca0012_polar.csv', index=False)
        
        # Export pressure distributions for key angles
        for i, alpha in enumerate(results['alpha']):
            solution = results['solutions'][i]
            
            # Combine upper and lower surface data
            pressure_data = pd.DataFrame({
                'x/c': np.concatenate([solution['x_upper'], solution['x_lower']]),
                'y/c': np.concatenate([
                    solution['airfoil'].y[solution['airfoil'].y >= 0],
                    solution['airfoil'].y[solution['airfoil'].y < 0]
                ]),
                'cp': np.concatenate([solution['cp_upper'], solution['cp_lower']])
            })
            
            # Sort by x/c
            pressure_data = pressure_data.sort_values('x/c')
            
            # Export as CSV
            pressure_data.to_csv(export_path / f'naca0012_pressure_alpha_{alpha}.csv', index=False)
    
    # Print summary of results
    print("\nNACA 0012 Validation Summary:")
    print(f"{'Alpha':^10} | {'CL':^10} | {'CD':^10} | {'CM':^10}")
    print(f"{'-'*10} | {'-'*10} | {'-'*10} | {'-'*10}")
    
    for i, alpha in enumerate(results['alpha']):
        print(f"{alpha:^10.1f} | {results['cl'][i]:^10.4f} | {results['cd'][i]:^10.6f} | {results['cm'][i]:^10.4f}")
    
    return results

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="PyMISES NACA 0012 Validation")
    parser.add_argument('--alpha-min', type=float, default=-4.0, help="Minimum angle of attack")
    parser.add_argument('--alpha-max', type=float, default=10.0, help="Maximum angle of attack")
    parser.add_argument('--alpha-step', type=float, default=2.0, help="Angle of attack step size")
    parser.add_argument('--export-dir', type=str, default='results', help="Directory to export results to")
    parser.add_argument('--no-compare', action='store_true', help="Skip comparison with experimental data")
    
    args = parser.parse_args()
    
    # Create alpha range
    alpha_range = np.arange(args.alpha_min, args.alpha_max + 0.1, args.alpha_step).tolist()
    
    # Run validation
    results = run_naca0012_validation(
        alpha_range=alpha_range,
        export_dir=args.export_dir,
        compare_with_exp=not args.no_compare
    )
    
    # Show plots
    plt.show()
