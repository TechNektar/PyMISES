"""
PyMISES - VKI Turbine Cascade Validation Example

This example validates the PyMISES solver against experimental data for the VKI turbine cascade,
which is a standard test case for turbine blade aerodynamics.

Reference data from:
1. Sieverding, C. H., "The VKI Turbine Blade: Pressure Distribution, and Boundary Layer and Wake
   Measurements," VKI Technical Note 174, 1974.
2. Arts, T., Lambert de Rouvroit, M., and Rutherford, A. W., "Aero-Thermal Investigation of a
   Highly Loaded Transonic Linear Turbine Guide Vane Cascade," VKI Technical Note 174, 1990.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd
from pathlib import Path

# Add the PyMISES directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from pymises.core.geometry import BladeGeometry, CascadeGeometry
from pymises.core.grid import GridGenerator
from pymises.core.euler import EulerSolver
from pymises.core.boundary_layer import BoundaryLayerFactory
from pymises.core.coupling import CoupledSolver
from pymises.core.newton import NewtonSolver
from pymises.boundary_conditions.wall import InviscidWallBC, ViscousWallBC
from pymises.boundary_conditions.periodicity import PeriodicityBC
from pymises.boundary_conditions.farfield import SubsonicInflow, SubsonicOutflow
from pymises.postprocessing.visualize import (
    plot_pressure,
    plot_boundary_layer,
    plot_mach_contours,
    plot_convergence_history
)

def load_experimental_data(data_dir):
    """
    Load experimental data for VKI turbine cascade.
    
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
    pressure_file = data_path / 'vki_turbine_pressure.csv'
    if pressure_file.exists():
        pressure_data = pd.read_csv(pressure_file)
    else:
        # Create dummy data if file doesn't exist
        print(f"Warning: Experimental pressure data file not found: {pressure_file}")
        pressure_data = pd.DataFrame({
            's/c': np.linspace(0, 1, 50),
            'cp': np.zeros(50)
        })
    
    # Load wake data
    wake_file = data_path / 'vki_turbine_wake.csv'
    if wake_file.exists():
        wake_data = pd.read_csv(wake_file)
    else:
        # Create dummy data if file doesn't exist
        print(f"Warning: Experimental wake data file not found: {wake_file}")
        wake_data = pd.DataFrame({
            'y/pitch': np.linspace(-0.5, 0.5, 20),
            'p_total_loss': np.zeros(20),
            'velocity_ratio': np.ones(20)
        })
    
    return {
        'pressure': pressure_data,
        'wake': wake_data
    }

def run_vki_turbine_validation(export_dir=None, compare_with_exp=True):
    """
    Run validation for VKI turbine cascade.
    
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
    print("Running VKI Turbine Cascade Validation")
    
    # Flow conditions
    mach_in = 0.15       # Subsonic inlet Mach number
    mach_out = 0.85      # Transonic outlet Mach number
    inlet_angle = 0.0    # Axial inlet flow (degrees)
    reynolds = 6.0e5     # Reynolds number based on chord
    
    # Cascade geometry parameters
    pitch = 0.85         # Pitch-to-chord ratio
    stagger = 55.0       # Stagger angle (degrees)
    
    # Convert angles to radians
    inlet_angle_rad = np.radians(inlet_angle)
    stagger_rad = np.radians(stagger)
    
    # 1. Create blade geometry
    print("Loading VKI turbine blade geometry...")
    
    # Check if geometry file exists
    geometry_file = Path(__file__).parent / 'data' / 'vki_turbine.dat'
    if geometry_file.exists():
        blade = BladeGeometry.import_from_file(geometry_file)
    else:
        # Create simplified blade if VKI file not found
        print(f"Warning: VKI turbine geometry file not found: {geometry_file}")
        print("Using simplified blade geometry as a substitute")
        
        # Create a simple blade shape with camber and thickness
        x = np.linspace(0, 1, 101)
        thickness = 0.2 * (1 - x) * np.sqrt(x)  # Maximum thickness at 30% chord
        camber = 0.3 * x * (1 - x)              # Maximum camber at 50% chord
        
        # Create upper and lower surfaces
        y_upper = camber + thickness / 2
        y_lower = camber - thickness / 2
        
        # Combine into a single airfoil
        x_coords = np.concatenate([x, x[::-1]])
        y_coords = np.concatenate([y_upper, y_lower[::-1]])
        
        # Create blade geometry
        blade = BladeGeometry(np.column_stack((x_coords, y_coords)), name='simplified_vki')
    
    # Apply stagger and scaling
    blade.rotate(stagger_rad)
    blade.scale_to_chord(1.0)
    
    # Create cascade geometry
    cascade = CascadeGeometry(
        blade=blade,
        stagger_angle=stagger_rad,
        pitch=pitch,
        chord=1.0,
        n_blades=3  # Number of blades to include in the cascade
    )
    
    # 2. Generate computational grid
    print("Generating computational grid...")
    grid_gen = GridGenerator(cascade, {
        'ni': 151,  # Number of points in streamwise direction
        'nj': 51,   # Number of points in normal direction
        'far_field_distance': 2.0,  # Far field boundary distance in chord lengths
        'wake_length': 2.0,  # Wake length in chord lengths
        'le_clustering': 0.15,  # Leading edge clustering factor
        'te_clustering': 0.25,  # Trailing edge clustering factor
        'wall_clustering': 0.2,  # Wall clustering factor
        'clustering_method': 'tanh'  # Clustering method ('tanh', 'exp', 'sine')
    })
    
    # Generate cascade grid
    grid = grid_gen.generate_grid(grid_type='cascade')
    
    # 3. Set up boundary conditions
    
    # Wall boundary condition on blade surface
    # For cascade grid, the blade surface is typically along i=0 to i=ni/2
    ni, nj = grid.ni, grid.nj
    blade_indices = list(range(ni // 3, 2 * ni // 3))
    wall_bc = InviscidWallBC(blade_indices, normal_direction='inner')
    
    # Periodicity boundary condition
    # For cascade grid, the upper and lower periodic boundaries are along j=0 and j=nj-1
    upper_indices = [i * nj + (nj - 1) for i in range(ni)]
    lower_indices = [i * nj for i in range(ni)]
    
    # Create pitch vector based on stagger angle
    pitch_vector = np.array([pitch * np.sin(stagger_rad), pitch * np.cos(stagger_rad)])
    periodicity_bc = PeriodicityBC(upper_indices, lower_indices, pitch_vector=pitch_vector)
    
    # Inlet boundary condition
    # For cascade grid, the inlet boundary is typically along i=0
    inlet_indices = list(range(0, ni // 3))
    
    # Standard atmospheric conditions
    p_inf = 101325.0  # Pa
    T_inf = 288.15    # K
    
    # Total conditions at inlet
    gamma = 1.4  # Specific heat ratio
    p_total = p_inf * (1 + (gamma-1)/2 * mach_in**2)**(gamma/(gamma-1))
    T_total = T_inf * (1 + (gamma-1)/2 * mach_in**2)
    
    inlet_bc = SubsonicInflow(inlet_indices, p_total, T_total, inlet_angle_rad)
    
    # Outlet boundary condition
    # For cascade grid, the outlet boundary is typically along i=ni-1
    outlet_indices = list(range(2 * ni // 3, ni))
    
    # Static pressure at outlet based on isentropic relations
    p_out = p_inf * (1 + (gamma-1)/2 * mach_out**2)**(-gamma/(gamma-1))
    outlet_bc = SubsonicOutflow(outlet_indices, p_out)
    
    # 4. Initialize Euler solver
    print("Setting up Euler solver...")
    euler_solver = EulerSolver(grid)
    euler_solver.add_boundary_condition(wall_bc)
    euler_solver.add_boundary_condition(periodicity_bc)
    euler_solver.add_boundary_condition(inlet_bc)
    euler_solver.add_boundary_condition(outlet_bc)
    
    # Initialize the flow field
    euler_solver.initialize(
        mach=mach_in,
        alpha=inlet_angle_rad,
        p0=p_total,
        T0=T_total
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
        max_iter=50, 
        tolerance=1e-6,
        relaxation=0.7
    )
    
    # Update the Euler solver with the converged solution
    euler_solver.set_solution_from_vector(inviscid_solution)
    
    # 6. Set up boundary layer solver
    print("Setting up viscous (boundary layer) solution...")
    
    # Create boundary layer solver factory
    bl_factory = BoundaryLayerFactory(reynolds, transition_model='modified_ags')
    
    # Create viscous wall boundary condition (replaces inviscid wall BC)
    def displacement_thickness_provider(idx):
        # This will be updated by the coupled solver
        return 0.0
    
    viscous_wall_bc = ViscousWallBC(
        blade_indices, 
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
        max_iter=40, 
        tolerance=1e-5,
        relaxation=0.6
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
    
    # Compute cascade performance parameters
    performance = compute_cascade_performance(final_solution, euler_solver, inlet_angle_rad)
    final_solution.update(performance)
    
    # Add additional data to the solution for visualization
    final_solution.update({
        'euler_convergence': euler_convergence,
        'viscous_convergence': viscous_convergence,
        'grid': grid,
        'blade': blade,
        'cascade': cascade,
        'mach_in': mach_in,
        'mach_out': mach_out,
        'inlet_angle': inlet_angle,
        'reynolds': reynolds,
        'pitch': pitch,
        'stagger': stagger,
        'p_inf': p_inf,
        'T_inf': T_inf
    })
    
    # 8. Compare with experimental data if requested
    if compare_with_exp:
        print("Comparing with experimental data...")
        
        # Load experimental data
        data_dir = Path(__file__).parent / 'data'
        exp_data = load_experimental_data(data_dir)
        
        # Extract numerical results for comparison
        # Surface pressure
        blade_x = grid.x.flatten()[blade_indices]
        blade_y = grid.y.flatten()[blade_indices]
        p = final_solution['pressure'].flatten()[blade_indices]
        
        # Calculate pressure coefficient
        q_ref = 0.5 * 1.4 * p_inf * mach_in**2
        cp = (p - p_inf) / q_ref
        
        # Calculate surface distance (s/c)
        # Find leading edge point (minimum x)
        le_idx = np.argmin(blade_x)
        
        # Calculate surface distance from leading edge
        s = np.zeros_like(blade_x)
        
        # Suction side (upper surface)
        suction_indices = np.arange(le_idx, len(blade_x))
        for i in range(1, len(suction_indices)):
            idx = suction_indices[i]
            prev_idx = suction_indices[i-1]
            s[idx] = s[prev_idx] + np.sqrt((blade_x[idx] - blade_x[prev_idx])**2 + 
                                          (blade_y[idx] - blade_y[prev_idx])**2)
        
        # Pressure side (lower surface)
        pressure_indices = np.arange(le_idx, -1, -1)
        for i in range(1, len(pressure_indices)):
            idx = pressure_indices[i]
            prev_idx = pressure_indices[i-1]
            s[idx] = s[prev_idx] - np.sqrt((blade_x[idx] - blade_x[prev_idx])**2 + 
                                          (blade_y[idx] - blade_y[prev_idx])**2)
        
        # Normalize by chord
        s = s / 1.0  # chord length is 1.0
        
        # Create comparison plots
        fig_cp, ax_cp = plt.subplots(figsize=(10, 6))
        
        # Plot numerical results
        ax_cp.scatter(s, -cp, color='blue', marker='o', label='PyMISES')
        
        # Plot experimental data
        if 'pressure' in exp_data:
            exp_s = exp_data['pressure']['s/c']
            exp_cp = exp_data['pressure']['cp']
            
            ax_cp.scatter(exp_s, -exp_cp, color='red', marker='x', s=30, label='Experiment')
        
        ax_cp.set_xlabel('Surface Distance (s/c)')
        ax_cp.set_ylabel('-Cp')
        ax_cp.set_title('VKI Turbine Cascade - Pressure Coefficient Comparison')
        ax_cp.grid(True)
        ax_cp.legend()
        
        # Add text box with flow conditions
        textstr = f"Min = {mach_in:.3f}, Mout = {mach_out:.3f}\nRe = {reynolds:.1e}\nInlet Angle = {inlet_angle:.1f}°"
        props = dict(boxstyle='round', facecolor='white', alpha=0.8)
        ax_cp.text(0.05, 0.95, textstr, transform=ax_cp.transAxes, fontsize=10,
                  verticalalignment='top', bbox=props)
        
        # Wake comparison
        if 'wake' in exp_data and 'wake_data' in final_solution:
            fig_wake, ax_wake = plt.subplots(figsize=(10, 6))
            
            # Plot numerical results
            wake_y = final_solution['wake_data']['y/pitch']
            wake_loss = final_solution['wake_data']['p_total_loss']
            
            ax_wake.plot(wake_y, wake_loss, 'b-', label='PyMISES')
            
            # Plot experimental data
            exp_wake_y = exp_data['wake']['y/pitch']
            exp_wake_loss = exp_data['wake']['p_total_loss']
            
            ax_wake.plot(exp_wake_y, exp_wake_loss, 'ro-', label='Experiment')
            
            ax_wake.set_xlabel('Pitchwise Position (y/pitch)')
            ax_wake.set_ylabel('Total Pressure Loss Coefficient')
            ax_wake.set_title('VKI Turbine Cascade - Wake Profile Comparison')
            ax_wake.grid(True)
            ax_wake.legend()
        
        # Save figures if export directory is specified
        if export_dir is not None:
            export_path = Path(export_dir)
            export_path.mkdir(parents=True, exist_ok=True)
            
            fig_cp.savefig(export_path / 'vki_turbine_cp_comparison.png', dpi=300)
            
            if 'wake' in exp_data and 'wake_data' in final_solution:
                fig_wake.savefig(export_path / 'vki_turbine_wake_comparison.png', dpi=300)
    
    # 9. Create additional visualization plots
    
    # Mach contours
    fig_mach = plot_mach_contours(final_solution, grid)
    fig_mach.suptitle('VKI Turbine Cascade - Mach Number Contours')
    
    # Convergence history
    fig_conv = plot_convergence_history(euler_convergence, viscous_convergence)
    fig_conv.suptitle('VKI Turbine Cascade - Convergence History')
    
    # Save additional figures if export directory is specified
    if export_dir is not None:
        export_path = Path(export_dir)
        export_path.mkdir(parents=True, exist_ok=True)
        
        fig_mach.savefig(export_path / 'vki_turbine_mach_contours.png', dpi=300)
        fig_conv.savefig(export_path / 'vki_turbine_convergence.png', dpi=300)
    
    # 10. Print summary of results
    print("\nVKI Turbine Cascade Results:")
    print(f"Outlet flow angle: {np.degrees(final_solution['outlet_angle']):.2f}°")
    print(f"Total pressure loss coefficient: {final_solution['loss_coefficient']:.6f}")
    print(f"Blade loading coefficient: {final_solution['loading_coefficient']:.4f}")
    print(f"Zweifel coefficient: {final_solution['zweifel_coefficient']:.4f}")
    
    if 'transition_location' in final_solution:
        trans_loc = final_solution['transition_location']
        print(f"Transition location (s/c):")
        print(f"  Suction surface: {trans_loc.get('suction', 'N/A'):.4f}")
        print(f"  Pressure surface: {trans_loc.get('pressure', 'N/A'):.4f}")
    
    if 'separation_location' in final_solution:
        sep_loc = final_solution['separation_location']
        print(f"Separation location (s/c):")
        if sep_loc.get('suction', -1) > 0:
            print(f"  Suction surface: {sep_loc['suction']:.4f}")
        else:
            print(f"  Suction surface: None")
        
        if sep_loc.get('pressure', -1) > 0:
            print(f"  Pressure surface: {sep_loc['pressure']:.4f}")
        else:
            print(f"  Pressure surface: None")
    
    return final_solution

def compute_cascade_performance(solution, euler_solver, inlet_angle):
    """
    Compute performance parameters for the cascade.
    
    Parameters
    ----------
    solution : dict
        Solution dictionary
    euler_solver : EulerSolver
        Euler solver instance
    inlet_angle : float
        Inlet flow angle in radians
        
    Returns
    -------
    dict
        Dictionary containing performance parameters
    """
    # Extract flow data at inlet and outlet
    inlet_data = euler_solver.get_inlet_flow_data()
    outlet_data = euler_solver.get_outlet_flow_data()
    
    # Extract key parameters
    p_t1 = inlet_data['total_pressure']
    p_t2 = outlet_data['total_pressure']
    p_1 = inlet_data['static_pressure']
    p_2 = outlet_data['static_pressure']
    v_1 = inlet_data['velocity_magnitude']
    v_2 = outlet_data['velocity_magnitude']
    rho_1 = inlet_data['density']
    
    # Calculate outlet flow angle
    vx_2 = outlet_data['velocity_x']
    vy_2 = outlet_data['velocity_y']
    outlet_angle = np.arctan2(vy_2, vx_2)
    
    # Calculate total pressure loss coefficient
    loss_coefficient = (p_t1 - p_t2) / (p_t1 - p_1)
    
    # Calculate blade loading coefficient (Euler turbine equation)
    # For turbine: work = -v_theta_change * blade_speed
    # Loading coefficient = work / (0.5 * v_1^2)
    v_theta1 = v_1 * np.sin(inlet_angle)
    v_theta2 = v_2 * np.sin(outlet_angle)
    loading_coefficient = 2 * abs(v_theta2 - v_theta1) / v_1
    
    # Calculate Zweifel coefficient
    # Zweifel = 2 * (s/c) * cos(beta_2)^2 * (tan(beta_1) + tan(beta_2))
    # where beta_1 and beta_2 are inlet and outlet flow angles
    # s/c is pitch-to-chord ratio
    pitch_chord = solution['pitch']
    zweifel_coefficient = 2 * pitch_chord * np.cos(outlet_angle)**2 * (np.tan(inlet_angle) + np.tan(outlet_angle))
    
    # Extract wake data if available
    wake_data = {}
    if 'wake_momentum_thickness' in solution:
        # Calculate wake profile at a downstream location
        # This is a simplified approach - in a real implementation,
        # you would extract the actual wake profile from the solution
        y_pitch = np.linspace(-0.5, 0.5, 101)
        wake_center = 0.0
        wake_width = solution['wake_momentum_thickness'] * 10
        
        # Create a Gaussian wake profile
        p_total_loss = loss_coefficient * np.exp(-0.5 * ((y_pitch - wake_center) / wake_width)**2)
        
        wake_data = {
            'y/pitch': y_pitch,
            'p_total_loss': p_total_loss,
            'velocity_ratio': 1.0 - 0.5 * p_total_loss  # Simplified velocity deficit
        }
    
    performance = {
        'inlet_mach': inlet_data['mach'],
        'outlet_mach': outlet_data['mach'],
        'inlet_angle': inlet_angle,
        'outlet_angle': outlet_angle,
        'loss_coefficient': loss_coefficient,
        'loading_coefficient': loading_coefficient,
        'zweifel_coefficient': zweifel_coefficient,
        'pressure_ratio': p_2 / p_1,
        'total_pressure_ratio': p_t2 / p_t1
    }
    
    if wake_data:
        performance['wake_data'] = wake_data
    
    return performance

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="PyMISES VKI Turbine Cascade Validation")
    parser.add_argument('--export-dir', type=str, default='results',
                        help="Directory to export results to")
    parser.add_argument('--no-compare', action='store_true',
                        help="Skip comparison with experimental data")
    
    args = parser.parse_args()
    
    # Run validation
    solution = run_vki_turbine_validation(
        export_dir=args.export_dir,
        compare_with_exp=not args.no_compare
    )
    
    # Show plots
    plt.show()
