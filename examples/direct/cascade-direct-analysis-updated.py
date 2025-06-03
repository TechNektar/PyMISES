"""
PyMISES - Cascade Direct Analysis Example

This example demonstrates the use of PyMISES for the direct analysis
of a turbine or compressor cascade.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import logging

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
from pymises.postprocessing.visualize import plot_pressure, plot_grid, plot_boundary_layer, plot_mach_contours

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_cascade_analysis(blade_file, inlet_angle=30.0, outlet_angle=None, mach_in=0.5, 
                        p_ratio=0.9, reynolds=5e5, pitch=0.8, stagger=30.0):
    """
    Run direct analysis for a blade cascade.
    
    Parameters
    ----------
    blade_file : str
        File containing blade profile coordinates
    inlet_angle : float
        Inlet flow angle in degrees (from axial)
    outlet_angle : float, optional
        Outlet flow angle in degrees (from axial), if None, determined by the solver
    mach_in : float
        Inlet Mach number
    p_ratio : float
        Static pressure ratio (outlet/inlet)
    reynolds : float
        Reynolds number based on chord
    pitch : float
        Blade-to-blade spacing divided by chord
    stagger : float
        Stagger angle in degrees
        
    Returns
    -------
    dict
        Solution dictionary containing all results
    """
    print(f"Running cascade analysis with:")
    print(f"  Blade: {blade_file}")
    print(f"  Inlet angle: {inlet_angle}°")
    print(f"  Inlet Mach: {mach_in}")
    print(f"  Pressure ratio: {p_ratio}")
    print(f"  Reynolds number: {reynolds:.2e}")
    print(f"  Pitch/chord: {pitch}")
    print(f"  Stagger angle: {stagger}°")
    
    # Convert angles to radians
    inlet_angle_rad = np.radians(inlet_angle)
    stagger_rad = np.radians(stagger)
    
    # 1. Create blade geometry
    blade = BladeGeometry.import_from_file(blade_file)
    
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
        'ni': 101,  # Number of points in streamwise direction
        'nj': 41,   # Number of points in normal direction
        'far_field_distance': 5.0,  # Far field boundary distance in chord lengths
        'wake_length': 2.0,  # Wake length in chord lengths
        'le_clustering': 0.2,  # Leading edge clustering factor
        'te_clustering': 0.3,  # Trailing edge clustering factor
        'wall_clustering': 0.3,  # Wall clustering factor
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
    p_out = p_inf * p_ratio  # Static pressure at outlet
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
    inviscid_solution, convergence_history = newton_solver.solve(
        max_iter=30, 
        tolerance=1e-6,
        relaxation=0.7
    )
    
    # Update the Euler solver with the converged solution
    euler_solver.set_solution_from_vector(inviscid_solution)
    
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
        viscous_solution, viscous_convergence = coupled_newton.solve(
            max_iter=40, 
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
    
    # Compute cascade performance parameters
    performance = compute_cascade_performance(final_solution, euler_solver, inlet_angle_rad)
    final_solution.update(performance)
    
    # Add additional data to the solution for visualization
    final_solution.update({
        'grid': grid,
        'blade': blade,
        'cascade': cascade,
        'mach_in': mach_in,
        'inlet_angle': inlet_angle,
        'reynolds': reynolds,
        'pitch': pitch,
        'stagger': stagger,
        'p_inf': p_inf,
        'T_inf': T_inf,
        'convergence_history': convergence_history
    })
    
    # Create result plots
    fig1 = plot_pressure(final_solution, blade)
    fig1.savefig(f'cascade_pressure.png', dpi=300)
    
    fig2 = plot_grid(grid)
    fig2.savefig(f'cascade_grid.png', dpi=300)
    
    fig3 = plot_mach_contours(final_solution, grid)
    fig3.savefig(f'cascade_mach.png', dpi=300)
    
    if reynolds > 0 and 'displacement_thickness' in final_solution:
        fig4 = plot_boundary_layer(final_solution, blade)
        fig4.savefig(f'cascade_boundary_layer.png', dpi=300)
    
    # Print summary of results
    print("\nCascade Analysis Results:")
    print(f"Total pressure loss coefficient: {final_solution['loss_coefficient']:.6f}")
    print(f"Outlet flow angle: {np.degrees(final_solution['outlet_angle']):.2f}°")
    print(f"Diffusion factor: {final_solution['diffusion_factor']:.4f}")
    
    if reynolds > 0 and 'transition_location' in final_solution:
        print(f"Transition location (x/c):")
        print(f"  Suction surface: {final_solution['transition_location']['suction']:.4f}")
        print(f"  Pressure surface: {final_solution['transition_location']['pressure']:.4f}")
    
    if reynolds > 0 and 'separation_location' in final_solution:
        print(f"Separation location (x/c):")
        print(f"  Suction surface: {final_solution['separation_location']['suction']:.4f}")
        print(f"  Pressure surface: {final_solution['separation_location']['pressure']:.4f}")
    
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
    
    # Calculate diffusion factor (Lieblein's diffusion factor)
    # DF = 1 - V2/V1 + (V_theta1 - V_theta2)/(2*sigma*V1)
    # where sigma is solidity (chord/pitch)
    solidity = 1.0 / solution['pitch']  # chord/pitch (chord = 1.0)
    v_theta1 = v_1 * np.sin(inlet_angle)
    v_theta2 = v_2 * np.sin(outlet_angle)
    
    diffusion_factor = 1.0 - (v_2/v_1) + abs(v_theta1 - v_theta2) / (2.0 * solidity * v_1)
    
    # Calculate wake momentum thickness
    # This would be extracted from the boundary layer solution
    if 'wake_momentum_thickness' in solution:
        theta_wake = solution['wake_momentum_thickness']
    else:
        # Estimate from loss coefficient if not available
        theta_wake = loss_coefficient * (v_1/v_2) * 0.01  # simplified relation
    
    return {
        'inlet_mach': inlet_data['mach'],
        'outlet_mach': outlet_data['mach'],
        'inlet_angle': inlet_angle,
        'outlet_angle': outlet_angle,
        'loss_coefficient': loss_coefficient,
        'diffusion_factor': diffusion_factor,
        'wake_momentum_thickness': theta_wake,
        'pressure_ratio': p_2 / p_1,
        'total_pressure_ratio': p_t2 / p_t1
    }

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="PyMISES Cascade Direct Analysis Example")
    parser.add_argument('--blade', type=str, default='naca65.dat', help="Blade profile file")
    parser.add_argument('--inlet_angle', type=float, default=30.0, help="Inlet flow angle (degrees)")
    parser.add_argument('--mach_in', type=float, default=0.5, help="Inlet Mach number")
    parser.add_argument('--p_ratio', type=float, default=0.9, help="Static pressure ratio (outlet/inlet)")
    parser.add_argument('--reynolds', type=float, default=5e5, help="Reynolds number")
    parser.add_argument('--pitch', type=float, default=0.8, help="Pitch/chord ratio")
    parser.add_argument('--stagger', type=float, default=30.0, help="Stagger angle (degrees)")
    parser.add_argument('--inviscid', action='store_true', help="Run inviscid analysis only")
    
    args = parser.parse_args()
    
    # If inviscid flag is set, set Reynolds to 0
    if args.inviscid:
        args.reynolds = 0
    
    # Run the analysis
    solution = run_cascade_analysis(
        blade_file=args.blade,
        inlet_angle=args.inlet_angle,
        mach_in=args.mach_in,
        p_ratio=args.p_ratio,
        reynolds=args.reynolds,
        pitch=args.pitch,
        stagger=args.stagger
    )
    
    # Show plots
    plt.show()
