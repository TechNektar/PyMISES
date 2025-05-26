"""
PyMISES - Modal Design Optimization Example

This example demonstrates the use of PyMISES for aerodynamic optimization
of an airfoil using modal shape functions.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import logging
import time
from scipy.optimize import minimize

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
from pymises.boundary_conditions.inverse import ModalInverseDesignBC
from pymises.postprocessing.visualize import plot_pressure, plot_geometry_comparison

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ModalDesignOptimizer:
    """
    Class for airfoil optimization using modal shape functions.
    """
    
    def __init__(self, base_airfoil='naca0012', mach=0.5, alpha=2.0, reynolds=5e6,
                 n_modes=6, objective='cl_cd'):
        """
        Initialize modal design optimizer.
        
        Parameters
        ----------
        base_airfoil : str
            Name of the baseline airfoil
        mach : float
            Freestream Mach number
        alpha : float
            Angle of attack in degrees
        reynolds : float
            Reynolds number based on chord
        n_modes : int
            Number of modal shape functions
        objective : str
            Objective function ('cl', 'cd', 'cl_cd', 'cm')
        """
        self.base_airfoil_name = base_airfoil
        self.mach = mach
        self.alpha = alpha
        self.alpha_rad = np.radians(alpha)
        self.reynolds = reynolds
        self.n_modes = n_modes
        self.objective = objective
        
        # Load baseline airfoil
        if base_airfoil.lower().startswith('naca'):
            self.base_airfoil = AirfoilGeometry.create_naca(base_airfoil, n_points=101)
        else:
            self.base_airfoil = AirfoilGeometry.load_from_file(f'{base_airfoil}.dat')
        
        # Create modal shape functions
        self._create_modal_functions()
        
        # Initialize optimization history
        self.history = {
            'iterations': [],
            'objective': [],
            'cl': [],
            'cd': [],
            'cm': [],
            'design_vars': []
        }
        
        # Counter for function evaluations
        self.eval_count = 0
        
        # Best design found
        self.best_design = None
        self.best_objective = float('inf') if objective == 'cd' else -float('inf')
        self.best_airfoil = None
        self.best_solution = None
    
    def _create_modal_functions(self):
        """
        Create modal shape functions for airfoil parameterization.
        
        This uses Hicks-Henne functions as the modal basis.
        """
        # Create separate modes for thickness and camber
        self.n_thickness_modes = self.n_modes // 2
        self.n_camber_modes = self.n_modes - self.n_thickness_modes
        
        # Thickness mode locations (concentrated near leading edge)
        self.thickness_locs = np.logspace(-1, 0, self.n_thickness_modes) - 0.1
        self.thickness_locs = np.clip(self.thickness_locs, 0.05, 0.95)
        
        # Camber mode locations (distributed along chord)
        self.camber_locs = np.linspace(0.1, 0.9, self.n_camber_modes)
        
        # Total number of design variables
        self.n_total_vars = self.n_modes
    
    def _apply_modes(self, design_vars):
        """
        Apply modal shape functions to modify the airfoil.
        
        Parameters
        ----------
        design_vars : np.ndarray
            Design variables controlling modal amplitudes
            
        Returns
        -------
        AirfoilGeometry
            Modified airfoil geometry
        """
        # Split design variables for thickness and camber
        thickness_vars = design_vars[:self.n_thickness_modes]
        camber_vars = design_vars[self.n_thickness_modes:]
        
        # Get baseline coordinates
        x = self.base_airfoil.x
        y = self.base_airfoil.y
        
        # Extract baseline thickness and camber
        upper_idx = np.where(y >= 0)[0]
        lower_idx = np.where(y < 0)[0]
        
        # Compute baseline camber line
        camber = np.zeros_like(x)
        thickness = np.zeros_like(x)
        
        # For points with matching x-coordinates on upper and lower surfaces
        for i in range(len(x)):
            # Find matching point on opposite surface with same x-coordinate
            if y[i] >= 0:  # Upper surface
                # Find matching point on lower surface
                match_idx = np.argmin(np.abs(x[lower_idx] - x[i]))
                if np.abs(x[lower_idx[match_idx]] - x[i]) < 0.01:
                    # Compute camber and thickness
                    camber[i] = (y[i] + y[lower_idx[match_idx]]) / 2
                    thickness[i] = np.abs(y[i] - y[lower_idx[match_idx]]) / 2
            else:  # Lower surface
                # Find matching point on upper surface
                match_idx = np.argmin(np.abs(x[upper_idx] - x[i]))
                if np.abs(x[upper_idx[match_idx]] - x[i]) < 0.01:
                    # Compute camber and thickness
                    camber[i] = (y[upper_idx[match_idx]] + y[i]) / 2
                    thickness[i] = np.abs(y[upper_idx[match_idx]] - y[i]) / 2
        
        # Apply thickness modes
        new_thickness = np.copy(thickness)
        for i, loc in enumerate(self.thickness_locs):
            # Hicks-Henne bump function: y = sin(pi * x^t)^n
            # where t controls the location of the bump maximum
            t = np.log(0.5) / np.log(loc)
            bump = np.sin(np.pi * x**t)**2
            new_thickness += thickness_vars[i] * bump
        
        # Apply camber modes
        new_camber = np.copy(camber)
        for i, loc in enumerate(self.camber_locs):
            t = np.log(0.5) / np.log(loc)
            bump = np.sin(np.pi * x**t)**2
            new_camber += camber_vars[i] * bump
        
        # Reconstruct airfoil from modified thickness and camber
        new_y = np.zeros_like(y)
        for i in range(len(x)):
            if y[i] >= 0:  # Upper surface
                new_y[i] = new_camber[i] + new_thickness[i]
            else:  # Lower surface
                new_y[i] = new_camber[i] - new_thickness[i]
        
        # Create new airfoil geometry
        new_coords = np.column_stack((x, new_y))
        new_airfoil = AirfoilGeometry(coords=new_coords, name=f"{self.base_airfoil_name}_modal")
        
        return new_airfoil
    
    def analyze_airfoil(self, airfoil):
        """
        Analyze an airfoil using PyMISES.
        
        Parameters
        ----------
        airfoil : AirfoilGeometry
            Airfoil geometry to analyze
            
        Returns
        -------
        dict
            Solution dictionary with aerodynamic coefficients
        """
        # Generate computational grid
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
            mach_inf=self.mach, 
            alpha=self.alpha_rad,
            p0=p_inf * (1 + 0.2*self.mach**2)**3.5,  # Total pressure
            T0=T_inf * (1 + 0.2*self.mach**2),       # Total temperature
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
            mach=self.mach,
            alpha=self.alpha_rad,
            p0=p_inf * (1 + 0.2*self.mach**2)**3.5,  # Total pressure
            T0=T_inf * (1 + 0.2*self.mach**2)        # Total temperature
        )
        
        # Set up the Newton solver for the Euler equations
        newton_solver = NewtonSolver(
            residual_function=euler_solver.compute_residuals,
            jacobian_function=euler_solver.compute_jacobian,
            solution=euler_solver.get_solution_vector()
        )
        
        # Run the Newton solver
        inviscid_solution, _ = newton_solver.solve(
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
        
        # If Reynolds number is specified, run viscous analysis
        if self.reynolds > 0:
            # Create boundary layer solver factory
            bl_factory = BoundaryLayerFactory(self.reynolds, transition_model='modified_ags')
            
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
            viscous_solution, _ = coupled_newton.solve(
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
        
        # Compute aerodynamic coefficients
        forces = euler_solver.compute_forces()
        final_solution.update({
            'cl': forces['cl'],
            'cd': forces['cd'],
            'cm': forces['cm'],
            'grid': grid,
            'airfoil': airfoil
        })
        
        return final_solution
    
    def objective_function(self, design_vars):
        """
        Objective function for optimization.
        
        Parameters
        ----------
        design_vars : np.ndarray
            Design variables controlling modal amplitudes
            
        Returns
        -------
        float
            Objective function value (to be minimized)
        """
        # Apply design variables to create new airfoil
        new_airfoil = self._apply_modes(design_vars)
        
        # Analyze the airfoil
        try:
            solution = self.analyze_airfoil(new_airfoil)
            
            # Extract aerodynamic coefficients
            cl = solution['cl']
            cd = solution['cd']
            cm = solution['cm']
            
            # Calculate objective based on specified type
            if self.objective == 'cl':
                obj_value = -cl  # Negative for maximization
            elif self.objective == 'cd':
                obj_value = cd   # Minimize drag
            elif self.objective == 'cl_cd':
                obj_value = -cl/cd  # Negative for maximization of L/D
            elif self.objective == 'cm':
                obj_value = abs(cm)  # Minimize moment magnitude
            else:
                obj_value = -cl/cd  # Default to L/D
            
            # Update history
            self.history['iterations'].append(self.eval_count)
            self.history['objective'].append(obj_value)
            self.history['cl'].append(cl)
            self.history['cd'].append(cd)
            self.history['cm'].append(cm)
            self.history['design_vars'].append(design_vars.copy())
            
            # Check if this is the best design so far
            if ((self.objective == 'cd' and obj_value < self.best_objective) or
                (self.objective != 'cd' and -obj_value < self.best_objective)):
                self.best_objective = obj_value if self.objective == 'cd' else -obj_value
                self.best_design = design_vars.copy()
                self.best_airfoil = new_airfoil
                self.best_solution = solution
            
            # Print progress
            self.eval_count += 1
            print(f"Evaluation {self.eval_count}: CL = {cl:.4f}, CD = {cd:.6f}, L/D = {cl/cd:.2f}, Obj = {obj_value:.4f}")
            
            return obj_value
        
        except Exception as e:
            # If analysis fails, return a penalty value
            logger.error(f"Analysis failed: {str(e)}")
            return 1000.0 if self.objective == 'cd' else -0.001
    
    def run_optimization(self, max_iter=20, method='SLSQP'):
        """
        Run the optimization process.
        
        Parameters
        ----------
        max_iter : int
            Maximum number of iterations
        method : str
            Optimization method ('SLSQP', 'BFGS', etc.)
            
        Returns
        -------
        dict
            Optimization results
        """
        print(f"Starting modal design optimization for {self.base_airfoil_name}")
        print(f"Objective: {self.objective}, Mach = {self.mach}, Alpha = {self.alpha}°, Re = {self.reynolds:.2e}")
        print(f"Using {self.n_modes} modal shape functions")
        
        # Initial design variables (all zeros = baseline airfoil)
        x0 = np.zeros(self.n_total_vars)
        
        # Bounds for design variables
        # Thickness modes have tighter bounds to prevent negative thickness
        thickness_bounds = [(-0.02, 0.05) for _ in range(self.n_thickness_modes)]
        camber_bounds = [(-0.05, 0.05) for _ in range(self.n_camber_modes)]
        bounds = thickness_bounds + camber_bounds
        
        # Start timer
        start_time = time.time()
        
        # Run optimization
        result = minimize(
            self.objective_function,
            x0,
            method=method,
            bounds=bounds,
            options={
                'maxiter': max_iter,
                'disp': True
            }
        )
        
        # End timer
        end_time = time.time()
        optimization_time = end_time - start_time
        
        print("\nOptimization completed:")
        print(f"Total time: {optimization_time:.1f} seconds")
        print(f"Function evaluations: {self.eval_count}")
        print(f"Optimization status: {result.message}")
        
        # Get the best design
        if self.best_airfoil is not None:
            # Analyze the best airfoil one more time to get detailed results
            final_solution = self.analyze_airfoil(self.best_airfoil)
            
            # Print final results
            print("\nBest Design Results:")
            print(f"CL: {final_solution['cl']:.4f}")
            print(f"CD: {final_solution['cd']:.6f}")
            print(f"L/D: {final_solution['cl']/final_solution['cd']:.2f}")
            print(f"CM: {final_solution['cm']:.4f}")
            
            # Create result plots
            fig1 = plot_pressure(final_solution, self.best_airfoil)
            fig1.savefig(f'{self.base_airfoil_name}_modal_pressure.png', dpi=300)
            
            fig2 = plot_geometry_comparison(self.base_airfoil, self.best_airfoil)
            fig2.savefig(f'{self.base_airfoil_name}_modal_geometry.png', dpi=300)
            
            # Plot optimization history
            fig3, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
            ax1.plot(self.history['iterations'], self.history['cl'], 'b-', label='CL')
            ax1.set_ylabel('Lift Coefficient')
            ax1.grid(True)
            ax1.legend()
            
            ax2.plot(self.history['iterations'], self.history['cd'], 'r-', label='CD')
            ax2.set_xlabel('Evaluation')
            ax2.set_ylabel('Drag Coefficient')
            ax2.grid(True)
            ax2.legend()
            
            fig3.suptitle('Optimization History')
            fig3.tight_layout()
            fig3.savefig(f'{self.base_airfoil_name}_modal_history.png', dpi=300)
            
            # Save optimized airfoil to file
            self.best_airfoil.save_to_file(f'{self.base_airfoil_name}_modal.dat')
            print(f"Optimized geometry saved to {self.base_airfoil_name}_modal.dat")
            
            return {
                'best_design': self.best_design,
                'best_airfoil': self.best_airfoil,
                'best_solution': final_solution,
                'history': self.history,
                'optimization_time': optimization_time,
                'evaluations': self.eval_count,
                'result': result
            }
        else:
            print("Optimization failed to find a valid design")
            return None

def run_modal_optimization(airfoil_name='naca0012', mach=0.5, alpha=2.0, reynolds=5e6,
                         n_modes=6, objective='cl_cd', max_iter=20):
    """
    Run modal design optimization.
    
    Parameters
    ----------
    airfoil_name : str
        Name of the baseline airfoil
    mach : float
        Freestream Mach number
    alpha : float
        Angle of attack in degrees
    reynolds : float
        Reynolds number based on chord
    n_modes : int
        Number of modal shape functions
    objective : str
        Objective function ('cl', 'cd', 'cl_cd', 'cm')
    max_iter : int
        Maximum number of iterations
        
    Returns
    -------
    dict
        Optimization results
    """
    # Create optimizer
    optimizer = ModalDesignOptimizer(
        base_airfoil=airfoil_name,
        mach=mach,
        alpha=alpha,
        reynolds=reynolds,
        n_modes=n_modes,
        objective=objective
    )
    
    # Run optimization
    results = optimizer.run_optimization(max_iter=max_iter)
    
    return results

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="PyMISES Modal Design Optimization Example")
    parser.add_argument('--airfoil', type=str, default='naca0012', help="Baseline airfoil name")
    parser.add_argument('--mach', type=float, default=0.5, help="Freestream Mach number")
    parser.add_argument('--alpha', type=float, default=2.0, help="Angle of attack (degrees)")
    parser.add_argument('--reynolds', type=float, default=5e6, help="Reynolds number")
    parser.add_argument('--n_modes', type=int, default=6, help="Number of modal shape functions")
    parser.add_argument('--objective', type=str, default='cl_cd', choices=['cl', 'cd', 'cl_cd', 'cm'], 
                      help="Objective function")
    parser.add_argument('--max_iter', type=int, default=20, help="Maximum number of iterations")
    parser.add_argument('--inviscid', action='store_true', help="Run inviscid analysis only")
    
    args = parser.parse_args()
    
    # If inviscid flag is set, set Reynolds to 0
    if args.inviscid:
        args.reynolds = 0
    
    # Run the optimization
    results = run_modal_optimization(
        airfoil_name=args.airfoil,
        mach=args.mach,
        alpha=args.alpha,
        reynolds=args.reynolds,
        n_modes=args.n_modes,
        objective=args.objective,
        max_iter=args.max_iter
    )
    
    # Show plots
    plt.show()
