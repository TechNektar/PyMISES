"""
Command-line interface for PyMISES.

This module provides a command-line interface for running PyMISES simulations.
It handles argument parsing, configuration loading, and execution of analysis
and design cases.
"""

import argparse
import os
import sys
import json
import yaml
import logging
from typing import Dict, Any, Optional, List, Tuple

from pymises.core.geometry import AirfoilGeometry, BladeGeometry
from pymises.core.grid import GridGenerator
from pymises.core.euler import EulerSolver
from pymises.core.boundary_layer import BoundaryLayerFactory
from pymises.core.coupling import CoupledSolver
from pymises.core.newton import NewtonSolver
from pymises.boundary_conditions.wall import InviscidWallBC, ViscousWallBC
from pymises.boundary_conditions.farfield import VortexFarfieldBC
from pymises.postprocessing.visualize import plot_pressure, plot_mach_contours, plot_convergence_history
from pymises.postprocessing.export import export_solution_to_vtk, export_performance_report

# Set up module logger
logger = logging.getLogger(__name__)

def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments for PyMISES.
    
    Returns:
        Parsed command-line arguments
    """
    parser = argparse.ArgumentParser(
        description="PyMISES - Python Multiple-blade Interacting Streamtube Euler Solver",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Add command-line arguments
    parser.add_argument(
        "--config", "-c",
        type=str,
        help="Path to configuration file (JSON or YAML)",
        default=None
    )
    
    parser.add_argument(
        "--mode", "-m",
        type=str,
        choices=["direct", "inverse"],
        help="Analysis mode (direct or inverse design)",
        default="direct"
    )
    
    parser.add_argument(
        "--geometry", "-g",
        type=str,
        help="Path to geometry file (airfoil coordinates)",
        default=None
    )
    
    parser.add_argument(
        "--mach", 
        type=float,
        help="Freestream Mach number",
        default=0.3
    )
    
    parser.add_argument(
        "--reynolds",
        type=float,
        help="Reynolds number",
        default=1e6
    )
    
    parser.add_argument(
        "--alpha",
        type=float,
        help="Angle of attack in degrees",
        default=0.0
    )
    
    parser.add_argument(
        "--output-dir", "-o",
        type=str,
        help="Directory for output files",
        default="./results"
    )
    
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Generate plots"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="count",
        default=0,
        help="Increase verbosity (can be used multiple times)"
    )
    
    return parser.parse_args()

def load_configuration(args: argparse.Namespace) -> Dict[str, Any]:
    """
    Load configuration from file and/or command-line arguments.
    
    Args:
        args: Parsed command-line arguments
        
    Returns:
        Configuration dictionary
    """
    config = {
        "mode": args.mode,
        "flow_conditions": {
            "mach": args.mach,
            "reynolds": args.reynolds,
            "alpha": args.alpha
        },
        "output": {
            "directory": args.output_dir,
            "generate_plots": args.plot
        }
    }
    
    # If config file is specified, load it
    if args.config is not None:
        if not os.path.exists(args.config):
            logger.error(f"Configuration file not found: {args.config}")
            sys.exit(1)
        
        # Load configuration from file
        try:
            with open(args.config, "r") as f:
                if args.config.endswith(".json"):
                    file_config = json.load(f)
                elif args.config.endswith((".yaml", ".yml")):
                    import yaml
                    file_config = yaml.safe_load(f)
                else:
                    logger.error(f"Unsupported configuration file format: {args.config}")
                    sys.exit(1)
            
            # Merge configurations (file has precedence over command-line args)
            _deep_update(config, file_config)
        except Exception as e:
            logger.error(f"Error loading configuration file: {e}")
            sys.exit(1)
    
    # Override with explicit command-line args
    if args.geometry is not None:
        config["geometry"] = {"file": args.geometry}
    
    return config

def _deep_update(d: Dict[str, Any], u: Dict[str, Any]) -> Dict[str, Any]:
    """
    Recursively update a dictionary with another dictionary.
    
    Args:
        d: Dictionary to update
        u: Dictionary with updates
        
    Returns:
        Updated dictionary
    """
    for k, v in u.items():
        if isinstance(v, dict) and k in d and isinstance(d[k], dict):
            d[k] = _deep_update(d[k], v)
        else:
            d[k] = v
    return d

def setup_logging(verbosity: int) -> None:
    """
    Set up logging with appropriate verbosity level.
    
    Args:
        verbosity: Verbosity level (0=WARNING, 1=INFO, 2=DEBUG)
    """
    log_levels = {
        0: logging.WARNING,
        1: logging.INFO,
        2: logging.DEBUG
    }
    
    # Set log level (default to DEBUG for verbosity > 2)
    log_level = log_levels.get(verbosity, logging.DEBUG)
    
    # Configure logging
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

def load_geometry(config: Dict[str, Any]) -> AirfoilGeometry:
    """
    Load airfoil geometry from configuration.
    
    Args:
        config: Configuration dictionary
        
    Returns:
        Airfoil geometry object
    """
    geometry_config = config.get("geometry", {})
    
    if "file" in geometry_config:
        # Load from file
        file_path = geometry_config["file"]
        if not os.path.exists(file_path):
            logger.error(f"Geometry file not found: {file_path}")
            sys.exit(1)
        
        try:
            airfoil = AirfoilGeometry.load_from_file(file_path)
            logger.info(f"Loaded geometry from file: {file_path}")
            return airfoil
        except Exception as e:
            logger.error(f"Error loading geometry file: {e}")
            sys.exit(1)
    elif "naca" in geometry_config:
        # Create NACA airfoil
        naca_code = geometry_config["naca"]
        n_points = geometry_config.get("n_points", 201)
        
        airfoil = AirfoilGeometry.create_naca(naca_code, n_points=n_points)
        logger.info(f"Created NACA {naca_code} airfoil with {n_points} points")
        return airfoil
    else:
        # Default to NACA 0012
        logger.warning("No geometry specified, using default NACA 0012")
        return AirfoilGeometry.create_naca("0012")

def run_direct_analysis(config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Run direct analysis simulation.
    
    Args:
        config: Configuration dictionary
        
    Returns:
        Solution dictionary
    """
    # Load geometry
    airfoil = load_geometry(config)
    
    # Extract flow conditions
    flow_config = config.get("flow_conditions", {})
    mach = flow_config.get("mach", 0.3)
    reynolds = flow_config.get("reynolds", 1e6)
    alpha = flow_config.get("alpha", 0.0)
    alpha_rad = alpha * 3.14159 / 180.0  # Convert to radians
    
    # Generate grid
    grid_config = config.get("grid", {})
    grid_type = grid_config.get("type", "c-grid")
    n_normal = grid_config.get("n_normal", 81)
    n_wake = grid_config.get("n_wake", 41)
    far_field_radius = grid_config.get("far_field_radius", 20.0)
    
    logger.info(f"Generating {grid_type} grid with {n_normal}x{n_wake} points")
    grid_gen = GridGenerator(airfoil, grid_type=grid_type)
    grid = grid_gen.generate_grid(
        n_normal=n_normal,
        n_wake=n_wake,
        far_field_radius=far_field_radius
    )
    
    # Set up boundary conditions
    airfoil_indices = grid.get_boundary_indices("airfoil")
    wall_bc = InviscidWallBC(airfoil_indices, normal_direction="inner")
    
    farfield_indices = grid.get_boundary_indices("farfield")
    farfield_bc = VortexFarfieldBC(
        farfield_indices,
        mach,
        alpha_rad,
        p_inf=101325.0,  # Standard atmosphere
        T_inf=288.15,    # Standard atmosphere
        circulation=0.0, # Will be updated after initial solution
        airfoil_x=0.25,  # Quarter-chord position
        airfoil_y=0.0
    )
    
    # Set up Euler solver
    euler_solver = EulerSolver(grid)
    euler_solver.add_boundary_condition(wall_bc)
    euler_solver.add_boundary_condition(farfield_bc)
    
    # Initialize flow field
    euler_solver.initialize(
        mach=mach,
        alpha=alpha_rad,
        p0=101325.0 * (1 + 0.2*mach**2)**3.5,  # Total pressure
        T0=288.15 * (1 + 0.2*mach**2)          # Total temperature
    )
    
    # Set up Newton solver
    newton_config = config.get("newton", {})
    max_iter = newton_config.get("max_iter", 50)
    tolerance = newton_config.get("tolerance", 1e-6)
    relaxation = newton_config.get("relaxation", 0.7)
    
    logger.info(f"Starting inviscid solution with max_iter={max_iter}, tolerance={tolerance}")
    newton_solver = NewtonSolver(
        residual_function=euler_solver.compute_residuals,
        jacobian_function=euler_solver.compute_jacobian,
        solution=euler_solver.get_solution()
    )
    
    # Run the Newton solver for inviscid solution
    inviscid_solution, euler_convergence = newton_solver.solve(
        max_iter=max_iter,
        tolerance=tolerance,
        relaxation=relaxation
    )
    
    # Update the Euler solver with the converged solution
    euler_solver.set_solution(inviscid_solution)
    
    # Calculate lift and update circulation
    inviscid_forces = euler_solver.compute_forces()
    lift = inviscid_forces["cl"]
    circulation = lift * 2.0  # Simple relationship based on Kutta-Joukowski
    
    # Update far-field boundary condition with calculated circulation
    farfield_bc.circulation = circulation
    
    # Run another iteration of the Euler solver with the updated circulation
    inviscid_solution, _ = newton_solver.solve(
        max_iter=20,
        tolerance=tolerance,
        relaxation=relaxation
    )
    
    # Add viscous effects if requested
    if config.get("viscous", True):
        logger.info("Adding viscous effects with boundary layer coupling")
        
        # Create boundary layer solver factory
        bl_factory = BoundaryLayerFactory(
            reynolds,
            transition_model=config.get("transition_model", "modified_ags")
        )
        
        # Create viscous wall boundary condition
        def displacement_thickness_provider(idx):
            # This will be updated by the coupled solver
            return 0.0
        
        viscous_wall_bc = ViscousWallBC(
            airfoil_indices,
            normal_direction="inner",
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
        logger.info("Running viscous-inviscid interaction")
        viscous_solution, viscous_convergence = coupled_newton.solve(
            max_iter=50,
            tolerance=1e-5,
            relaxation=0.6
        )
        
        # Extract boundary layer properties
        bl_properties = coupled_solver.get_boundary_layer_properties()
        final_solution = viscous_solution.copy()
        final_solution.update(bl_properties)
    else:
        # Use inviscid solution
        final_solution = inviscid_solution
        viscous_convergence = []
    
    # Compute final aerodynamic coefficients
    forces = euler_solver.compute_forces()
    final_solution.update({
        "cl": forces["cl"],
        "cd": forces["cd"],
        "cm": forces["cm"],
        "euler_convergence": euler_convergence,
        "viscous_convergence": viscous_convergence,
        "config": config  # Store the configuration for reference
    })
    
    return final_solution

def save_results(solution: Dict[str, Any], config: Dict[str, Any]) -> None:
    """
    Save simulation results to disk.
    
    Args:
        solution: Solution dictionary
        config: Configuration dictionary
    """
    output_config = config.get("output", {})
    output_dir = output_config.get("directory", "./results")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Export solution to VTK
    vtk_file = os.path.join(output_dir, "solution.vtk")
    export_solution_to_vtk(solution, vtk_file)
    logger.info(f"Exported solution to VTK: {vtk_file}")
    
    # Export performance report
    report_file = os.path.join(output_dir, "performance.txt")
    export_performance_report(solution, report_file)
    logger.info(f"Exported performance report: {report_file}")
    
    # Save convergence history
    convergence_file = os.path.join(output_dir, "convergence.csv")
    with open(convergence_file, "w") as f:
        f.write("iteration,euler_residual,viscous_residual\n")
        for i, (euler_res, viscous_res) in enumerate(
            zip(
                solution.get("euler_convergence", []),
                solution.get("viscous_convergence", []) + [float("nan")] * 
                (len(solution.get("euler_convergence", [])) - len(solution.get("viscous_convergence", [])))
            )
        ):
            f.write(f"{i+1},{euler_res},{viscous_res if not np.isnan(viscous_res) else ''}\n")
    logger.info(f"Saved convergence history: {convergence_file}")
    
    # Generate plots if requested
    if output_config.get("generate_plots", False):
        # Pressure distribution
        pressure_file = os.path.join(output_dir, "pressure.png")
        fig = plot_pressure(solution)
        fig.savefig(pressure_file, dpi=300)
        logger.info(f"Generated pressure plot: {pressure_file}")
        
        # Mach contours
        mach_file = os.path.join(output_dir, "mach.png")
        fig = plot_mach_contours(solution)
        fig.savefig(mach_file, dpi=300)
        logger.info(f"Generated Mach contour plot: {mach_file}")
        
        # Convergence history
        convergence_file = os.path.join(output_dir, "convergence.png")
        fig = plot_convergence_history(solution)
        fig.savefig(convergence_file, dpi=300)
        logger.info(f"Generated convergence plot: {convergence_file}")

def run_cli() -> None:
    """
    Main entry point for the command-line interface.
    """
    # Parse command-line arguments
    args = parse_arguments()
    
    # Set up logging
    setup_logging(args.verbose)
    
    # Load configuration
    config = load_configuration(args)
    
    # Run simulation based on mode
    if config["mode"] == "direct":
        logger.info("Running direct analysis...")
        solution = run_direct_analysis(config)
    elif config["mode"] == "inverse":
        logger.error("Inverse design mode not implemented yet")
        sys.exit(1)
    else:
        logger.error(f"Unsupported mode: {config['mode']}")
        sys.exit(1)
    
    # Save results
    save_results(solution, config)
    
    # Print summary
    print("\nSolution Summary:")
    print(f"  CL = {solution.get('cl', 0.0):.4f}")
    print(f"  CD = {solution.get('cd', 0.0):.6f}")
    print(f"  CM = {solution.get('cm', 0.0):.4f}")
    
    # If solution has transition and separation information, print it
    if "transition_location" in solution:
        trans_loc = solution["transition_location"]
        print("\nTransition Location (x/c):")
        if "upper" in trans_loc and trans_loc["upper"] is not None:
            print(f"  Upper Surface: {trans_loc['upper']:.4f}")
        else:
            print("  Upper Surface: N/A")
            
        if "lower" in trans_loc and trans_loc["lower"] is not None:
            print(f"  Lower Surface: {trans_loc['lower']:.4f}")
        else:
            print("  Lower Surface: N/A")
    
    if "separation_location" in solution:
        sep_loc = solution["separation_location"]
        print("\nSeparation Location (x/c):")
        if "upper" in sep_loc and sep_loc["upper"] is not None:
            print(f"  Upper Surface: {sep_loc['upper']:.4f}")
        else:
            print("  Upper Surface: None")
            
        if "lower" in sep_loc and sep_loc["lower"] is not None:
            print(f"  Lower Surface: {sep_loc['lower']:.4f}")
        else:
            print("  Lower Surface: None")
    
    logger.info("Analysis completed successfully")

if __name__ == "__main__":
    try:
        import numpy as np  # Required for save_results function
        run_cli()
    except KeyboardInterrupt:
        logger.info("Execution interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Unhandled exception: {e}", exc_info=True)
        sys.exit(1)
