"""
Streamlit web application for PyMISES.

This module provides a web-based user interface for running PyMISES simulations
using the Streamlit framework.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
import json
import tempfile
import logging
from typing import Dict, Any, Optional, List, Tuple, Union

from pymises.core.geometry import AirfoilGeometry, BladeGeometry
from pymises.core.grid import GridGenerator
from pymises.core.euler import EulerSolver
from pymises.core.boundary_layer import BoundaryLayerFactory
from pymises.core.coupling import CoupledSolver
from pymises.core.newton import NewtonSolver
from pymises.boundary_conditions.wall import InviscidWallBC, ViscousWallBC
from pymises.boundary_conditions.farfield import VortexFarfieldBC
from pymises.postprocessing.visualize import plot_pressure, plot_mach_contours, plot_boundary_layer, plot_convergence_history
from pymises.postprocessing.performance import calculate_airfoil_forces

# Set up module logger
logger = logging.getLogger(__name__)

# Configure logger
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

def build_ui_components() -> None:
    """Build the Streamlit UI components."""
    # Set up page config
    st.set_page_config(
        page_title="PyMISES - Streamtube Euler Solver",
        page_icon="✈️",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    # Add header and description
    st.title("PyMISES - Python Multiple-blade Interacting Streamtube Euler Solver")
    st.markdown("""
    A Python implementation of the MISES solver for aerodynamic analysis and design of airfoils and turbine cascades.
    
    This web application allows you to:
    - Analyze airfoils with coupled Euler/boundary-layer methods
    - Visualize pressure distributions, Mach contours, and boundary layer properties
    - Design airfoils through inverse methods
    """)
    
    # Create sidebar for inputs
    with st.sidebar:
        st.header("Analysis Configuration")
        
        # Analysis mode
        analysis_mode = st.selectbox(
            "Analysis Mode",
            options=["Direct Analysis", "Inverse Design"],
            index=0
        )
        
        # Geometry selection
        st.subheader("Geometry")
        geometry_source = st.radio(
            "Geometry Source",
            options=["Upload File", "NACA Airfoil", "Sample Airfoil"],
            index=2
        )
        
        if geometry_source == "Upload File":
            airfoil_file = st.file_uploader(
                "Upload Airfoil Geometry",
                type=["dat", "txt", "csv"],
                help="Airfoil coordinates in Selig or Lednicer format"
            )
        elif geometry_source == "NACA Airfoil":
            naca_code = st.text_input(
                "NACA Code",
                value="0012",
                help="4-digit or 5-digit NACA airfoil code"
            )
            n_points = st.slider(
                "Number of Points",
                min_value=51,
                max_value=401,
                value=201,
                step=50,
                help="Number of points for airfoil discretization"
            )
        else:  # Sample Airfoil
            sample_airfoil = st.selectbox(
                "Sample Airfoil",
                options=["NACA 0012", "NACA 4412", "RAE 2822", "S809", "SC-1095"],
                index=0
            )
        
        # Flow conditions
        st.subheader("Flow Conditions")
        mach = st.slider(
            "Mach Number",
            min_value=0.0,
            max_value=0.95,
            value=0.3,
            step=0.05,
            help="Freestream Mach number"
        )
        
        reynolds = st.number_input(
            "Reynolds Number",
            min_value=1e4,
            max_value=1e8,
            value=1e6,
            format="%.2e",
            help="Reynolds number based on chord length"
        )
        
        alpha = st.slider(
            "Angle of Attack (degrees)",
            min_value=-10.0,
            max_value=20.0,
            value=0.0,
            step=0.5,
            help="Angle of attack in degrees"
        )
        
        # Advanced options
        st.subheader("Advanced Options")
        show_advanced = st.checkbox("Show Advanced Options", value=False)
        
        if show_advanced:
            with st.expander("Grid Generation"):
                grid_type = st.selectbox(
                    "Grid Type",
                    options=["C-Grid", "O-Grid"],
                    index=0
                )
                
                n_normal = st.slider(
                    "Normal Grid Points",
                    min_value=41,
                    max_value=161,
                    value=81,
                    step=20,
                    help="Number of grid points in the normal direction"
                )
                
                n_wake = st.slider(
                    "Wake Grid Points",
                    min_value=21,
                    max_value=81,
                    value=41,
                    step=10,
                    help="Number of grid points in the wake region"
                )
                
                far_field_radius = st.slider(
                    "Far-Field Radius",
                    min_value=5.0,
                    max_value=50.0,
                    value=20.0,
                    step=5.0,
                    help="Far-field boundary distance from the airfoil"
                )
            
            with st.expander("Solver Options"):
                viscous = st.checkbox(
                    "Include Viscous Effects",
                    value=True,
                    help="Enable boundary layer coupling"
                )
                
                if viscous:
                    transition_model = st.selectbox(
                        "Transition Model",
                        options=["Modified Abu-Ghannam/Shaw", "Envelope e^n Method"],
                        index=0
                    )
                    
                    turbulence_level = st.slider(
                        "Turbulence Level",
                        min_value=0.001,
                        max_value=0.1,
                        value=0.01,
                        format="%.3f",
                        help="Freestream turbulence intensity (0-1)"
                    )
                
                max_iter = st.slider(
                    "Maximum Iterations",
                    min_value=10,
                    max_value=200,
                    value=50,
                    step=10,
                    help="Maximum number of Newton iterations"
                )
                
                tolerance = st.number_input(
                    "Convergence Tolerance",
                    min_value=1e-8,
                    max_value=1e-3,
                    value=1e-6,
                    format="%.1e",
                    help="Residual convergence tolerance"
                )
                
                relaxation = st.slider(
                    "Relaxation Factor",
                    min_value=0.1,
                    max_value=1.0,
                    value=0.7,
                    step=0.1,
                    help="Relaxation factor for Newton iterations"
                )
        
        # Run analysis button
        run_button = st.button(
            "Run Analysis",
            type="primary",
            help="Start the analysis with the current settings"
        )
        
        # If run button is clicked, prepare configuration
        if run_button:
            # Collect configuration from UI inputs
            config = {
                "analysis_mode": analysis_mode,
                "geometry_source": geometry_source.lower().replace(" ", "_"),
                "mach": mach,
                "reynolds": reynolds,
                "alpha": alpha
            }
            
            # Add geometry-specific options
            if geometry_source == "Upload File" and airfoil_file is not None:
                config["airfoil_file"] = airfoil_file
            elif geometry_source == "NACA Airfoil":
                config["naca_code"] = naca_code
                config["n_points"] = n_points
            elif geometry_source == "Sample Airfoil":
                config["sample_airfoil"] = sample_airfoil
            
            # Add advanced options if shown
            if show_advanced:
                if grid_type == "C-Grid":
                    config["grid_type"] = "c-grid"
                else:
                    config["grid_type"] = "o-grid"
                
                config["n_normal"] = n_normal
                config["n_wake"] = n_wake
                config["far_field_radius"] = far_field_radius
                config["viscous"] = viscous
                
                if viscous:
                    if transition_model == "Modified Abu-Ghannam/Shaw":
                        config["transition_model"] = "modified_ags"
                    else:
                        config["transition_model"] = "envelope_en"
                    
                    config["turbulence_level"] = turbulence_level
                
                config["max_iter"] = max_iter
                config["tolerance"] = tolerance
                config["relaxation"] = relaxation
            
            # Show a progress bar during analysis
            progress_bar = st.sidebar.progress(0)
            status_text = st.sidebar.empty()
            
            try:
                # Run the analysis stages
                status_text.text("Setting up geometry and grid...")
                progress_bar.progress(10)
                
                status_text.text("Running Euler solver...")
                progress_bar.progress(30)
                
                status_text.text("Computing boundary layer...")
                progress_bar.progress(60)
                
                status_text.text("Finalizing results...")
                progress_bar.progress(90)
                
                # Run the actual analysis (this would call run_analysis in a real app)
                solution = run_analysis(config)
                
                # Store solution in session state
                st.session_state.solution = solution
                st.session_state.config = config
                
                # Complete progress
                progress_bar.progress(100)
                status_text.text("Analysis complete!")
                
                # Force a rerun to show results
                st.experimental_rerun()
            except Exception as e:
                st.sidebar.error(f"Error during analysis: {str(e)}")
                logger.error(f"Analysis error: {e}", exc_info=True)
            finally:
                # Clean up progress elements
                progress_bar.empty()
                status_text.empty()
    
    # Main content area
    if "solution" not in st.session_state:
        st.info("Configure the analysis parameters in the sidebar and click 'Run Analysis' to begin.")
        
        # Sample images
        col1, col2 = st.columns(2)
        with col1:
            st.markdown("### Example Pressure Distribution")
            st.image("https://raw.githubusercontent.com/your-repository/pymises/main/docs/images/pressure_example.png", 
                    caption="Pressure distribution example", use_column_width=True)
        
        with col2:
            st.markdown("### Example Mach Contours")
            st.image("https://raw.githubusercontent.com/your-repository/pymises/main/docs/images/mach_example.png", 
                    caption="Mach contours example", use_column_width=True)
    else:
        # Display results
        st.header("Analysis Results")
        
        # Create tabs for different results
        tab1, tab2, tab3, tab4 = st.tabs(["Pressure Distribution", "Mach Contours", "Boundary Layer", "Convergence"])
        
        with tab1:
            st.markdown("### Pressure Distribution")
            fig_pressure = plot_pressure(st.session_state.solution)
            st.pyplot(fig_pressure)
        
        with tab2:
            st.markdown("### Mach Contours")
            fig_mach = plot_mach_contours(st.session_state.solution)
            st.pyplot(fig_mach)
        
        with tab3:
            st.markdown("### Boundary Layer Properties")
            fig_bl = plot_boundary_layer(st.session_state.solution)
            st.pyplot(fig_bl)
        
        with tab4:
            st.markdown("### Convergence History")
            fig_conv = plot_convergence_history(st.session_state.solution)
            st.pyplot(fig_conv)
        
        # Display performance metrics
        st.header("Performance Metrics")
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Lift Coefficient (CL)", f"{st.session_state.solution.get('cl', 0.0):.4f}")
        
        with col2:
            st.metric("Drag Coefficient (CD)", f"{st.session_state.solution.get('cd', 0.0):.6f}")
        
        with col3:
            st.metric("Moment Coefficient (CM)", f"{st.session_state.solution.get('cm', 0.0):.4f}")
        
        # Display transition and separation information if available
        if "transition_location" in st.session_state.solution:
            st.subheader("Boundary Layer Transition")
            trans_loc = st.session_state.solution["transition_location"]
            
            col1, col2 = st.columns(2)
            with col1:
                if "upper" in trans_loc and trans_loc["upper"] is not None:
                    st.metric("Upper Surface (x/c)", f"{trans_loc['upper']:.4f}")
                else:
                    st.metric("Upper Surface (x/c)", "N/A")
            
            with col2:
                if "lower" in trans_loc and trans_loc["lower"] is not None:
                    st.metric("Lower Surface (x/c)", f"{trans_loc['lower']:.4f}")
                else:
                    st.metric("Lower Surface (x/c)", "N/A")
        
        if "separation_location" in st.session_state.solution:
            st.subheader("Boundary Layer Separation")
            sep_loc = st.session_state.solution["separation_location"]
            
            col1, col2 = st.columns(2)
            with col1:
                if "upper" in sep_loc and sep_loc["upper"] is not None:
                    st.metric("Upper Surface (x/c)", f"{sep_loc['upper']:.4f}")
                else:
                    st.metric("Upper Surface (x/c)", "None")
            
            with col2:
                if "lower" in sep_loc and sep_loc["lower"] is not None:
                    st.metric("Lower Surface (x/c)", f"{sep_loc['lower']:.4f}")
                else:
                    st.metric("Lower Surface (x/c)", "None")
        
        # Add download buttons for results
        st.header("Download Results")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            if st.button("Download Pressure Distribution (CSV)"):
                # This would create a downloadable CSV in a real app
                pass
        
        with col2:
            if st.button("Download Solution (VTK)"):
                # This would create a downloadable VTK file in a real app
                pass
        
        with col3:
            if st.button("Download Performance Report (TXT)"):
                # This would create a downloadable report in a real app
                pass

def run_analysis(config: Dict[str, Any]) -> Dict[str, Any]:
    """Run PyMISES analysis with the given configuration.
    
    Args:
        config: Analysis configuration dictionary
        
    Returns:
        Solution dictionary
    """
    # This is a simplified version of the analysis code
    # In practice, this would use the same code as the CLI runner
    
    # Load or create geometry
    if "geometry_source" in config:
        if config["geometry_source"] == "upload_file" and "airfoil_file" in config:
            # Save uploaded file to temporary file
            with tempfile.NamedTemporaryFile(suffix=".dat", delete=False) as tmp:
                tmp.write(config["airfoil_file"].read())
                tmp_path = tmp.name
            
            try:
                airfoil = AirfoilGeometry.load_from_file(tmp_path)
            finally:
                os.unlink(tmp_path)  # Delete temporary file
        
        elif config["geometry_source"] == "naca_airfoil" and "naca_code" in config:
            airfoil = AirfoilGeometry.create_naca(
                config["naca_code"],
                n_points=config.get("n_points", 201)
            )
        
        elif config["geometry_source"] == "sample_airfoil" and "sample_airfoil" in config:
            # Determine which sample airfoil to use
            sample_name = config["sample_airfoil"]
            if sample_name == "NACA 0012":
                airfoil = AirfoilGeometry.create_naca("0012")
            elif sample_name == "NACA 4412":
                airfoil = AirfoilGeometry.create_naca("4412")
            elif sample_name == "RAE 2822":
                # Path to sample airfoil data (assuming it's included with PyMISES)
                sample_path = os.path.join(
                    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "examples", "validation", "data", "rae2822.dat"
                )
                if os.path.exists(sample_path):
                    airfoil = AirfoilGeometry.load_from_file(sample_path)
                else:
                    st.error(f"Sample airfoil data not found: {sample_path}")
                    airfoil = AirfoilGeometry.create_naca("0012")  # Fallback
            else:
                # Fallback to NACA 0012 for other samples not yet implemented
                airfoil = AirfoilGeometry.create_naca("0012")
        else:
            # Default to NACA 0012
            airfoil = AirfoilGeometry.create_naca("0012")
    else:
        # Default to NACA 0012
        airfoil = AirfoilGeometry.create_naca("0012")
    
    # Generate grid
    grid_type = config.get("grid_type", "c-grid").lower()
    n_normal = config.get("n_normal", 81)
    n_wake = config.get("n_wake", 41)
    far_field_radius = config.get("far_field_radius", 20.0)
    
    grid_gen = GridGenerator(airfoil, grid_type=grid_type)
    grid = grid_gen.generate_grid(
        n_normal=n_normal,
        n_wake=n_wake,
        far_field_radius=far_field_radius
    )
    
    # Set up flow conditions
    mach = config.get("mach", 0.3)
    reynolds = config.get("reynolds", 1e6)
    alpha = config.get("alpha", 0.0)
    alpha_rad = alpha * np.pi / 180.0  # Convert to radians
    
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
    max_iter = config.get("max_iter", 50)
    tolerance = config.get("tolerance", 1e-6)
    relaxation = config.get("relaxation", 0.7)
    
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

def launch_app():
    """Launch the Streamlit application."""
    build_ui_components()

if __name__ == "__main__":
    launch_app()