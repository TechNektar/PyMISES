"""
PyMISES - Visualization Module

This module provides functions for visualizing simulation results, including
pressure distributions, streamlines, Mach contours, and boundary layer properties.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.ticker as ticker
from typing import Dict, List, Tuple, Optional, Union

from pymises.core.geometry import AirfoilGeometry

def plot_pressure(solution: Dict, geometry, fig=None, ax=None):
    """
    Plot pressure distribution around an airfoil or blade.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing pressure and grid data
    geometry : AirfoilGeometry or BladeGeometry
        Geometry object containing airfoil or blade information
    fig : matplotlib.figure.Figure, optional
        Figure to plot on, if None, a new figure is created
    ax : matplotlib.axes.Axes, optional
        Axes to plot on, if None, a new axes is created
        
    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the plot
    """
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    
    # Extract data
    pressure = solution['pressure']
    grid_x = solution['grid_x']
    grid_y = solution['grid_y']
    
    # Get airfoil surface points
    surface_indices = geometry.get_surface_indices()
    
    # Compute pressure coefficient
    # Cp = (p - p_inf) / (0.5 * rho_inf * V_inf^2)
    if 'mach' in solution and 'p_inf' in solution:
        mach = solution['mach']
        p_inf = solution['p_inf']
        cp = (pressure[surface_indices] - p_inf) / (0.5 * p_inf * 1.4 * mach**2)
    else:
        # Normalize pressure for visualization if no reference values available
        p_min = np.min(pressure[surface_indices])
        p_max = np.max(pressure[surface_indices])
        cp = (pressure[surface_indices] - p_min) / (p_max - p_min)
    
    # Get surface coordinates
    x = grid_x[surface_indices]
    y = grid_y[surface_indices]
    
    # Sort points by x-coordinate for airfoil or cascade
    if isinstance(geometry, AirfoilGeometry):
        # For airfoil, separate upper and lower surfaces
        # Assuming leading edge is at min(x) and trailing edge at max(x)
        le_idx = np.argmin(x)
        upper_idx = np.arange(le_idx, len(x))
        lower_idx = np.arange(0, le_idx)
        
        # Plot upper surface
        ax.plot(x[upper_idx], -cp[upper_idx], 'b-', label='Upper Surface')
        
        # Plot lower surface
        ax.plot(x[lower_idx], -cp[lower_idx], 'r-', label='Lower Surface')
    else:
        # For cascade, plot along arc length
        arc_length = np.zeros_like(x)
        arc_length[1:] = np.sqrt(np.diff(x)**2 + np.diff(y)**2)
        arc_length = np.cumsum(arc_length)
        arc_length = arc_length / arc_length[-1]  # Normalize to [0, 1]
        
        # Plot pressure vs. arc length
        ax.plot(arc_length, -cp, 'b-', label='Pressure')
    
    # Set up plot
    ax.set_xlabel('x/c or s/s_max')
    ax.set_ylabel('-Cp')
    ax.set_title('Pressure Distribution')
    ax.grid(True)
    ax.legend()
    
    return fig

def plot_grid(grid, solution=None, fig=None, ax=None):
    """
    Plot computational grid with optional solution overlay.
    
    Parameters
    ----------
    grid : Grid
        Computational grid object
    solution : Dict, optional
        Solution dictionary containing flow data
    fig : matplotlib.figure.Figure, optional
        Figure to plot on, if None, a new figure is created
    ax : matplotlib.axes.Axes, optional
        Axes to plot on, if None, a new axes is created
        
    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the plot
    """
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))
    
    # Extract grid coordinates
    x = grid.x
    y = grid.y
    ni = grid.ni
    nj = grid.nj
    
    # Plot grid lines
    for i in range(ni):
        ax.plot(x[i, :], y[i, :], 'k-', linewidth=0.5, alpha=0.3)
    
    for j in range(nj):
        ax.plot(x[:, j], y[:, j], 'k-', linewidth=0.5, alpha=0.3)
    
    # Overlay solution if provided
    if solution is not None and 'mach' in solution:
        # Create structured grid for contour plot
        x_grid = x.reshape(ni, nj)
        y_grid = y.reshape(ni, nj)
        mach = solution['mach'].reshape(ni, nj)
        
        # Plot Mach contours
        contour = ax.contourf(x_grid, y_grid, mach, levels=20, cmap=cm.jet, alpha=0.6)
        fig.colorbar(contour, ax=ax, label='Mach Number')
    
    # Set up plot
    ax.set_aspect('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Computational Grid')
    
    return fig

def plot_mach_contours(solution, grid, fig=None, ax=None):
    """
    Plot Mach number contours.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing Mach number data
    grid : Grid
        Computational grid object
    fig : matplotlib.figure.Figure, optional
        Figure to plot on, if None, a new figure is created
    ax : matplotlib.axes.Axes, optional
        Axes to plot on, if None, a new axes is created
        
    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the plot
    """
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))
    
    # Extract grid coordinates
    x = grid.x
    y = grid.y
    ni = grid.ni
    nj = grid.nj
    
    # Check if solution contains velocity or Mach number
    if 'mach' in solution:
        mach = solution['mach']
    elif 'velocity_x' in solution and 'velocity_y' in solution and 'speed_of_sound' in solution:
        vx = solution['velocity_x']
        vy = solution['velocity_y']
        a = solution['speed_of_sound']
        mach = np.sqrt(vx**2 + vy**2) / a
    else:
        raise ValueError("Solution must contain Mach number or velocity and speed of sound")
    
    # Reshape for contour plot
    x_grid = x.reshape(ni, nj)
    y_grid = y.reshape(ni, nj)
    mach_grid = mach.reshape(ni, nj)
    
    # Create contour plot
    contour = ax.contourf(x_grid, y_grid, mach_grid, levels=20, cmap=cm.jet, alpha=0.9)
    fig.colorbar(contour, ax=ax, label='Mach Number')
    
    # Add some contour lines
    contour_lines = ax.contour(x_grid, y_grid, mach_grid, colors='black', linewidths=0.5,
                              levels=[0.3, 0.5, 0.7, 0.9, 1.0, 1.1, 1.3])
    ax.clabel(contour_lines, inline=True, fontsize=8, fmt='%.1f')
    
    # Overlay grid outline
    ax.plot(x_grid[:, 0], y_grid[:, 0], 'k-', linewidth=1)
    ax.plot(x_grid[:, -1], y_grid[:, -1], 'k-', linewidth=1)
    ax.plot(x_grid[0, :], y_grid[0, :], 'k-', linewidth=1)
    ax.plot(x_grid[-1, :], y_grid[-1, :], 'k-', linewidth=1)
    
    # Set up plot
    ax.set_aspect('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Mach Number Contours')
    
    return fig

def plot_streamlines(solution, grid, density=1.0, fig=None, ax=None):
    """
    Plot streamlines of the flow.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing velocity data
    grid : Grid
        Computational grid object
    density : float, optional
        Density of streamlines, higher values mean more streamlines
    fig : matplotlib.figure.Figure, optional
        Figure to plot on, if None, a new figure is created
    ax : matplotlib.axes.Axes, optional
        Axes to plot on, if None, a new axes is created
        
    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the plot
    """
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))
    
    # Extract grid coordinates
    x = grid.x
    y = grid.y
    ni = grid.ni
    nj = grid.nj
    
    # Extract velocity components
    vx = solution['velocity_x']
    vy = solution['velocity_y']
    
    # Reshape for streamplot
    x_grid = x.reshape(ni, nj)
    y_grid = y.reshape(ni, nj)
    vx_grid = vx.reshape(ni, nj)
    vy_grid = vy.reshape(ni, nj)
    
    # Create streamplot
    streamplot = ax.streamplot(x_grid.T, y_grid.T, vx_grid.T, vy_grid.T,
                              density=density, color='blue', linewidth=0.8, arrowsize=1.2)
    
    # Overlay grid outline
    ax.plot(x_grid[:, 0], y_grid[:, 0], 'k-', linewidth=1)
    ax.plot(x_grid[:, -1], y_grid[:, -1], 'k-', linewidth=1)
    ax.plot(x_grid[0, :], y_grid[0, :], 'k-', linewidth=1)
    ax.plot(x_grid[-1, :], y_grid[-1, :], 'k-', linewidth=1)
    
    # Set up plot
    ax.set_aspect('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Flow Streamlines')
    
    return fig

def plot_boundary_layer(solution, geometry, fig=None, ax=None):
    """
    Plot boundary layer properties along the surface.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing boundary layer data
    geometry : AirfoilGeometry or BladeGeometry
        Geometry object containing airfoil or blade information
    fig : matplotlib.figure.Figure, optional
        Figure to plot on, if None, a new figure is created
    ax : matplotlib.axes.Axes, optional
        Axes to plot on, if None, a new axes is created
        
    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the plot
    """
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    
    # Check if boundary layer data exists
    if 'displacement_thickness' not in solution or 'momentum_thickness' not in solution:
        raise ValueError("Solution missing boundary layer data")
    
    # Extract boundary layer parameters
    delta_star = solution['displacement_thickness']
    theta = solution['momentum_thickness']
    
    # Compute shape factor H
    H = delta_star / theta
    
    # Get surface coordinates
    surface_indices = geometry.get_surface_indices()
    x = solution['grid_x'][surface_indices]
    
    # Plot shape factor
    ax.plot(x, H, 'b-', label='Shape Factor (H)')
    
    # Create twin axis for displacement and momentum thickness
    ax2 = ax.twinx()
    ax2.plot(x, delta_star, 'r-', label='δ*')
    ax2.plot(x, theta, 'g-', label='θ')
    
    # Highlight transition and separation points if available
    if 'transition_location' in solution:
        trans_loc = solution['transition_location']
        if isinstance(trans_loc, dict):
            # For airfoil with upper/lower surfaces
            if 'upper' in trans_loc:
                ax.axvline(x=trans_loc['upper'], color='purple', linestyle='--')
                ax.text(trans_loc['upper'], ax.get_ylim()[1]*0.9, 'Trans. (U)', 
                       rotation=90, verticalalignment='top')
            if 'lower' in trans_loc:
                ax.axvline(x=trans_loc['lower'], color='purple', linestyle='--')
                ax.text(trans_loc['lower'], ax.get_ylim()[1]*0.9, 'Trans. (L)', 
                       rotation=90, verticalalignment='top')
        else:
            # For cascade with single surface
            ax.axvline(x=trans_loc, color='purple', linestyle='--')
            ax.text(trans_loc, ax.get_ylim()[1]*0.9, 'Transition', 
                   rotation=90, verticalalignment='top')
    
    if 'separation_location' in solution:
        sep_loc = solution['separation_location']
        if isinstance(sep_loc, dict):
            # For airfoil with upper/lower surfaces
            if 'upper' in sep_loc and sep_loc['upper'] > 0:
                ax.axvline(x=sep_loc['upper'], color='red', linestyle='--')
                ax.text(sep_loc['upper'], ax.get_ylim()[1]*0.9, 'Sep. (U)', 
                       rotation=90, verticalalignment='top')
            if 'lower' in sep_loc and sep_loc['lower'] > 0:
                ax.axvline(x=sep_loc['lower'], color='red', linestyle='--')
                ax.text(sep_loc['lower'], ax.get_ylim()[1]*0.9, 'Sep. (L)', 
                       rotation=90, verticalalignment='top')
        else:
            # For cascade with single surface
            if sep_loc > 0:
                ax.axvline(x=sep_loc, color='red', linestyle='--')
                ax.text(sep_loc, ax.get_ylim()[1]*0.9, 'Separation', 
                       rotation=90, verticalalignment='top')
    
    # Set up plot
    ax.set_xlabel('x/c')
    ax.set_ylabel('Shape Factor (H)')
    ax2.set_ylabel('Thickness (δ*, θ)')
    ax.set_title('Boundary Layer Properties')
    ax.grid(True)
    
    # Combine legends
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
    
    return fig

def plot_geometry_comparison(original_geometry, modified_geometry, fig=None, ax=None):
    """
    Plot comparison between original and modified geometries.
    
    Parameters
    ----------
    original_geometry : AirfoilGeometry or BladeGeometry
        Original geometry object
    modified_geometry : AirfoilGeometry or BladeGeometry
        Modified geometry object
    fig : matplotlib.figure.Figure, optional
        Figure to plot on, if None, a new figure is created
    ax : matplotlib.axes.Axes, optional
        Axes to plot on, if None, a new axes is created
        
    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the plot
    """
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    
    # Extract coordinates
    orig_x = original_geometry.x
    orig_y = original_geometry.y
    mod_x = modified_geometry.x
    mod_y = modified_geometry.y
    
    # Plot geometries
    ax.plot(orig_x, orig_y, 'b-', label='Original', linewidth=1.5)
    ax.plot(mod_x, mod_y, 'r--', label='Modified', linewidth=1.5)
    
    # Calculate and plot difference
    if len(orig_x) == len(mod_x) and np.allclose(orig_x, mod_x):
        # If x-coordinates match, plot direct difference
        ax.plot(orig_x, (mod_y - orig_y) * 10, 'g-', label='Difference (×10)', linewidth=1.0)
    
    # Set up plot
    ax.set_aspect('equal')
    ax.set_xlabel('x/c')
    ax.set_ylabel('y/c')
    ax.set_title('Geometry Comparison')
    ax.grid(True)
    ax.legend()
    
    return fig

def plot_convergence_history(convergence_history, fig=None, ax=None):
    """
    Plot convergence history of a simulation.
    
    Parameters
    ----------
    convergence_history : Dict or array-like
        Convergence history data
    fig : matplotlib.figure.Figure, optional
        Figure to plot on, if None, a new figure is created
    ax : matplotlib.axes.Axes, optional
        Axes to plot on, if None, a new axes is created
        
    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the plot
    """
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    
    # Extract residuals
    if isinstance(convergence_history, dict):
        iterations = convergence_history.get('iterations', np.arange(len(next(iter(convergence_history.values())))))
        
        for key, values in convergence_history.items():
            if key != 'iterations':
                ax.semilogy(iterations, values, label=key)
    else:
        # If it's a simple array or list, just plot it
        ax.semilogy(np.arange(len(convergence_history)), convergence_history, label='Residual')
    
    # Set up plot
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Residual')
    ax.set_title('Convergence History')
    ax.grid(True)
    if len(ax.get_legend_handles_labels()[0]) > 1:
        ax.legend()
    
    return fig

def create_performance_report(solution, output_file=None):
    """
    Create a performance report for the simulation.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing performance data
    output_file : str, optional
        File to save the report to, if None, report is printed to console
        
    Returns
    -------
    str
        Report text
    """
    # Generate report text
    report = "Performance Report\n"
    report += "=================\n\n"
    
    # Flow conditions
    report += "Flow Conditions\n"
    report += "--------------\n"
    if 'mach' in solution:
        report += f"Mach Number: {solution['mach']:.3f}\n"
    if 'reynolds' in solution:
        report += f"Reynolds Number: {solution['reynolds']:.2e}\n"
    if 'alpha' in solution:
        report += f"Angle of Attack: {np.degrees(solution['alpha']):.2f}°\n"
    report += "\n"
    
    # Aerodynamic coefficients
    report += "Aerodynamic Coefficients\n"
    report += "-----------------------\n"
    if 'cl' in solution:
        report += f"Lift Coefficient (CL): {solution['cl']:.4f}\n"
    if 'cd' in solution:
        report += f"Drag Coefficient (CD): {solution['cd']:.6f}\n"
    if 'cm' in solution:
        report += f"Moment Coefficient (CM): {solution['cm']:.4f}\n"
    if 'cl' in solution and 'cd' in solution:
        report += f"Lift-to-Drag Ratio (L/D): {solution['cl']/solution['cd']:.2f}\n"
    report += "\n"
    
    # Boundary layer properties
    report += "Boundary Layer Properties\n"
    report += "------------------------\n"
    if 'transition_location' in solution:
        trans_loc = solution['transition_location']
        if isinstance(trans_loc, dict):
            # For airfoil with upper/lower surfaces
            report += f"Transition Location (Upper): {trans_loc.get('upper', 'N/A'):.4f} x/c\n"
            report += f"Transition Location (Lower): {trans_loc.get('lower', 'N/A'):.4f} x/c\n"
        else:
            # For cascade with single surface
            report += f"Transition Location: {trans_loc:.4f} x/c\n"
    
    if 'separation_location' in solution:
        sep_loc = solution['separation_location']
        if isinstance(sep_loc, dict):
            # For airfoil with upper/lower surfaces
            upper_sep = sep_loc.get('upper', -1)
            lower_sep = sep_loc.get('lower', -1)
            if upper_sep > 0:
                report += f"Separation Location (Upper): {upper_sep:.4f} x/c\n"
            else:
                report += "Separation Location (Upper): None\n"
            if lower_sep > 0:
                report += f"Separation Location (Lower): {lower_sep:.4f} x/c\n"
            else:
                report += "Separation Location (Lower): None\n"
        else:
            # For cascade with single surface
            if sep_loc > 0:
                report += f"Separation Location: {sep_loc:.4f} x/c\n"
            else:
                report += "Separation Location: None\n"
    
    # Additional cascade performance metrics
    if 'loss_coefficient' in solution:
        report += "\nCascade Performance\n"
        report += "-------------------\n"
        report += f"Total Pressure Loss Coefficient: {solution['loss_coefficient']:.6f}\n"
    if 'diffusion_factor' in solution:
        report += f"Diffusion Factor: {solution['diffusion_factor']:.4f}\n"
    if 'outlet_angle' in solution:
        report += f"Outlet Flow Angle: {np.degrees(solution['outlet_angle']):.2f}°\n"
    
    # Save to file if requested
    if output_file is not None:
        with open(output_file, 'w') as f:
            f.write(report)
    
    return report