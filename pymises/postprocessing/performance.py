"""
PyMISES - Performance Metrics Module

This module provides functions for calculating performance metrics from simulation results,
including aerodynamic coefficients, boundary layer parameters, and cascade performance.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Union

def calculate_airfoil_forces(solution: Dict, geometry) -> Dict:
    """
    Calculate aerodynamic forces and coefficients for an airfoil.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing pressure and grid data
    geometry : AirfoilGeometry
        Airfoil geometry object
        
    Returns
    -------
    Dict
        Dictionary containing lift, drag, and moment coefficients
    """
    # Extract solution variables
    pressure = solution['pressure']
    grid_x = solution['grid_x']
    grid_y = solution['grid_y']
    
    # Get airfoil surface indices
    surface_indices = geometry.get_surface_indices()
    
    # Extract flow conditions
    if 'mach' in solution and 'p_inf' in solution and 'alpha' in solution:
        mach = solution['mach']
        p_inf = solution['p_inf']
        alpha = solution['alpha']
    else:
        # Default values if not specified
        mach = 0.5
        p_inf = 101325.0  # Pa
        alpha = 0.0  # radians
    
    # Dynamic pressure
    q_inf = 0.5 * p_inf * 1.4 * mach**2
    
    # Allocate arrays for surface data
    x = grid_x[surface_indices]
    y = grid_y[surface_indices]
    p = pressure[surface_indices]
    
    # Calculate surface normals
    nx = np.zeros_like(x)
    ny = np.zeros_like(y)
    ds = np.zeros_like(x)
    
    # Calculate surface lengths and normals
    for i in range(len(x)):
        if i == 0:
            # Forward difference at leading edge
            dx = x[i+1] - x[i]
            dy = y[i+1] - y[i]
        elif i == len(x) - 1:
            # Backward difference at trailing edge
            dx = x[i] - x[i-1]
            dy = y[i] - y[i-1]
        else:
            # Central difference everywhere else
            dx = 0.5 * (x[i+1] - x[i-1])
            dy = 0.5 * (y[i+1] - y[i-1])
        
        # Surface length element
        ds[i] = np.sqrt(dx**2 + dy**2)
        
        # Outward normal (rotate tangent vector 90 degrees)
        nx[i] = dy / ds[i]
        ny[i] = -dx / ds[i]
    
    # Calculate pressure forces in x and y directions
    fx = np.sum((p - p_inf) * nx * ds)
    fy = np.sum((p - p_inf) * ny * ds)
    
    # Rotate forces to get lift and drag in wind axes
    cl = (fy * np.cos(alpha) - fx * np.sin(alpha)) / q_inf
    cd = (fx * np.cos(alpha) + fy * np.sin(alpha)) / q_inf
    
    # Calculate moment around quarter-chord
    xref = 0.25  # Quarter-chord reference point
    yref = 0.0
    
    moment = 0.0
    for i in range(len(x)):
        # Force components at this point
        fx_i = (p[i] - p_inf) * nx[i] * ds[i]
        fy_i = (p[i] - p_inf) * ny[i] * ds[i]
        
        # Moment arm from reference point
        dx = x[i] - xref
        dy = y[i] - yref
        
        # Contribution to moment (positive nose-up)
        moment += dx * fy_i - dy * fx_i
    
    # Normalize by dynamic pressure and chord
    cm = -moment / q_inf  # Negative sign for conventional aerodynamic sign convention
    
    return {
        'cl': cl,
        'cd': cd,
        'cm': cm,
        'fx': fx,
        'fy': fy
    }

def calculate_cascade_performance(solution: Dict, geometry, pitch: float) -> Dict:
    """
    Calculate performance metrics for a cascade.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing flow data
    geometry : BladeGeometry
        Blade geometry object
    pitch : float
        Blade-to-blade pitch
        
    Returns
    -------
    Dict
        Dictionary containing cascade performance metrics
    """
    # Extract inlet and outlet flow properties
    if 'inlet_data' in solution and 'outlet_data' in solution:
        inlet_data = solution['inlet_data']
        outlet_data = solution['outlet_data']
    else:
        # Extract from raw solution (simplified)
        inlet_data = extract_inlet_data(solution)
        outlet_data = extract_outlet_data(solution)
    
    # Calculate key performance parameters
    
    # Total pressure loss coefficient
    # ω = (p_t1 - p_t2) / (p_t1 - p_1)
    p_t1 = inlet_data['total_pressure']
    p_t2 = outlet_data['total_pressure']
    p_1 = inlet_data['static_pressure']
    loss_coefficient = (p_t1 - p_t2) / (p_t1 - p_1)
    
    # Flow angles
    alpha_1 = inlet_data['flow_angle']  # Inlet flow angle (radians)
    alpha_2 = outlet_data['flow_angle']  # Outlet flow angle (radians)
    
    # Velocities
    v_1 = inlet_data['velocity_magnitude']
    v_2 = outlet_data['velocity_magnitude']
    
    # Calculate tangential velocity components
    v_theta1 = v_1 * np.sin(alpha_1)
    v_theta2 = v_2 * np.sin(alpha_2)
    
    # Calculate solidity (chord/pitch)
    chord = geometry.get_chord()
    solidity = chord / pitch
    
    # Diffusion factor (Lieblein's definition)
    # DF = 1 - (V_2/V_1) + (|V_θ1 - V_θ2|)/(2 * solidity * V_1)
    diffusion_factor = 1.0 - (v_2/v_1) + abs(v_theta1 - v_theta2) / (2.0 * solidity * v_1)
    
    # Deflection (change in flow angle)
    deflection = abs(alpha_1 - alpha_2)
    
    # Zweifel loading coefficient
    # Z = (2 * solidity * cos(alpha_2) / cos(alpha_1)) * (tan(alpha_1) - tan(alpha_2))
    zweifel = (2 * solidity * np.cos(alpha_2) / np.cos(alpha_1)) * (np.tan(alpha_1) - np.tan(alpha_2))
    
    # Calculate static pressure rise coefficient
    # Cp = (p_2 - p_1) / (0.5 * rho_1 * V_1^2)
    p_2 = outlet_data['static_pressure']
    rho_1 = inlet_data['density']
    pressure_rise_coefficient = (p_2 - p_1) / (0.5 * rho_1 * v_1**2)
    
    # Return all performance metrics
    return {
        'loss_coefficient': loss_coefficient,
        'diffusion_factor': diffusion_factor,
        'deflection': deflection,
        'zweifel_coefficient': zweifel,
        'pressure_rise_coefficient': pressure_rise_coefficient,
        'inlet_flow_angle': alpha_1,
        'outlet_flow_angle': alpha_2,
        'inlet_mach': inlet_data['mach'],
        'outlet_mach': outlet_data['mach'],
        'pressure_ratio': p_2 / p_1,
        'total_pressure_ratio': p_t2 / p_t1
    }

def extract_inlet_data(solution: Dict) -> Dict:
    """
    Extract inlet flow data from solution.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing flow data
        
    Returns
    -------
    Dict
        Dictionary containing inlet flow properties
    """
    # This is a simplified implementation
    # In a real implementation, we would extract data at the inlet boundary
    
    # Get inlet boundary indices
    if 'inlet_indices' in solution:
        inlet_indices = solution['inlet_indices']
    else:
        # Simplified approach: assume inlet is at i=0
        ni = solution['grid_ni']
        nj = solution['grid_nj']
        inlet_indices = [i * nj for i in range(ni)]
    
    # Extract flow variables at inlet
    rho = np.mean(solution['density'][inlet_indices])
    vx = np.mean(solution['velocity_x'][inlet_indices])
    vy = np.mean(solution['velocity_y'][inlet_indices])
    p = np.mean(solution['pressure'][inlet_indices])
    
    # Calculate derived quantities
    v_magnitude = np.sqrt(vx**2 + vy**2)
    flow_angle = np.arctan2(vy, vx)
    
    # Calculate Mach number
    if 'speed_of_sound' in solution:
        a = np.mean(solution['speed_of_sound'][inlet_indices])
    else:
        # Estimate speed of sound
        gamma = 1.4  # Ratio of specific heats for air
        a = np.sqrt(gamma * p / rho)
    
    mach = v_magnitude / a
    
    # Calculate total pressure
    gamma = 1.4
    total_pressure = p * (1 + (gamma-1)/2 * mach**2)**(gamma/(gamma-1))
    
    return {
        'density': rho,
        'velocity_x': vx,
        'velocity_y': vy,
        'velocity_magnitude': v_magnitude,
        'flow_angle': flow_angle,
        'static_pressure': p,
        'total_pressure': total_pressure,
        'mach': mach,
        'speed_of_sound': a
    }

def extract_outlet_data(solution: Dict) -> Dict:
    """
    Extract outlet flow data from solution.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing flow data
        
    Returns
    -------
    Dict
        Dictionary containing outlet flow properties
    """
    # This is a simplified implementation
    # In a real implementation, we would extract data at the outlet boundary
    
    # Get outlet boundary indices
    if 'outlet_indices' in solution:
        outlet_indices = solution['outlet_indices']
    else:
        # Simplified approach: assume outlet is at i=ni-1
        ni = solution['grid_ni']
        nj = solution['grid_nj']
        outlet_indices = [(ni-1) * nj + j for j in range(nj)]
    
    # Extract flow variables at outlet
    rho = np.mean(solution['density'][outlet_indices])
    vx = np.mean(solution['velocity_x'][outlet_indices])
    vy = np.mean(solution['velocity_y'][outlet_indices])
    p = np.mean(solution['pressure'][outlet_indices])
    
    # Calculate derived quantities
    v_magnitude = np.sqrt(vx**2 + vy**2)
    flow_angle = np.arctan2(vy, vx)
    
    # Calculate Mach number
    if 'speed_of_sound' in solution:
        a = np.mean(solution['speed_of_sound'][outlet_indices])
    else:
        # Estimate speed of sound
        gamma = 1.4  # Ratio of specific heats for air
        a = np.sqrt(gamma * p / rho)
    
    mach = v_magnitude / a
    
    # Calculate total pressure
    gamma = 1.4
    total_pressure = p * (1 + (gamma-1)/2 * mach**2)**(gamma/(gamma-1))
    
    return {
        'density': rho,
        'velocity_x': vx,
        'velocity_y': vy,
        'velocity_magnitude': v_magnitude,
        'flow_angle': flow_angle,
        'static_pressure': p,
        'total_pressure': total_pressure,
        'mach': mach,
        'speed_of_sound': a
    }

def calculate_boundary_layer_parameters(solution: Dict) -> Dict:
    """
    Calculate boundary layer parameters from solution.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing boundary layer data
        
    Returns
    -------
    Dict
        Dictionary containing boundary layer parameters
    """
    # Extract boundary layer data
    if 'displacement_thickness' not in solution or 'momentum_thickness' not in solution:
        raise ValueError("Solution missing boundary layer data")
    
    delta_star = solution['displacement_thickness']
    theta = solution['momentum_thickness']
    
    # Calculate shape factor
    H = delta_star / theta
    
    # Identify transition and separation points
    # In the boundary layer, separation is typically identified where H > 3.0
    # and transition where there's a sudden drop in H
    
    # For airfoil, we separate upper and lower surfaces
    if 'surface_type' in solution:
        surface_type = solution['surface_type']
        upper_indices = np.where(surface_type == 'upper')[0]
        lower_indices = np.where(surface_type == 'lower')[0]
        
        # Calculate maximum H along each surface
        max_H_upper = np.max(H[upper_indices]) if len(upper_indices) > 0 else 0
        max_H_lower = np.max(H[lower_indices]) if len(lower_indices) > 0 else 0
        
        # Find potential separation points (H > 3.0)
        sep_upper_idx = np.where((H[upper_indices] > 3.0) & (H[upper_indices] >= 0.95 * max_H_upper))[0]
        sep_lower_idx = np.where((H[lower_indices] > 3.0) & (H[lower_indices] >= 0.95 * max_H_lower))[0]
        
        # Get x locations
        x = solution['x'] if 'x' in solution else np.arange(len(delta_star)) / (len(delta_star) - 1)
        
        # Determine separation points
        if len(sep_upper_idx) > 0:
            x_sep_upper = x[upper_indices[sep_upper_idx[0]]]
        else:
            x_sep_upper = -1  # No separation
        
        if len(sep_lower_idx) > 0:
            x_sep_lower = x[lower_indices[sep_lower_idx[0]]]
        else:
            x_sep_lower = -1  # No separation
        
        # For transition, we look for rapid drops in H
        # This is a simplified approach; a real implementation would use the transition flag
        # from the boundary layer solver
        
        # Return parameters
        return {
            'shape_factor': H,
            'max_shape_factor_upper': max_H_upper,
            'max_shape_factor_lower': max_H_lower,
            'separation_location': {
                'upper': x_sep_upper,
                'lower': x_sep_lower
            }
        }
    else:
        # For a single surface (e.g., cascade blade)
        max_H = np.max(H)
        
        # Find potential separation point
        sep_idx = np.where((H > 3.0) & (H >= 0.95 * max_H))[0]
        
        # Get x location
        x = solution['x'] if 'x' in solution else np.arange(len(delta_star)) / (len(delta_star) - 1)
        
        # Determine separation point
        if len(sep_idx) > 0:
            x_sep = x[sep_idx[0]]
        else:
            x_sep = -1  # No separation
        
        # Return parameters
        return {
            'shape_factor': H,
            'max_shape_factor': max_H,
            'separation_location': x_sep
        }

def calculate_energy_dissipation(solution: Dict) -> Dict:
    """
    Calculate energy dissipation in the flow field.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing flow data
        
    Returns
    -------
    Dict
        Dictionary containing energy dissipation metrics
    """
    # Extract flow variables
    rho = solution['density']
    vx = solution['velocity_x']
    vy = solution['velocity_y']
    p = solution['pressure']
    
    # Calculate velocity magnitude
    v_magnitude = np.sqrt(vx**2 + vy**2)
    
    # Calculate total pressure
    if 'speed_of_sound' in solution:
        a = solution['speed_of_sound']
        mach = v_magnitude / a
    else:
        # Estimate Mach number
        gamma = 1.4
        a = np.sqrt(gamma * p / rho)
        mach = v_magnitude / a
    
    # Calculate total pressure
    gamma = 1.4
    p_total = p * (1 + (gamma-1)/2 * mach**2)**(gamma/(gamma-1))
    
    # Calculate entropy generation (approximation)
    # ds = (gamma-1)/(gamma) * R * ln(p2/p1 * (rho1/rho2)^gamma)
    # For small changes, we can use gradients
    # Here we just compute total pressure loss, which is related to entropy generation
    
    # Get reference values (usually from inlet)
    p_total_ref = np.max(p_total)
    
    # Calculate total pressure loss relative to reference
    total_pressure_loss = (p_total_ref - p_total) / p_total_ref
    
    # Calculate volume-averaged total pressure loss
    avg_total_pressure_loss = np.mean(total_pressure_loss)
    max_total_pressure_loss = np.max(total_pressure_loss)
    
    return {
        'total_pressure': p_total,
        'total_pressure_loss': total_pressure_loss,
        'avg_total_pressure_loss': avg_total_pressure_loss,
        'max_total_pressure_loss': max_total_pressure_loss
    }

def analyze_wake(solution: Dict, wake_location: float = 1.0) -> Dict:
    """
    Analyze wake properties at a specific downstream location.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing flow data
    wake_location : float
        Downstream location (x/c) for wake analysis
        
    Returns
    -------
    Dict
        Dictionary containing wake properties
    """
    # Extract grid information
    if 'grid_x' not in solution or 'grid_y' not in solution:
        raise ValueError("Solution missing grid data")
    
    x = solution['grid_x']
    y = solution['grid_y']
    
    # Find points at the specified downstream location
    tolerance = 0.01  # Tolerance for finding points at wake_location
    wake_indices = np.where(np.abs(x - wake_location) < tolerance)[0]
    
    if len(wake_indices) == 0:
        raise ValueError(f"No points found at downstream location x/c = {wake_location}")
    
    # Extract flow variables at wake
    rho_wake = solution['density'][wake_indices]
    vx_wake = solution['velocity_x'][wake_indices]
    vy_wake = solution['velocity_y'][wake_indices]
    p_wake = solution['pressure'][wake_indices]
    
    # Calculate velocity magnitude at wake
    v_wake = np.sqrt(vx_wake**2 + vy_wake**2)
    
    # Get y-coordinates at wake
    y_wake = y[wake_indices]
    
    # Sort points by y-coordinate
    sort_idx = np.argsort(y_wake)
    y_wake = y_wake[sort_idx]
    v_wake = v_wake[sort_idx]
    vx_wake = vx_wake[sort_idx]
    
    # Find wake center (minimum velocity)
    wake_center_idx = np.argmin(v_wake)
    wake_center_y = y_wake[wake_center_idx]
    wake_deficit = 1.0 - v_wake[wake_center_idx] / np.max(v_wake)
    
    # Calculate wake width (distance between half-deficit points)
    half_deficit = wake_deficit / 2.0
    half_deficit_indices = np.where(1.0 - v_wake / np.max(v_wake) >= half_deficit)[0]
    
    if len(half_deficit_indices) >= 2:
        wake_width = y_wake[half_deficit_indices[-1]] - y_wake[half_deficit_indices[0]]
    else:
        wake_width = 0.0
    
    # Calculate wake momentum thickness
    # θ = ∫(1 - u/U)(u/U) dy
    y_spacing = np.diff(y_wake)
    y_spacing = np.append(y_spacing, y_spacing[-1])  # Extend for last point
    
    # Normalize velocity by maximum
    u_ratio = v_wake / np.max(v_wake)
    integrand = (1.0 - u_ratio) * u_ratio
    
    # Integrate using trapezoidal rule
    wake_theta = np.sum(integrand * y_spacing)
    
    return {
        'wake_center_y': wake_center_y,
        'wake_deficit': wake_deficit,
        'wake_width': wake_width,
        'wake_momentum_thickness': wake_theta,
        'wake_y': y_wake,
        'wake_velocity': v_wake,
        'wake_velocity_x': vx_wake
    }