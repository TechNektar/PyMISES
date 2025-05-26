"""
PyMISES - Simple Panel Method Example

This example demonstrates how to use a simple panel method to analyze an airfoil.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import logging
from scipy.linalg import solve

# Add the PyMISES directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from pymises.core.geometry import AirfoilGeometry

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def panel_method(airfoil, alpha_deg):
    """
    Simple panel method for airfoil analysis.
    
    Parameters
    ----------
    airfoil : AirfoilGeometry
        Airfoil geometry
    alpha_deg : float
        Angle of attack in degrees
        
    Returns
    -------
    dict
        Dictionary containing results
    """
    # Convert angle of attack to radians
    alpha = np.radians(alpha_deg)
    
    # Get airfoil coordinates
    x = airfoil.x
    y = airfoil.y
    
    # Number of panels
    n = len(x) - 1
    
    # Panel endpoints
    x_p = x
    y_p = y
    
    # Panel midpoints
    x_m = 0.5 * (x_p[:-1] + x_p[1:])
    y_m = 0.5 * (y_p[:-1] + y_p[1:])
    
    # Panel lengths and angles
    dx = x_p[1:] - x_p[:-1]
    dy = y_p[1:] - y_p[:-1]
    s = np.sqrt(dx**2 + dy**2)
    theta = np.arctan2(dy, dx)
    
    # Normal vectors
    nx = np.sin(theta)
    ny = -np.cos(theta)
    
    # Influence coefficient matrices
    A = np.zeros((n, n))
    B = np.zeros((n, n))
    
    # Compute influence coefficients
    for i in range(n):
        for j in range(n):
            if i == j:
                # Self-influence
                A[i, j] = 0.0
                B[i, j] = 0.5
            else:
                # Compute influence of panel j on panel i
                xr = x_m[i] - x_p[j]
                yr = y_m[i] - y_p[j]
                xs = x_m[i] - x_p[j+1]
                ys = y_m[i] - y_p[j+1]
                
                r1 = np.sqrt(xr**2 + yr**2)
                r2 = np.sqrt(xs**2 + ys**2)
                
                beta = np.arctan2(yr*xs - ys*xr, xr*xs + yr*ys)
                
                # Source influence
                A[i, j] = (beta / (2 * np.pi))
                
                # Vortex influence
                B[i, j] = (np.log(r1/r2) / (2 * np.pi))
    
    # Right-hand side
    RHS = np.zeros(n)
    for i in range(n):
        RHS[i] = -(nx[i] * np.cos(alpha) + ny[i] * np.sin(alpha))
    
    # Add Kutta condition
    A_kutta = np.zeros(n+1)
    A_kutta[0] = 1.0
    A_kutta[n-1] = 1.0
    
    # Augment matrices
    A_aug = np.zeros((n+1, n+1))
    A_aug[:n, :n] = A
    A_aug[:n, n] = B.sum(axis=1)
    A_aug[n, :] = A_kutta
    
    RHS_aug = np.zeros(n+1)
    RHS_aug[:n] = RHS
    
    # Solve system
    solution = solve(A_aug, RHS_aug)
    
    # Extract source and vortex strengths
    sigma = solution[:n]
    gamma = solution[n]
    
    # Compute tangential velocity
    V_t = np.zeros(n)
    for i in range(n):
        V_t[i] = np.cos(alpha - theta[i]) + np.sum(sigma * A[i, :]) + gamma * np.sum(B[i, :])
    
    # Compute pressure coefficient
    Cp = 1.0 - V_t**2
    
    # Compute forces
    cl = 2 * gamma
    cm = 0.0
    for i in range(n):
        cm -= Cp[i] * (x_m[i] - 0.25) * s[i]
    
    # Return results
    return {
        'x': x_m,
        'y': y_m,
        'Cp': Cp,
        'cl': cl,
        'cm': cm,
        'alpha': alpha_deg
    }

def run_panel_method_validation(airfoil_name='naca0012', alpha_range=None):
    """
    Run panel method validation for an airfoil across a range of angles of attack.
    
    Parameters
    ----------
    airfoil_name : str
        Name of the airfoil ('naca0012', 'naca4412', etc.)
    alpha_range : list or None
        List of angles of attack to analyze
        If None, defaults to [-4, -2, 0, 2, 4, 6, 8, 10]
        
    Returns
    -------
    dict
        Dictionary containing all results
    """
    print(f"Running panel method validation for {airfoil_name}")
    
    # Default alpha range if not specified
    if alpha_range is None:
        alpha_range = [-4, -2, 0, 2, 4, 6, 8, 10]
    
    # 1. Create airfoil geometry
    if airfoil_name.lower().startswith('naca'):
        # Use NACA airfoil generator
        airfoil = AirfoilGeometry.create_naca(airfoil_name, n_points=101)
    else:
        # Load coordinates from file
        try:
            airfoil = AirfoilGeometry.import_from_file(f'{airfoil_name}.dat')
        except:
            print(f"Could not load {airfoil_name}.dat, using NACA 0012 instead")
            airfoil = AirfoilGeometry.create_naca('0012', n_points=101)
    
    # Initialize results dictionary
    results = {
        'airfoil': airfoil_name,
        'alpha': [],
        'cl': [],
        'cm': [],
        'solutions': []
    }
    
    # Run analysis for each angle of attack
    for alpha in alpha_range:
        print(f"Analyzing α = {alpha}°")
        
        # Run panel method
        solution = panel_method(airfoil, alpha)
        
        # Store results
        results['alpha'].append(alpha)
        results['cl'].append(solution['cl'])
        results['cm'].append(solution['cm'])
        results['solutions'].append(solution)
    
    # Create lift curve plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot lift curve
    ax1.plot(results['alpha'], results['cl'], 'bo-')
    ax1.set_xlabel('Angle of Attack (deg)')
    ax1.set_ylabel('Lift Coefficient (CL)')
    ax1.set_title(f'{airfoil_name.upper()} - Lift Curve')
    ax1.grid(True)
    
    # Plot moment curve
    ax2.plot(results['alpha'], results['cm'], 'ro-')
    ax2.set_xlabel('Angle of Attack (deg)')
    ax2.set_ylabel('Moment Coefficient (CM)')
    ax2.set_title(f'{airfoil_name.upper()} - Moment Curve')
    ax2.grid(True)
    
    fig.tight_layout()
    fig.savefig(f'{airfoil_name}_panel_method_forces.png', dpi=300)
    
    # Create pressure distribution plots for selected angles
    key_angles = [0, 4, 8]
    fig2, axes = plt.subplots(len(key_angles), 1, figsize=(8, 10), sharex=True)
    
    for i, alpha in enumerate(key_angles):
        if alpha in results['alpha']:
            idx = results['alpha'].index(alpha)
            solution = results['solutions'][idx]
            
            # Plot pressure distribution
            axes[i].plot(solution['x'], -solution['Cp'], 'b-')
            axes[i].set_ylabel('-Cp')
            axes[i].set_title(f'α = {alpha}°')
            axes[i].grid(True)
            axes[i].set_ylim(-1.5, 1.5)
    
    axes[-1].set_xlabel('x/c')
    fig2.suptitle(f'{airfoil_name.upper()} - Pressure Distributions')
    fig2.tight_layout()
    fig2.savefig(f'{airfoil_name}_panel_method_pressure.png', dpi=300)
    
    # Print summary of results
    print("\nPanel Method Validation Summary:")
    print(f"{'Alpha':^10} | {'CL':^10} | {'CM':^10}")
    print(f"{'-'*10} | {'-'*10} | {'-'*10}")
    
    for i, alpha in enumerate(results['alpha']):
        print(f"{alpha:^10.1f} | {results['cl'][i]:^10.4f} | {results['cm'][i]:^10.4f}")
    
    # Show plots
    plt.show()
    
    return results

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="PyMISES Simple Panel Method Example")
    parser.add_argument('--airfoil', type=str, default='naca0012', help="Airfoil name")
    parser.add_argument('--alpha-min', type=float, default=-4.0, help="Minimum angle of attack")
    parser.add_argument('--alpha-max', type=float, default=10.0, help="Maximum angle of attack")
    parser.add_argument('--alpha-step', type=float, default=2.0, help="Angle of attack step size")
    
    args = parser.parse_args()
    
    # Create alpha range
    alpha_range = np.arange(args.alpha_min, args.alpha_max + 0.1, args.alpha_step).tolist()
    
    # Run panel method validation
    results = run_panel_method_validation(
        airfoil_name=args.airfoil,
        alpha_range=alpha_range
    )
