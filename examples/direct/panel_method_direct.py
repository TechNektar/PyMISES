"""
PyMISES - Panel Method Direct Example

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

class PanelMethod:
    """
    Simple panel method for airfoil analysis.
    """
    
    def __init__(self, airfoil):
        """
        Initialize the panel method.
        
        Parameters
        ----------
        airfoil : AirfoilGeometry
            Airfoil geometry
        """
        self.airfoil = airfoil
        
        # Get airfoil coordinates
        self.x = airfoil.x
        self.y = airfoil.y
        
        # Number of panels
        self.n = len(self.x) - 1
        
        # Panel endpoints
        self.x_p = self.x
        self.y_p = self.y
        
        # Panel midpoints
        self.x_m = 0.5 * (self.x_p[:-1] + self.x_p[1:])
        self.y_m = 0.5 * (self.y_p[:-1] + self.y_p[1:])
        
        # Panel lengths and angles
        self.dx = self.x_p[1:] - self.x_p[:-1]
        self.dy = self.y_p[1:] - self.y_p[:-1]
        self.s = np.sqrt(self.dx**2 + self.dy**2)
        self.theta = np.arctan2(self.dy, self.dx)
        
        # Normal vectors
        self.nx = np.sin(self.theta)
        self.ny = -np.cos(self.theta)
        
        # Influence coefficient matrices
        self.A = np.zeros((self.n, self.n))
        self.B = np.zeros((self.n, self.n))
        
        # Compute influence coefficients
        self._compute_influence_coefficients()
    
    def _compute_influence_coefficients(self):
        """
        Compute influence coefficients for the panel method.
        """
        for i in range(self.n):
            for j in range(self.n):
                if i == j:
                    # Self-influence
                    self.A[i, j] = 0.0
                    self.B[i, j] = 0.5
                else:
                    # Compute influence of panel j on panel i
                    xr = self.x_m[i] - self.x_p[j]
                    yr = self.y_m[i] - self.y_p[j]
                    xs = self.x_m[i] - self.x_p[j+1]
                    ys = self.y_m[i] - self.y_p[j+1]
                    
                    r1 = np.sqrt(xr**2 + yr**2)
                    r2 = np.sqrt(xs**2 + ys**2)
                    
                    beta = np.arctan2(yr*xs - ys*xr, xr*xs + yr*ys)
                    
                    # Source influence
                    self.A[i, j] = (beta / (2 * np.pi))
                    
                    # Vortex influence
                    self.B[i, j] = (np.log(r1/r2) / (2 * np.pi))
    
    def solve(self, alpha_deg):
        """
        Solve the panel method for a given angle of attack.
        
        Parameters
        ----------
        alpha_deg : float
            Angle of attack in degrees
            
        Returns
        -------
        dict
            Dictionary containing results
        """
        # Convert angle of attack to radians
        alpha = np.radians(alpha_deg)
        
        # Right-hand side
        RHS = np.zeros(self.n)
        for i in range(self.n):
            RHS[i] = -(self.nx[i] * np.cos(alpha) + self.ny[i] * np.sin(alpha))
        
        # Add Kutta condition
        A_kutta = np.zeros(self.n+1)
        A_kutta[0] = 1.0
        A_kutta[self.n-1] = 1.0
        
        # Augment matrices
        A_aug = np.zeros((self.n+1, self.n+1))
        A_aug[:self.n, :self.n] = self.A
        A_aug[:self.n, self.n] = self.B.sum(axis=1)
        A_aug[self.n, :] = A_kutta
        
        RHS_aug = np.zeros(self.n+1)
        RHS_aug[:self.n] = RHS
        
        # Solve system
        solution = solve(A_aug, RHS_aug)
        
        # Extract source and vortex strengths
        sigma = solution[:self.n]
        gamma = solution[self.n]
        
        # Compute tangential velocity
        V_t = np.zeros(self.n)
        for i in range(self.n):
            V_t[i] = np.cos(alpha - self.theta[i]) + np.sum(sigma * self.A[i, :]) + gamma * np.sum(self.B[i, :])
        
        # Compute pressure coefficient
        Cp = 1.0 - V_t**2
        
        # Compute forces
        cl = 2 * gamma
        cm = 0.0
        for i in range(self.n):
            cm -= Cp[i] * (self.x_m[i] - 0.25) * self.s[i]
        
        # Return results
        return {
            'x': self.x_m,
            'y': self.y_m,
            'Cp': Cp,
            'cl': cl,
            'cm': cm,
            'alpha': alpha_deg,
            'V_t': V_t,
            'sigma': sigma,
            'gamma': gamma
        }

def run_panel_method_direct(airfoil_name='naca0012', alpha=4.0):
    """
    Run panel method direct analysis for an airfoil.
    
    Parameters
    ----------
    airfoil_name : str
        Name of the airfoil ('naca0012', 'naca4412', etc.)
    alpha : float
        Angle of attack in degrees
        
    Returns
    -------
    dict
        Dictionary containing results
    """
    print(f"Running panel method direct analysis for {airfoil_name} at α={alpha}°")
    
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
    
    # 2. Create panel method
    panel_method = PanelMethod(airfoil)
    
    # 3. Solve for the specified angle of attack
    solution = panel_method.solve(alpha)
    
    # 4. Create visualization plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot airfoil geometry
    ax1.plot(airfoil.x, airfoil.y, 'k-')
    ax1.set_xlabel('x/c')
    ax1.set_ylabel('y/c')
    ax1.set_title(f'{airfoil_name.upper()} Geometry')
    ax1.grid(True)
    ax1.set_aspect('equal')
    
    # Plot pressure distribution
    ax2.plot(solution['x'], -solution['Cp'], 'b-')
    ax2.set_xlabel('x/c')
    ax2.set_ylabel('-Cp')
    ax2.set_title(f'Pressure Distribution, α={alpha}°')
    ax2.grid(True)
    
    fig.tight_layout()
    fig.savefig(f'{airfoil_name}_panel_method.png', dpi=300)
    
    # Print summary of results
    print("\nPanel Method Results:")
    print(f"Lift coefficient (CL): {solution['cl']:.4f}")
    print(f"Moment coefficient (CM): {solution['cm']:.4f}")
    
    # Show plots
    plt.show()
    
    return solution

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="PyMISES Panel Method Direct Example")
    parser.add_argument('--airfoil', type=str, default='naca0012', help="Airfoil name")
    parser.add_argument('--alpha', type=float, default=4.0, help="Angle of attack in degrees")
    
    args = parser.parse_args()
    
    # Run panel method direct analysis
    solution = run_panel_method_direct(
        airfoil_name=args.airfoil,
        alpha=args.alpha
    )
