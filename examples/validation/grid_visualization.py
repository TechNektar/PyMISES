"""
PyMISES - Grid Visualization Example

This example demonstrates how to generate and visualize a grid around an airfoil.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import logging

# Add the PyMISES directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from pymises.core.geometry import AirfoilGeometry
from pymises.core.grid import GridGenerator

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def visualize_grid(airfoil_name='naca0012', grid_type='o-grid', ni=51, nj=31):
    """
    Generate and visualize a grid around an airfoil.
    
    Parameters
    ----------
    airfoil_name : str
        Name of the airfoil ('naca0012', 'naca4412', etc.)
    grid_type : str
        Type of grid to generate ('o-grid', 'c-grid', etc.)
    ni : int
        Number of points in the streamwise direction
    nj : int
        Number of points in the normal direction
    """
    print(f"Generating {grid_type} for {airfoil_name} with {ni}x{nj} points")
    
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
    
    # 2. Generate computational grid
    grid_gen = GridGenerator(airfoil, {
        'ni': ni,  # Number of points in streamwise direction
        'nj': nj,  # Number of points in normal direction
        'far_field_distance': 10.0,  # Far field boundary distance in chord lengths
        'le_clustering': 0.2,  # Leading edge clustering factor
        'te_clustering': 0.3,  # Trailing edge clustering factor
        'wall_clustering': 0.3,  # Wall clustering factor
        'clustering_method': 'tanh'  # Clustering method ('tanh', 'exp', 'sine')
    })
    
    # Generate grid
    grid = grid_gen.generate_grid(grid_type=grid_type)
    
    # 3. Visualize grid
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot grid lines
    for j in range(nj):
        ax.plot(grid.x[:, j], grid.y[:, j], 'k-', linewidth=0.5)
    
    for i in range(ni):
        ax.plot(grid.x[i, :], grid.y[i, :], 'k-', linewidth=0.5)
    
    # Plot airfoil
    ax.plot(airfoil.x, airfoil.y, 'r-', linewidth=2)
    
    # Set axis limits and labels
    ax.set_xlim(-0.2, 1.2)
    ax.set_ylim(-0.7, 0.7)
    ax.set_xlabel('x/c')
    ax.set_ylabel('y/c')
    ax.set_title(f'{airfoil_name.upper()} - {grid_type.upper()} ({ni}x{nj})')
    ax.set_aspect('equal')
    ax.grid(True)
    
    # Save figure
    fig.savefig(f'{airfoil_name}_{grid_type}_{ni}x{nj}.png', dpi=300)
    
    # Show figure
    plt.show()
    
    return grid

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="PyMISES Grid Visualization Example")
    parser.add_argument('--airfoil', type=str, default='naca0012', help="Airfoil name")
    parser.add_argument('--grid-type', type=str, default='o-grid', choices=['o-grid', 'c-grid'], help="Grid type")
    parser.add_argument('--ni', type=int, default=51, help="Number of points in streamwise direction")
    parser.add_argument('--nj', type=int, default=31, help="Number of points in normal direction")
    
    args = parser.parse_args()
    
    # Visualize grid
    grid = visualize_grid(
        airfoil_name=args.airfoil,
        grid_type=args.grid_type,
        ni=args.ni,
        nj=args.nj
    )
