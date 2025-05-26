"""
PyMISES - Python Multiple-Blade Interacting Streamtube Euler Solver

A Python implementation of the MISES solver originally developed by Mark Drela at MIT.
This package provides tools for aerodynamic analysis and design of airfoils and turbine cascades
using a coupled Euler/boundary-layer approach.

Key features:
- Streamline-based Euler solver with accurate shock capturing
- Integral boundary layer solver with transition prediction
- Viscous-inviscid coupling for separated flows
- Inverse design capabilities
- Comprehensive post-processing tools

For more information, see the documentation at [TBD].
"""

__version__ = '0.1.0'
__author__ = 'PyMISES Development Team'

# Import key components for easier access
from pymises.core.geometry import BladeGeometry, AirfoilGeometry
from pymises.core.grid import GridGenerator, StreamlineGrid
from pymises.core.euler import EulerSolver
from pymises.core.boundary_layer import BoundaryLayerSolver, BoundaryLayerFactory
from pymises.core.coupling import CoupledSolver
from pymises.core.newton import NewtonSolver

# Define package exports
__all__ = [
    # Version info
    '__version__',
    '__author__',
    
    # Core components
    'BladeGeometry',
    'AirfoilGeometry',
    'GridGenerator',
    'StreamlineGrid',
    'EulerSolver',
    'BoundaryLayerSolver',
    'BoundaryLayerFactory',
    'CoupledSolver',
    'NewtonSolver'
]

# Optional imports for nicer user experience
try:
    # Import basic boundary conditions
    from pymises.boundary_conditions.wall import InviscidWallBC, ViscousWallBC
    from pymises.boundary_conditions.farfield import VortexFarfieldBC
    
    __all__.extend([
        'InviscidWallBC',
        'ViscousWallBC',
        'VortexFarfieldBC'
    ])
except ImportError:
    pass

try:
    # Import common post-processing functions
    from pymises.postprocessing.visualize import plot_pressure, plot_mach_contours
    from pymises.postprocessing.performance import calculate_airfoil_forces
    
    __all__.extend([
        'plot_pressure',
        'plot_mach_contours',
        'calculate_airfoil_forces'
    ])
except ImportError:
    pass