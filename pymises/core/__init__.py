"""
Core solvers and components for PyMISES.

This module contains the primary components of the PyMISES solver:
- Euler equation solver
- Boundary layer solver
- Viscous-inviscid coupling
- Newton iteration framework
- Grid generation and management
- Airfoil/cascade geometry handling
"""

# Import key components
from pymises.core.euler import EulerSolver
from pymises.core.boundary_layer import BoundaryLayerSolver, BoundaryLayerFactory
from pymises.core.coupling import CoupledSolver
from pymises.core.newton import NewtonSolver
from pymises.core.grid import GridGenerator, StreamlineGrid
from pymises.core.geometry import BladeGeometry, AirfoilGeometry

# Define exports
__all__ = [
    'EulerSolver',
    'BoundaryLayerSolver',
    'BoundaryLayerFactory',
    'CoupledSolver',
    'NewtonSolver',
    'GridGenerator',
    'StreamlineGrid',
    'BladeGeometry',
    'AirfoilGeometry'
]