"""
PyMISES - Boundary Conditions Package

This package provides implementations of various boundary conditions for the PyMISES CFD solver.
"""

# Wall boundary conditions
from pymises.boundary_conditions.wall import (
    WallBoundaryCondition,
    InviscidWallBC,
    ViscousWallBC
)

# Farfield boundary conditions
from pymises.boundary_conditions.farfield import (
    FarfieldBC,
    SubsonicInflow,
    SubsonicOutflow,
    VortexFarfieldBC,
    CharacteristicBC
)

# Periodicity boundary conditions
from pymises.boundary_conditions.periodicity import (
    PeriodicityBC,
    PhaseLagPeriodicityBC
)

# Inverse design boundary conditions
from pymises.boundary_conditions.inverse import (
    InverseDesignBC,
    PressureSpecificationBC,
    GeometricConstraintBC,
    MixedInverseBC,
    ModalInverseDesignBC
)

# Define package exports
__all__ = [
    # Wall boundary conditions
    'WallBoundaryCondition',
    'InviscidWallBC',
    'ViscousWallBC',
    
    # Farfield boundary conditions
    'FarfieldBC',
    'SubsonicInflow',
    'SubsonicOutflow',
    'VortexFarfieldBC',
    'CharacteristicBC',
    
    # Periodicity boundary conditions
    'PeriodicityBC',
    'PhaseLagPeriodicityBC',
    
    # Inverse design boundary conditions
    'InverseDesignBC',
    'PressureSpecificationBC',
    'GeometricConstraintBC',
    'MixedInverseBC',
    'ModalInverseDesignBC'
]