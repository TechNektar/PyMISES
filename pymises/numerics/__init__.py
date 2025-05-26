"""
PyMISES Numerical Methods Package

This package provides numerical methods for solving equations, calculating
derivatives, and handling matrix operations efficiently.
"""

from pymises.core.newton import NewtonSolver
from pymises.numerics.jacobian import JacobianCalculator, JacobianAssembler
from pymises.numerics.linear_solver import LinearSolver
from pymises.utils.validation import validate_jacobian

__all__ = [
    'NewtonSolver',
    'JacobianCalculator',
    'JacobianAssembler',
    'LinearSolver',
    'validate_jacobian'
]