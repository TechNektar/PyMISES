"""
Physics models for PyMISES.

This module contains the physical models used in PyMISES:
- Thermodynamic relations
- Laminar closure models
- Turbulent closure models
- Transition prediction
- Artificial dissipation models
"""

# Import key components
from pymises.physics.thermo import Thermodynamics
from pymises.physics.laminar import LaminarClosure
from pymises.physics.turbulent import TurbulentClosure
from pymises.physics.transition import TransitionModel, ModifiedAGSTransition, EnvelopeEnMethod, TransitionPredictor, create_transition_model
from pymises.physics.dissipation import ArtificialDissipation, create_dissipation_model

# Define exports
__all__ = [
    # Thermodynamic class
    'Thermodynamics',
    
    # Closure models
    'LaminarClosure',
    'TurbulentClosure',
    
    # Transition models
    'TransitionModel',
    'ModifiedAGSTransition',
    'EnvelopeEnMethod',
    'TransitionPredictor',
    'create_transition_model',
    
    # Artificial dissipation
    'ArtificialDissipation',
    'create_dissipation_model'
]