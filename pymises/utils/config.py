"""
Configuration management for PyMISES.

This module provides utilities for managing configuration settings
throughout the PyMISES software. It allows for loading configuration
from files, setting default values, and accessing configuration in
a consistent way.
"""

import os
import json
import logging
from typing import Dict, Any, Optional, Union

logger = logging.getLogger(__name__)

class Config:
    """
    Configuration manager for PyMISES.
    
    This class handles loading, saving, and accessing configuration settings
    used throughout PyMISES. It provides sensible defaults and validation for
    required parameters.
    
    Attributes:
        settings (Dict[str, Any]): Configuration settings dictionary.
        config_path (str): Path to the configuration file.
    """
    
    # Default configuration values
    _defaults = {
        # Euler solver settings
        "euler": {
            "artificial_dissipation": {
                "threshold_mach": 0.95,  # Threshold Mach number for artificial dissipation
                "dissipation_coeff": 0.5,  # Dissipation coefficient
            },
            "relaxation_factor": 0.8,  # Relaxation factor for Newton iterations
            "max_iterations": 100,  # Maximum number of iterations
            "convergence_tolerance": 1e-6,  # Convergence tolerance
        },
        
        # Boundary layer settings
        "boundary_layer": {
            "laminar_params": {
                "falkner_skan_coeffs": [1.515, 0.076, 0.040],  # Falkner-Skan profile coefficients
            },
            "turbulent_params": {
                "lag_constant": 5.6,  # Lag constant for non-equilibrium effects
            },
            "transition": {
                "amplification_constant": 0.10,  # Amplification rate constant (A)
                "ramp_width": 0.30,  # Ramp width parameter (B)
                "turbulence_level": 0.01,  # Default turbulence level (1%)
            },
        },
        
        # Grid settings
        "grid": {
            "nodes_i": 101,  # Number of nodes in i-direction
            "nodes_j": 33,  # Number of nodes in j-direction
            "far_field_distance": 20.0,  # Far-field distance in chord lengths
            "leading_edge_clustering": 0.2,  # Clustering factor at leading edge
            "trailing_edge_clustering": 0.3,  # Clustering factor at trailing edge
        },
        
        # Newton solver settings
        "newton": {
            "max_iterations": 20,  # Maximum number of Newton iterations
            "convergence_tolerance": 1e-6,  # Convergence tolerance
            "min_relaxation": 0.1,  # Minimum relaxation factor
            "max_relaxation": 1.0,  # Maximum relaxation factor
            "adaptive_relaxation": True,  # Use adaptive relaxation
        },
        
        # Output settings
        "output": {
            "verbose": True,  # Verbose output
            "save_interim_results": False,  # Save intermediate results
            "visualization": True,  # Generate visualization
        },
    }
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize the configuration manager.
        
        Args:
            config_path: Path to the configuration file (optional).
                         If not provided, will use default values.
        """
        self.settings = self._defaults.copy()
        self.config_path = config_path
        
        if config_path and os.path.exists(config_path):
            self.load(config_path)
    
    def load(self, config_path: str) -> None:
        """
        Load configuration from a JSON file.
        
        Args:
            config_path: Path to the configuration file.
        """
        try:
            with open(config_path, 'r') as f:
                user_config = json.load(f)
            
            # Update the default settings with user-provided settings
            self._update_nested_dict(self.settings, user_config)
            logger.info(f"Loaded configuration from {config_path}")
        except Exception as e:
            logger.error(f"Error loading configuration: {str(e)}")
            logger.info("Using default configuration")
    
    def save(self, config_path: Optional[str] = None) -> None:
        """
        Save the current configuration to a JSON file.
        
        Args:
            config_path: Path where to save the configuration file.
                         If not provided, will use the path from initialization.
        """
        save_path = config_path or self.config_path
        if not save_path:
            logger.warning("No configuration path specified, cannot save")
            return
        
        try:
            with open(save_path, 'w') as f:
                json.dump(self.settings, f, indent=4)
            logger.info(f"Saved configuration to {save_path}")
        except Exception as e:
            logger.error(f"Error saving configuration: {str(e)}")
    
    def get(self, section: str, key: Optional[str] = None, default: Any = None) -> Any:
        """
        Get a configuration value.
        
        Args:
            section: Section name in the configuration.
            key: Key within the section. If None, returns the entire section.
            default: Default value to return if the key is not found.
        
        Returns:
            The configuration value, or default if not found.
        """
        if section not in self.settings:
            return default
        
        if key is None:
            return self.settings[section]
        
        return self.settings[section].get(key, default)
    
    def set(self, section: str, key: str, value: Any) -> None:
        """
        Set a configuration value.
        
        Args:
            section: Section name in the configuration.
            key: Key within the section.
            value: Value to set.
        """
        if section not in self.settings:
            self.settings[section] = {}
        
        self.settings[section][key] = value
    
    def update(self, updates: Dict[str, Any]) -> None:
        """
        Update multiple configuration values at once.
        
        Args:
            updates: Dictionary containing configuration updates.
        """
        self._update_nested_dict(self.settings, updates)
    
    def _update_nested_dict(self, d: Dict[str, Any], u: Dict[str, Any]) -> None:
        """
        Update a nested dictionary with values from another nested dictionary.
        
        Args:
            d: Dictionary to update (modified in-place).
            u: Dictionary with updates.
        """
        for k, v in u.items():
            if isinstance(v, dict) and k in d and isinstance(d[k], dict):
                self._update_nested_dict(d[k], v)
            else:
                d[k] = v
    
    def __getitem__(self, key: str) -> Any:
        """
        Allow dictionary-like access to configuration sections.
        
        Args:
            key: Section name.
            
        Returns:
            The configuration section.
        """
        return self.settings.get(key, {})
    
    def __setitem__(self, key: str, value: Any) -> None:
        """
        Allow dictionary-like setting of configuration sections.
        
        Args:
            key: Section name.
            value: Section value (typically a dict).
        """
        self.settings[key] = value

# Create a global configuration instance with default values
default_config = Config()

def get_config(config_path: Optional[str] = None) -> Config:
    """
    Get a configuration instance.
    
    Args:
        config_path: Path to the configuration file. 
                    If provided, will load configuration from this file.
                    
    Returns:
        A Configuration instance.
    """
    if config_path:
        return Config(config_path)
    return default_config
