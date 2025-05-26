"""
Artificial dissipation models for PyMISES.

This module provides implementations of artificial dissipation methods
used in the Euler solver to stabilize solutions in regions with shocks.
"""

import numpy as np
from typing import Dict, Tuple, Union, Optional, Callable, Any

from pymises.utils.logger import get_logger
from pymises.physics.thermo import Thermodynamics

logger = get_logger(__name__)

class ArtificialDissipation:
    """
    Base class for artificial dissipation models in PyMISES.
    
    This class provides the interface for artificial dissipation
    used in the Euler solver to capture shocks and prevent numerical
    instabilities.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the artificial dissipation model.
        
        Args:
            config: Configuration dictionary with dissipation parameters.
        """
        self.config = config or {}
    
    def compute_coefficient(self, mach: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Compute the artificial dissipation coefficient.
        
        Args:
            mach: Local Mach number.
            
        Returns:
            Dissipation coefficient.
        """
        raise NotImplementedError("Subclasses must implement compute_coefficient")
    
    def compute_jacobian(self, mach: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Compute the Jacobian of the artificial dissipation coefficient wrt Mach number.
        
        Args:
            mach: Local Mach number.
            
        Returns:
            Jacobian of dissipation coefficient wrt Mach number.
        """
        raise NotImplementedError("Subclasses must implement compute_jacobian")
    
    def apply(self, residuals: np.ndarray, solution: Dict[str, np.ndarray], 
             grid: Dict[str, np.ndarray]) -> np.ndarray:
        """
        Apply artificial dissipation to the residuals.
        
        Args:
            residuals: Current residual values.
            solution: Current solution variables.
            grid: Computational grid.
            
        Returns:
            Modified residuals with artificial dissipation applied.
        """
        raise NotImplementedError("Subclasses must implement apply")


class BulkViscosityDissipation(ArtificialDissipation):
    """
    Bulk viscosity-like artificial dissipation model.
    
    This model applies dissipation proportional to a nonlinear function
    of the local Mach number, activating only in transonic and supersonic
    regions to capture shocks.
    
    Attributes:
        threshold_mach: Mach number threshold above which dissipation is applied.
        dissipation_coeff: Coefficient controlling dissipation strength.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the bulk viscosity dissipation model.
        
        Args:
            config: Configuration dictionary with the following keys:
                - threshold_mach: Mach number threshold (default: 0.95)
                - dissipation_coeff: Dissipation coefficient (default: 0.5)
        """
        super().__init__(config)
        
        self.threshold_mach = self.config.get('threshold_mach', 0.95)
        self.dissipation_coeff = self.config.get('dissipation_coeff', 0.5)
        
        logger.info(f"Initialized bulk viscosity dissipation with threshold_mach={self.threshold_mach}, "
                   f"dissipation_coeff={self.dissipation_coeff}")
    
    def compute_coefficient(self, mach: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Compute the artificial dissipation coefficient based on local Mach number.
        
        Args:
            mach: Local Mach number.
            
        Returns:
            Dissipation coefficient.
        """
        # Only apply dissipation in regions where Mach number exceeds the threshold
        if np.isscalar(mach):
            if mach < self.threshold_mach:
                return 0.0
            else:
                return self.dissipation_coeff * (1.0 - (self.threshold_mach/mach)**2) / (2.0 * mach**2)
        else:
            result = np.zeros_like(mach)
            mask = mach >= self.threshold_mach
            result[mask] = self.dissipation_coeff * (1.0 - (self.threshold_mach/mach[mask])**2) / (2.0 * mach[mask]**2)
            return result
    
    def compute_jacobian(self, mach: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Compute the Jacobian of the artificial dissipation coefficient wrt Mach number.
        
        Args:
            mach: Local Mach number.
            
        Returns:
            Jacobian of dissipation coefficient wrt Mach number.
        """
        # Derivative of the dissipation coefficient with respect to Mach number
        if np.isscalar(mach):
            if mach < self.threshold_mach:
                return 0.0
            else:
                term1 = -self.dissipation_coeff / (2.0 * mach**3)
                term2 = self.dissipation_coeff * self.threshold_mach**2 / mach**4
                return term1 + term2
        else:
            result = np.zeros_like(mach)
            mask = mach >= self.threshold_mach
            term1 = -self.dissipation_coeff / (2.0 * mach[mask]**3)
            term2 = self.dissipation_coeff * self.threshold_mach**2 / mach[mask]**4
            result[mask] = term1 + term2
            return result
    
    def apply(self, residuals: np.ndarray, solution: Dict[str, np.ndarray], 
             grid: Dict[str, np.ndarray]) -> np.ndarray:
        """
        Apply artificial dissipation to the residuals.
        
        Args:
            residuals: Current residual values.
            solution: Current solution variables.
            grid: Computational grid.
            
        Returns:
            Modified residuals with artificial dissipation applied.
        """
        # Extract relevant variables from solution
        mach = solution['mach']
        density = solution['density']
        velocity = solution['velocity']
        
        # Compute dissipation coefficient based on local Mach number
        mu = self.compute_coefficient(mach)
        
        # Grid dimensions
        i_max, j_max = mach.shape
        
        # Initialize modified residuals
        modified_residuals = residuals.copy()
        
        # Apply dissipation as a bulk viscosity term
        for i in range(1, i_max-1):
            for j in range(1, j_max-1):
                if mu[i, j] > 0:
                    # Compute divergence of velocity field
                    div_v = (velocity[i+1, j, 0] - velocity[i-1, j, 0]) / (grid['x'][i+1, j] - grid['x'][i-1, j]) + \
                            (velocity[i, j+1, 1] - velocity[i, j-1, 1]) / (grid['y'][i, j+1] - grid['y'][i, j-1])
                    
                    # Scale by local density and dissipation coefficient
                    diss_term = mu[i, j] * density[i, j] * div_v
                    
                    # Add to residuals (for momentum and energy equations)
                    # Note: This is a simplified treatment and would need to be extended
                    # for the actual conservation form of the Euler equations
                    modified_residuals[1][i, j] += diss_term * velocity[i, j, 0]  # x-momentum
                    modified_residuals[2][i, j] += diss_term * velocity[i, j, 1]  # y-momentum
                    modified_residuals[3][i, j] += diss_term  # energy
        
        return modified_residuals


class LaplacianDissipation(ArtificialDissipation):
    """
    Laplacian-based artificial dissipation model.
    
    This model applies dissipation proportional to the Laplacian of the
    solution variables, with scaling based on local Mach number to 
    concentrate dissipation near shocks.
    
    Attributes:
        threshold_mach: Mach number threshold above which dissipation is applied.
        dissipation_coeff: Coefficient controlling dissipation strength.
        pressure_weight: Weight for pressure-based dissipation.
        entropy_weight: Weight for entropy-based dissipation.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the Laplacian dissipation model.
        
        Args:
            config: Configuration dictionary with the following keys:
                - threshold_mach: Mach number threshold (default: 0.95)
                - dissipation_coeff: Dissipation coefficient (default: 0.5)
                - pressure_weight: Weight for pressure-based dissipation (default: 1.0)
                - entropy_weight: Weight for entropy-based dissipation (default: 0.1)
        """
        super().__init__(config)
        
        self.threshold_mach = self.config.get('threshold_mach', 0.95)
        self.dissipation_coeff = self.config.get('dissipation_coeff', 0.5)
        self.pressure_weight = self.config.get('pressure_weight', 1.0)
        self.entropy_weight = self.config.get('entropy_weight', 0.1)
        
        logger.info(f"Initialized Laplacian dissipation with threshold_mach={self.threshold_mach}, "
                   f"dissipation_coeff={self.dissipation_coeff}")
    
    def compute_coefficient(self, mach: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Compute the artificial dissipation coefficient based on local Mach number.
        
        Args:
            mach: Local Mach number.
            
        Returns:
            Dissipation coefficient.
        """
        # Similar to bulk viscosity model, but can be tuned differently
        if np.isscalar(mach):
            if mach < self.threshold_mach:
                return 0.0
            else:
                return self.dissipation_coeff * (1.0 - (self.threshold_mach/mach)**2)
        else:
            result = np.zeros_like(mach)
            mask = mach >= self.threshold_mach
            result[mask] = self.dissipation_coeff * (1.0 - (self.threshold_mach/mach[mask])**2)
            return result
    
    def compute_jacobian(self, mach: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Compute the Jacobian of the artificial dissipation coefficient wrt Mach number.
        
        Args:
            mach: Local Mach number.
            
        Returns:
            Jacobian of dissipation coefficient wrt Mach number.
        """
        if np.isscalar(mach):
            if mach < self.threshold_mach:
                return 0.0
            else:
                return 2.0 * self.dissipation_coeff * self.threshold_mach**2 / mach**3
        else:
            result = np.zeros_like(mach)
            mask = mach >= self.threshold_mach
            result[mask] = 2.0 * self.dissipation_coeff * self.threshold_mach**2 / mach[mask]**3
            return result
    
    def apply(self, residuals: np.ndarray, solution: Dict[str, np.ndarray], 
             grid: Dict[str, np.ndarray]) -> np.ndarray:
        """
        Apply artificial dissipation to the residuals.
        
        Args:
            residuals: Current residual values.
            solution: Current solution variables.
            grid: Computational grid.
            
        Returns:
            Modified residuals with artificial dissipation applied.
        """
        # Extract relevant variables from solution
        mach = solution['mach']
        pressure = solution['pressure']
        density = solution['density']
        energy = solution.get('energy', None)
        
        # Compute dissipation coefficient based on local Mach number
        mu = self.compute_coefficient(mach)
        
        # Grid dimensions
        i_max, j_max = mach.shape
        
        # Initialize modified residuals
        modified_residuals = residuals.copy()
        
        # Apply dissipation based on pressure and entropy gradients
        for i in range(1, i_max-1):
            for j in range(1, j_max-1):
                if mu[i, j] > 0:
                    # Compute Laplacian of pressure (simplified, not accounting for grid metrics)
                    laplacian_p = (pressure[i+1, j] + pressure[i-1, j] + pressure[i, j+1] + pressure[i, j-1] - 
                                  4 * pressure[i, j])
                    
                    # Compute Laplacian of density (proxy for entropy variation)
                    laplacian_rho = (density[i+1, j] + density[i-1, j] + density[i, j+1] + density[i, j-1] - 
                                    4 * density[i, j])
                    
                    # Combined dissipation term
                    diss_term = mu[i, j] * (self.pressure_weight * laplacian_p / pressure[i, j] + 
                                           self.entropy_weight * laplacian_rho / density[i, j])
                    
                    # Add to residuals
                    modified_residuals[0][i, j] += diss_term  # continuity
                    modified_residuals[1][i, j] += diss_term  # x-momentum
                    modified_residuals[2][i, j] += diss_term  # y-momentum
                    
                    # Energy equation if available
                    if energy is not None and len(modified_residuals) > 3:
                        modified_residuals[3][i, j] += diss_term
        
        return modified_residuals


class AdaptiveDissipation(ArtificialDissipation):
    """
    Adaptive artificial dissipation model.
    
    This model dynamically adjusts dissipation based on local solution
    features such as pressure gradients and entropy changes, providing
    targeted dissipation only where needed.
    
    Attributes:
        base_coefficient: Base dissipation coefficient.
        shock_sensor_weight: Weight for shock sensor term.
        min_coefficient: Minimum dissipation coefficient.
        max_coefficient: Maximum dissipation coefficient.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the adaptive dissipation model.
        
        Args:
            config: Configuration dictionary with the following keys:
                - base_coefficient: Base dissipation coefficient (default: 0.3)
                - shock_sensor_weight: Weight for shock sensor (default: 5.0)
                - min_coefficient: Minimum dissipation coefficient (default: 0.0)
                - max_coefficient: Maximum dissipation coefficient (default: 1.0)
        """
        super().__init__(config)
        
        self.base_coefficient = self.config.get('base_coefficient', 0.3)
        self.shock_sensor_weight = self.config.get('shock_sensor_weight', 5.0)
        self.min_coefficient = self.config.get('min_coefficient', 0.0)
        self.max_coefficient = self.config.get('max_coefficient', 1.0)
        
        logger.info(f"Initialized adaptive dissipation with base_coefficient={self.base_coefficient}, "
                   f"shock_sensor_weight={self.shock_sensor_weight}")
    
    def compute_shock_sensor(self, pressure: np.ndarray) -> np.ndarray:
        """
        Compute a shock sensor based on pressure gradients.
        
        Args:
            pressure: Pressure field.
            
        Returns:
            Shock sensor field.
        """
        i_max, j_max = pressure.shape
        sensor = np.zeros_like(pressure)
        
        # Compute second derivative of pressure normalized by local pressure
        for i in range(1, i_max-1):
            for j in range(1, j_max-1):
                # Second derivative in x-direction
                d2p_dx2 = (pressure[i+1, j] - 2*pressure[i, j] + pressure[i-1, j])
                
                # Second derivative in y-direction
                d2p_dy2 = (pressure[i, j+1] - 2*pressure[i, j] + pressure[i, j-1])
                
                # Normalize by local pressure
                sensor[i, j] = abs(d2p_dx2 + d2p_dy2) / pressure[i, j]
        
        # Normalize sensor to [0, 1] range
        if np.max(sensor) > 0:
            sensor = sensor / np.max(sensor)
        
        return sensor
    
    def compute_coefficient(self, mach: Union[float, np.ndarray], 
                          shock_sensor: Optional[np.ndarray] = None) -> Union[float, np.ndarray]:
        """
        Compute the artificial dissipation coefficient based on local Mach number and shock sensor.
        
        Args:
            mach: Local Mach number.
            shock_sensor: Optional shock sensor field.
            
        Returns:
            Dissipation coefficient.
        """
        # Base dissipation scaled by Mach number
        if np.isscalar(mach):
            base_diss = self.base_coefficient * min(1.0, max(0.0, mach - 0.7) / 0.3)
            
            # If shock sensor is provided, use it to adjust dissipation
            if shock_sensor is not None:
                coef = base_diss + self.shock_sensor_weight * shock_sensor
            else:
                coef = base_diss
                
            # Clamp to range
            return max(self.min_coefficient, min(self.max_coefficient, coef))
        else:
            # Base dissipation for array
            base_diss = np.zeros_like(mach)
            mask = mach > 0.7
            base_diss[mask] = self.base_coefficient * np.minimum(1.0, np.maximum(0.0, mach[mask] - 0.7) / 0.3)
            
            # Adjust with shock sensor if provided
            if shock_sensor is not None:
                coef = base_diss + self.shock_sensor_weight * shock_sensor
            else:
                coef = base_diss
                
            # Clamp to range
            return np.maximum(self.min_coefficient, np.minimum(self.max_coefficient, coef))
    
    def compute_jacobian(self, mach: Union[float, np.ndarray],
                        shock_sensor: Optional[np.ndarray] = None) -> Union[float, np.ndarray]:
        """
        Compute the Jacobian of the artificial dissipation coefficient wrt Mach number.
        
        Args:
            mach: Local Mach number.
            shock_sensor: Optional shock sensor field.
            
        Returns:
            Jacobian of dissipation coefficient wrt Mach number.
        """
        # Derivative of base dissipation with respect to Mach number
        if np.isscalar(mach):
            if mach < 0.7 or mach > 1.0:
                return 0.0
            else:
                return self.base_coefficient / 0.3
        else:
            result = np.zeros_like(mach)
            mask = (mach >= 0.7) & (mach <= 1.0)
            result[mask] = self.base_coefficient / 0.3
            return result
    
    def apply(self, residuals: np.ndarray, solution: Dict[str, np.ndarray], 
             grid: Dict[str, np.ndarray]) -> np.ndarray:
        """
        Apply artificial dissipation to the residuals.
        
        Args:
            residuals: Current residual values.
            solution: Current solution variables.
            grid: Computational grid.
            
        Returns:
            Modified residuals with artificial dissipation applied.
        """
        # Extract relevant variables from solution
        mach = solution['mach']
        pressure = solution['pressure']
        
        # Compute shock sensor
        shock_sensor = self.compute_shock_sensor(pressure)
        
        # Compute dissipation coefficient
        mu = self.compute_coefficient(mach, shock_sensor)
        
        # Grid dimensions
        i_max, j_max = mach.shape
        
        # Initialize modified residuals
        modified_residuals = residuals.copy()
        
        # Apply dissipation as a weighted Laplacian
        # This is a simplified approach; actual implementation would
        # account for grid metrics and conservation form
        for var_idx in range(len(residuals)):
            var = solution.get(f'var_{var_idx}', None)
            if var is None:
                continue
                
            for i in range(1, i_max-1):
                for j in range(1, j_max-1):
                    # Compute Laplacian of variable
                    laplacian = (var[i+1, j] + var[i-1, j] + var[i, j+1] + var[i, j-1] - 4 * var[i, j])
                    
                    # Apply weighted dissipation
                    modified_residuals[var_idx][i, j] += mu[i, j] * laplacian
        
        return modified_residuals


def create_dissipation_model(model_type: str, config: Optional[Dict[str, Any]] = None) -> ArtificialDissipation:
    """
    Factory function to create an artificial dissipation model.
    
    Args:
        model_type: Type of dissipation model ('bulk_viscosity', 'laplacian', or 'adaptive').
        config: Configuration dictionary with model parameters.
        
    Returns:
        Instantiated artificial dissipation model.
        
    Raises:
        ValueError: If the specified model type is not supported.
    """
    model_map = {
        'bulk_viscosity': BulkViscosityDissipation,
        'laplacian': LaplacianDissipation,
        'adaptive': AdaptiveDissipation
    }
    
    if model_type not in model_map:
        raise ValueError(f"Unsupported dissipation model type: {model_type}. "
                        f"Supported types are: {list(model_map.keys())}")
    
    return model_map[model_type](config)
