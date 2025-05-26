"""
Improved artificial dissipation models for PyMISES.

This module provides enhanced implementations of artificial dissipation methods
used in the Euler solver to stabilize solutions in regions with shocks.
Key improvements:
- Matrix dissipation for better shock capturing
- Entropy-based sensors for detecting physical discontinuities
- Automatic coefficient scaling based on local flow conditions
- Better handling of boundaries and corner regions
"""

import numpy as np
from typing import Dict, Tuple, Union, Optional, Callable, Any
import scipy.sparse as sp

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
        self.thermo = None
        # Initialize thermodynamics if provided in the config
        if 'gamma' in self.config and 'gas_constant' in self.config:
            self.thermo = Thermodynamics(
                self.config.get('gamma', 1.4),
                self.config.get('gas_constant', 287.058)
            )
    
    def compute_coefficient(self, mach: Union[float, np.ndarray], 
                         pressure: Optional[Union[float, np.ndarray]] = None,
                         density: Optional[Union[float, np.ndarray]] = None) -> Union[float, np.ndarray]:
        """
        Compute the artificial dissipation coefficient.
        
        Args:
            mach: Local Mach number.
            pressure: Optional local pressure for pressure-based switches.
            density: Optional local density for density-based switches.
            
        Returns:
            Dissipation coefficient.
        """
        raise NotImplementedError("Subclasses must implement compute_coefficient")
    
    def compute_jacobian(self, mach: Union[float, np.ndarray],
                       pressure: Optional[Union[float, np.ndarray]] = None,
                       density: Optional[Union[float, np.ndarray]] = None) -> Union[float, np.ndarray]:
        """
        Compute the Jacobian of the artificial dissipation coefficient wrt Mach number.
        
        Args:
            mach: Local Mach number.
            pressure: Optional local pressure for pressure-based switches.
            density: Optional local density for density-based switches.
            
        Returns:
            Jacobian of dissipation coefficient wrt Mach number.
        """
        raise NotImplementedError("Subclasses must implement compute_jacobian")
    
    def apply_dissipation(self, residuals: np.ndarray, 
                        solution: Dict[str, np.ndarray], 
                        grid: Dict[str, np.ndarray],
                        grid_metrics: Dict[str, np.ndarray]) -> np.ndarray:
        """
        Apply artificial dissipation to the residuals.
        
        Args:
            residuals: Current residual values.
            solution: Current solution variables.
            grid: Computational grid.
            grid_metrics: Grid metrics for coordinate transformations.
            
        Returns:
            Modified residuals with artificial dissipation applied.
        """
        raise NotImplementedError("Subclasses must implement apply_dissipation")


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
        
        # Additional parameters for robustness
        self.min_coefficient = self.config.get('min_coefficient', 0.01)
        self.max_coefficient = self.config.get('max_coefficient', 2.0)
        
        logger.info(f"Initialized bulk viscosity dissipation with threshold_mach={self.threshold_mach}, "
                   f"dissipation_coeff={self.dissipation_coeff}")
    
    def compute_coefficient(self, mach: Union[float, np.ndarray],
                          pressure: Optional[Union[float, np.ndarray]] = None,
                          density: Optional[Union[float, np.ndarray]] = None) -> Union[float, np.ndarray]:
        """
        Compute the artificial dissipation coefficient based on local Mach number.
        
        Args:
            mach: Local Mach number.
            pressure: Optional local pressure (unused in this model).
            density: Optional local density (unused in this model).
            
        Returns:
            Dissipation coefficient.
        """
        # Only apply dissipation in regions where Mach number exceeds the threshold
        if np.isscalar(mach):
            if mach < self.threshold_mach:
                return self.min_coefficient
            else:
                # Smooth activation around threshold
                activation = min(1.0, (mach - self.threshold_mach) / 0.1)
                return min(self.max_coefficient, 
                          self.min_coefficient + activation * self.dissipation_coeff * (1.0 - (self.threshold_mach/mach)**2) / (2.0 * mach**2))
        else:
            # Array version
            result = np.ones_like(mach) * self.min_coefficient
            mask = mach >= self.threshold_mach
            
            if np.any(mask):
                # Smooth activation around threshold
                activation = np.minimum(1.0, (mach[mask] - self.threshold_mach) / 0.1)
                result[mask] = np.minimum(
                    self.max_coefficient, 
                    self.min_coefficient + activation * self.dissipation_coeff * (1.0 - (self.threshold_mach/mach[mask])**2) / (2.0 * mach[mask]**2)
                )
            
            return result
    
    def compute_jacobian(self, mach: Union[float, np.ndarray],
                         pressure: Optional[Union[float, np.ndarray]] = None,
                         density: Optional[Union[float, np.ndarray]] = None) -> Union[float, np.ndarray]:
        """
        Compute the Jacobian of the artificial dissipation coefficient wrt Mach number.
        
        Args:
            mach: Local Mach number.
            pressure: Optional local pressure (unused in this model).
            density: Optional local density (unused in this model).
            
        Returns:
            Jacobian of dissipation coefficient wrt Mach number.
        """
        # Derivative of the dissipation coefficient with respect to Mach number
        if np.isscalar(mach):
            if mach < self.threshold_mach - 0.1:
                return 0.0
            elif mach < self.threshold_mach:
                # Smooth activation region
                return self.dissipation_coeff * 10.0  # Derivative of activation ramp
            else:
                # Check if at max coefficient
                if self.min_coefficient + self.dissipation_coeff * (1.0 - (self.threshold_mach/mach)**2) / (2.0 * mach**2) > self.max_coefficient:
                    return 0.0
                
                term1 = -self.dissipation_coeff / (2.0 * mach**3)
                term2 = self.dissipation_coeff * self.threshold_mach**2 / mach**4
                return term1 + term2
        else:
            # Array version
            result = np.zeros_like(mach)
            
            # Activation region
            activation_region = (mach >= self.threshold_mach - 0.1) & (mach < self.threshold_mach)
            result[activation_region] = self.dissipation_coeff * 10.0
            
            # Above threshold region
            above_threshold = mach >= self.threshold_mach
            
            if np.any(above_threshold):
                # Check for max coefficient
                coeff_values = self.min_coefficient + self.dissipation_coeff * (1.0 - (self.threshold_mach/mach[above_threshold])**2) / (2.0 * mach[above_threshold]**2)
                not_maxed = coeff_values < self.max_coefficient
                
                # Only calculate derivatives for points not at max
                if np.any(not_maxed):
                    mach_subset = mach[above_threshold][not_maxed]
                    term1 = -self.dissipation_coeff / (2.0 * mach_subset**3)
                    term2 = self.dissipation_coeff * self.threshold_mach**2 / mach_subset**4
                    
                    # Place values back in result array at correct indices
                    indices = np.where(above_threshold)[0][not_maxed]
                    result[indices] = term1 + term2
            
            return result
    
    def apply_dissipation(self, residuals: np.ndarray, 
                        solution: Dict[str, np.ndarray], 
                        grid: Dict[str, np.ndarray],
                        grid_metrics: Dict[str, np.ndarray]) -> np.ndarray:
        """
        Apply artificial dissipation to the residuals.
        
        Args:
            residuals: Current residual values, shape (n_equations, ni, nj).
            solution: Current solution variables.
            grid: Computational grid.
            grid_metrics: Grid metrics for coordinate transformations.
            
        Returns:
            Modified residuals with artificial dissipation applied.
        """
        # Extract relevant variables from solution
        mach = solution.get('mach')
        if mach is None:
            # Calculate Mach number if not provided
            if 'velocity' in solution and 'pressure' in solution and 'density' in solution:
                velocity = solution['velocity']
                pressure = solution['pressure']
                density = solution['density']
                
                # Calculate velocity magnitude
                if len(velocity.shape) == 3:  # (ni, nj, 2)
                    v_mag = np.sqrt(velocity[:, :, 0]**2 + velocity[:, :, 1]**2)
                else:  # Separate components
                    v_mag = np.sqrt(solution.get('velocity_x', 0)**2 + solution.get('velocity_y', 0)**2)
                
                # Calculate sound speed
                if self.thermo is None:
                    # Default values if thermo not initialized
                    gamma = 1.4
                    R = 287.058
                else:
                    gamma = self.thermo.gamma
                    R = self.thermo.R
                
                temperature = pressure / (density * R)
                sound_speed = np.sqrt(gamma * R * temperature)
                
                # Calculate Mach number
                mach = v_mag / sound_speed
            else:
                logger.warning("Cannot compute Mach number for dissipation calculation")
                return residuals
        
        # Compute dissipation coefficient based on local Mach number
        mu = self.compute_coefficient(mach)
        
        # Grid dimensions
        ni, nj = mach.shape
        
        # Number of equations
        n_equations = len(residuals)
        
        # Initialize modified residuals
        modified_residuals = residuals.copy()
        
        # Get spacing
        ds = grid_metrics.get('ds', np.ones((ni, nj)))
        dn = grid_metrics.get('dn', np.ones((ni, nj)))
        
        # Apply bulk viscosity dissipation to each equation
        for eq in range(n_equations):
            residual = residuals[eq]
            
            for i in range(1, ni-1):
                for j in range(1, nj-1):
                    if mu[i, j] > 0:
                        # Get the relevant variable for this equation
                        var = solution.get(f'var_{eq}', None)
                        if var is None:
                            # Try to get standard variables based on equation index
                            if eq == 0:
                                var = solution.get('density', None)
                            elif eq == 1:
                                var = solution.get('velocity_x', None)
                            elif eq == 2:
                                var = solution.get('velocity_y', None)
                            elif eq == 3:
                                var = solution.get('pressure', None)
                            else:
                                continue
                        
                        if var is None:
                            continue
                        
                        # Compute Laplacian of the variable
                        # Streamwise second derivative
                        d2var_ds2 = (var[i+1, j] - 2*var[i, j] + var[i-1, j]) / ds[i, j]**2
                        
                        # Normal second derivative
                        d2var_dn2 = (var[i, j+1] - 2*var[i, j] + var[i, j-1]) / dn[i, j]**2
                        
                        # Total Laplacian
                        laplacian = d2var_ds2 + d2var_dn2
                        
                        # Apply dissipation
                        modified_residuals[eq][i, j] += mu[i, j] * laplacian
        
        return modified_residuals


class MatrixDissipation(ArtificialDissipation):
    """
    Matrix-based artificial dissipation model.
    
    This model applies dissipation using a matrix formulation that respects
    the characteristic structure of the Euler equations, providing better
    properties for shock capturing and stability.
    
    Attributes:
        threshold_mach: Mach number threshold above which dissipation is applied.
        k2: Coefficient for second-order dissipation (shock regions).
        k4: Coefficient for fourth-order dissipation (smooth regions).
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the matrix dissipation model.
        
        Args:
            config: Configuration dictionary with the following keys:
                - threshold_mach: Mach number threshold (default: 0.95)
                - k2: Second-order dissipation coefficient (default: 0.5)
                - k4: Fourth-order dissipation coefficient (default: 0.016)
        """
        super().__init__(config)
        
        self.threshold_mach = self.config.get('threshold_mach', 0.95)
        self.k2 = self.config.get('k2', 0.5)
        self.k4 = self.config.get('k4', 0.016)
        
        # Additional parameters for entropy-based switching
        self.use_entropy_fix = self.config.get('use_entropy_fix', True)
        self.entropy_fix_coeff = self.config.get('entropy_fix_coeff', 0.1)
        
        # Parameters for eigenvalue scaling
        self.eigenvalue_scale = self.config.get('eigenvalue_scale', 1.0)
        
        logger.info(f"Initialized matrix dissipation with k2={self.k2}, k4={self.k4}, "
                  f"entropy_fix={self.use_entropy_fix}")
    
    def compute_coefficient(self, mach: Union[float, np.ndarray],
                          pressure: Optional[Union[float, np.ndarray]] = None,
                          density: Optional[Union[float, np.ndarray]] = None) -> Union[float, np.ndarray]:
        """
        Compute the artificial dissipation coefficient based on local flow conditions.
        
        For matrix dissipation, this is a pressure-based switch that detects shocks.
        
        Args:
            mach: Local Mach number.
            pressure: Local pressure field.
            density: Local density field (unused in this model).
            
        Returns:
            Pressure-based shock sensor.
        """
        if pressure is None:
            # Return Mach-based coefficient if pressure not provided
            if np.isscalar(mach):
                return 0.0 if mach < self.threshold_mach else self.k2
            else:
                result = np.zeros_like(mach)
                result[mach >= self.threshold_mach] = self.k2
                return result
        
        # Calculate pressure-based shock sensor
        ni, nj = pressure.shape
        sensor = np.zeros_like(pressure)
        
        # Compute pressure-based sensor using Jameson's approach
        for i in range(1, ni-1):
            for j in range(1, nj-1):
                # Calculate normalized second derivative of pressure
                p_ij = pressure[i, j]
                p_im1j = pressure[i-1, j]
                p_ip1j = pressure[i+1, j]
                p_ijm1 = pressure[i, j-1]
                p_ijp1 = pressure[i, j+1]
                
                # Compute denominators with safety factor
                denom_i = 0.5 * (p_ip1j + p_im1j) + 1e-10
                denom_j = 0.5 * (p_ijp1 + p_ijm1) + 1e-10
                
                # Compute normalized pressure differences
                nu_i = abs(p_ip1j - 2*p_ij + p_im1j) / denom_i
                nu_j = abs(p_ijp1 - 2*p_ij + p_ijm1) / denom_j
                
                # Take maximum as sensor value
                sensor[i, j] = self.k2 * max(nu_i, nu_j)
        
        # Apply Mach number threshold to prevent excessive dissipation in low Mach regions
        if isinstance(mach, np.ndarray) and mach.shape == sensor.shape:
            sensor[mach < self.threshold_mach] *= 0.1  # Reduce dissipation in subsonic regions
        
        return sensor
    
    def compute_jacobian(self, mach: Union[float, np.ndarray],
                       pressure: Optional[Union[float, np.ndarray]] = None, 
                       density: Optional[Union[float, np.ndarray]] = None) -> Union[float, np.ndarray]:
        """
        Compute the Jacobian of the artificial dissipation coefficient.
        
        For matrix dissipation with pressure-based switching, this is very complex.
        This implementation provides a simplified version.
        
        Args:
            mach: Local Mach number.
            pressure: Local pressure field.
            density: Local density field (unused in this model).
            
        Returns:
            Simplified Jacobian of the shock sensor wrt Mach and pressure.
        """
        # This is a simplified version that just returns non-zero values where the
        # sensor is active to ensure the Jacobian structure is correct
        if pressure is None:
            # Return Mach-based Jacobian if pressure not provided
            if np.isscalar(mach):
                return 0.0 if mach < self.threshold_mach else 0.01
            else:
                result = np.zeros_like(mach)
                result[mach >= self.threshold_mach] = 0.01
                return result
        
        # Return a simplified Jacobian based on sensor
        sensor = self.compute_coefficient(mach, pressure)
        return 0.01 * (sensor > 0.01)
    
    def compute_eigenvalues(self, solution: Dict[str, np.ndarray], 
                          grid_metrics: Dict[str, np.ndarray], 
                          i: int, j: int, direction: str) -> np.ndarray:
        """
        Compute the eigenvalues of the flux Jacobian matrix.
        
        Args:
            solution: Current solution variables.
            grid_metrics: Grid metrics for coordinate transformations.
            i: i-index.
            j: j-index.
            direction: Direction ('i' or 'j').
            
        Returns:
            Array of eigenvalues.
        """
        # Extract solution variables at the given point
        density = solution['density'][i, j]
        pressure = solution['pressure'][i, j]
        if 'velocity' in solution and len(solution['velocity'].shape) == 3:
            velocity = solution['velocity'][i, j]
            u = velocity[0]
            v = velocity[1]
        else:
            u = solution['velocity_x'][i, j]
            v = solution['velocity_y'][i, j]
        
        # Extract grid metrics
        if direction == 'i':
            # Streamwise direction
            metrics = grid_metrics['s_metrics'][i, j]
            nx, ny = metrics
        else:
            # Normal direction
            metrics = grid_metrics['n_metrics'][i, j]
            nx, ny = metrics
        
        # Normalize metrics
        metrics_norm = np.sqrt(nx**2 + ny**2)
        if metrics_norm > 0:
            nx /= metrics_norm
            ny /= metrics_norm
        
        # Calculate the contravariant velocity
        contravariant_velocity = u * nx + v * ny
        
        # Calculate sound speed
        if self.thermo is None:
            # Default values if thermo not initialized
            gamma = 1.4
        else:
            gamma = self.thermo.gamma
        
        sound_speed = np.sqrt(gamma * pressure / density)
        
        # Compute the eigenvalues
        eigenvalues = np.array([
            contravariant_velocity,
            contravariant_velocity,
            contravariant_velocity + sound_speed,
            contravariant_velocity - sound_speed
        ])
        
        # Apply entropy fix to prevent vanishing eigenvalues
        if self.use_entropy_fix:
            # Minimum eigenvalue magnitude
            eps = self.entropy_fix_coeff * sound_speed
            
            # Apply entropy fix
            for i in range(len(eigenvalues)):
                if abs(eigenvalues[i]) < eps:
                    eigenvalues[i] = np.sign(eigenvalues[i]) * eps
        
        return eigenvalues * self.eigenvalue_scale
    
    def apply_dissipation(self, residuals: np.ndarray, 
                        solution: Dict[str, np.ndarray], 
                        grid: Dict[str, np.ndarray],
                        grid_metrics: Dict[str, np.ndarray]) -> np.ndarray:
        """
        Apply matrix-based artificial dissipation to the residuals.
        
        Args:
            residuals: Current residual values.
            solution: Current solution variables.
            grid: Computational grid.
            grid_metrics: Grid metrics for coordinate transformations.
            
        Returns:
            Modified residuals with artificial dissipation applied.
        """
        # Extract solution variables
        mach = solution.get('mach')
        pressure = solution.get('pressure')
        
        if mach is None or pressure is None:
            logger.warning("Cannot compute matrix dissipation without mach and pressure fields")
            return residuals
        
        # Compute shock sensor
        sensor = self.compute_coefficient(mach, pressure)
        
        # Grid dimensions
        ni, nj = mach.shape
        
        # Number of equations
        n_equations = len(residuals)
        
        # Initialize modified residuals
        modified_residuals = residuals.copy()
        
        # Get conservative variables
        u_conserv = self._solution_to_conservative(solution, n_equations)
        
        # Apply dissipation in i-direction (streamwise)
        for i in range(2, ni-2):
            for j in range(2, nj-2):
                # Compute spectral radius (maximum eigenvalue)
                lambda_i = np.max(np.abs(self.compute_eigenvalues(solution, grid_metrics, i, j, 'i')))
                
                # Pressure-based sensors for switching between 2nd/4th order
                eps2_i = sensor[i, j]
                eps4_i = max(0, self.k4 - eps2_i)
                
                # Apply to each equation
                for eq in range(n_equations):
                    # Second-order dissipation (captures shocks)
                    d2_i = eps2_i * lambda_i * (u_conserv[i+1, j, eq] - u_conserv[i, j, eq] - 
                                               u_conserv[i, j, eq] + u_conserv[i-1, j, eq])
                    
                    # Fourth-order dissipation (background dissipation for stability)
                    d4_i = eps4_i * lambda_i * (u_conserv[i+2, j, eq] - 3*u_conserv[i+1, j, eq] + 
                                               3*u_conserv[i, j, eq] - u_conserv[i-1, j, eq])
                    
                    # Add dissipation to residual
                    modified_residuals[eq][i, j] += d2_i - d4_i
        
        # Apply dissipation in j-direction (normal)
        for i in range(2, ni-2):
            for j in range(2, nj-2):
                # Compute spectral radius
                lambda_j = np.max(np.abs(self.compute_eigenvalues(solution, grid_metrics, i, j, 'j')))
                
                # Pressure-based sensors
                eps2_j = sensor[i, j]
                eps4_j = max(0, self.k4 - eps2_j)
                
                # Apply to each equation
                for eq in range(n_equations):
                    # Second-order dissipation
                    d2_j = eps2_j * lambda_j * (u_conserv[i, j+1, eq] - u_conserv[i, j, eq] - 
                                               u_conserv[i, j, eq] + u_conserv[i, j-1, eq])
                    
                    # Fourth-order dissipation
                    d4_j = eps4_j * lambda_j * (u_conserv[i, j+2, eq] - 3*u_conserv[i, j+1, eq] + 
                                               3*u_conserv[i, j, eq] - u_conserv[i, j-1, eq])
                    
                    # Add dissipation to residual
                    modified_residuals[eq][i, j] += d2_j - d4_j
        
        return modified_residuals
    
    def _solution_to_conservative(self, solution: Dict[str, np.ndarray], n_equations: int) -> np.ndarray:
        """
        Convert solution variables to conservative variables for dissipation.
        
        Args:
            solution: Current solution variables.
            n_equations: Number of equations in the system.
            
        Returns:
            Array of conservative variables.
        """
        # Extract grid dimensions
        if 'density' in solution:
            ni, nj = solution['density'].shape
        elif 'pressure' in solution:
            ni, nj = solution['pressure'].shape
        else:
            logger.error("Cannot determine grid dimensions from solution variables")
            return np.array([])
        
        # Initialize conservative variables array
        u_conserv = np.zeros((ni, nj, n_equations))
        
        # Extract primitive variables
        density = solution['density']
        
        if 'velocity' in solution and len(solution['velocity'].shape) == 3:
            velocity = solution['velocity']
            u = velocity[:, :, 0]
            v = velocity[:, :, 1]
        else:
            u = solution['velocity_x']
            v = solution['velocity_y']
        
        pressure = solution['pressure']
        
        # Fill in conservative variables
        if n_equations >= 1:
            u_conserv[:, :, 0] = density  # Continuity
        
        if n_equations >= 3:
            u_conserv[:, :, 1] = density * u  # x-momentum
            u_conserv[:, :, 2] = density * v  # y-momentum
        
        if n_equations >= 4:
            # Energy
            if 'energy' in solution:
                u_conserv[:, :, 3] = solution['energy']
            else:
                # Compute energy from primitive variables
                gamma = 1.4 if self.thermo is None else self.thermo.gamma
                u_conserv[:, :, 3] = pressure / (gamma - 1) + 0.5 * density * (u**2 + v**2)
        
        return u_conserv
