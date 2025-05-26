"""
PyMISES - Periodicity Boundary Conditions

This module implements periodic boundary conditions for cascade simulations.
Enforces flow periodicity between upper and lower boundaries of the cascade domain.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Union, Callable

class PeriodicityBC:
    """Periodic boundary condition for cascade simulations."""
    
    def __init__(self, upper_indices: List[int], lower_indices: List[int], 
                 pitch: float = None, pitch_vector: List[float] = None):
        """
        Initialize periodicity boundary condition.
        
        Parameters
        ----------
        upper_indices : List[int]
            Indices of grid points on the upper periodic boundary
        lower_indices : List[int]
            Indices of grid points on the lower periodic boundary
        pitch : float, optional
            Blade-to-blade pitch (spacing) of the cascade
        pitch_vector : List[float], optional
            Vector representing the pitch direction and magnitude [dx, dy]
            If provided, this takes precedence over the pitch parameter
        """
        if len(upper_indices) != len(lower_indices):
            raise ValueError("upper_indices and lower_indices must have the same length")
        
        self.upper_indices = upper_indices
        self.lower_indices = lower_indices
        
        # Handle pitch specification
        if pitch_vector is not None:
            # Use pitch vector if provided
            self.pitch_vector = pitch_vector
            # Calculate scalar pitch as magnitude of the vector
            self.pitch = np.sqrt(pitch_vector[0]**2 + pitch_vector[1]**2)
        elif pitch is not None:
            # Use scalar pitch if provided
            self.pitch = pitch
            # Default pitch vector (vertical)
            self.pitch_vector = [0.0, self.pitch]
        else:
            raise ValueError("Either pitch or pitch_vector must be provided")
            
        self._validate_inputs()
    
    def _validate_inputs(self) -> None:
        """Validate input parameters."""
        if not isinstance(self.upper_indices, list) or not all(isinstance(i, int) for i in self.upper_indices):
            raise ValueError("upper_indices must be a list of integers")
        
        if not isinstance(self.lower_indices, list) or not all(isinstance(i, int) for i in self.lower_indices):
            raise ValueError("lower_indices must be a list of integers")
        
        if not isinstance(self.pitch, (int, float)) or self.pitch <= 0:
            raise ValueError("pitch must be a positive number")
            
        if not isinstance(self.pitch_vector, list) or len(self.pitch_vector) != 2:
            raise ValueError("pitch_vector must be a list of 2 numbers")
            
        if not all(isinstance(val, (int, float)) for val in self.pitch_vector):
            raise ValueError("pitch_vector must contain numeric values")
    
    def apply(self, solution: Dict) -> Dict:
        """
        Apply periodicity boundary condition to the solution.
        
        For cascade periodicity:
        - The flow at corresponding points on upper and lower boundaries should be identical
        - For simplicity, we average the values and apply to both boundaries
        
        Parameters
        ----------
        solution : Dict
            Current solution dictionary
            
        Returns
        -------
        Dict
            Updated solution dictionary
        """
        # Make a copy to avoid modifying the original solution
        updated_solution = solution.copy()
        
        # Get grid dimensions
        grid_size = len(solution['grid_x'])
        
        # Process each pair of periodic points
        for i, (upper_idx, lower_idx) in enumerate(zip(self.upper_indices, self.lower_indices)):
            # Skip if either index is out of bounds
            if upper_idx >= grid_size or lower_idx >= grid_size:
                continue
                
            # Average density
            avg_density = 0.5 * (solution['density'][upper_idx] + solution['density'][lower_idx])
            updated_solution['density'][upper_idx] = avg_density
            updated_solution['density'][lower_idx] = avg_density
            
            # Average velocities (need to account for possible direction flip in cascade)
            # This is a simplified approach; in a real implementation, we would need
            # to transform velocities based on the cascade stagger angle
            avg_vx = 0.5 * (solution['velocity_x'][upper_idx] + solution['velocity_x'][lower_idx])
            avg_vy = 0.5 * (solution['velocity_y'][upper_idx] + solution['velocity_y'][lower_idx])
            
            updated_solution['velocity_x'][upper_idx] = avg_vx
            updated_solution['velocity_x'][lower_idx] = avg_vx
            updated_solution['velocity_y'][upper_idx] = avg_vy
            updated_solution['velocity_y'][lower_idx] = avg_vy
            
            # Average pressure
            avg_pressure = 0.5 * (solution['pressure'][upper_idx] + solution['pressure'][lower_idx])
            updated_solution['pressure'][upper_idx] = avg_pressure
            updated_solution['pressure'][lower_idx] = avg_pressure
        
        return updated_solution
    
    def get_residual_contributions(self, solution: Dict) -> np.ndarray:
        """
        Get residual contributions from this boundary condition.
        
        Parameters
        ----------
        solution : Dict
            Current solution dictionary
            
        Returns
        -------
        np.ndarray
            Residual contributions, one for each pair of periodic points
        """
        # Residuals are the differences between corresponding points
        residuals = np.zeros(len(self.upper_indices) * 4)  # 4 variables per point
        
        for i, (upper_idx, lower_idx) in enumerate(zip(self.upper_indices, self.lower_indices)):
            # Density residual
            residuals[i*4] = solution['density'][upper_idx] - solution['density'][lower_idx]
            
            # Velocity residuals
            residuals[i*4+1] = solution['velocity_x'][upper_idx] - solution['velocity_x'][lower_idx]
            residuals[i*4+2] = solution['velocity_y'][upper_idx] - solution['velocity_y'][lower_idx]
            
            # Pressure residual
            residuals[i*4+3] = solution['pressure'][upper_idx] - solution['pressure'][lower_idx]
        
        return residuals
    
    def get_jacobian_contributions(self, solution: Dict) -> Tuple[List[Tuple], List[Tuple], List[float]]:
        """
        Get Jacobian contributions from this boundary condition.
        
        Parameters
        ----------
        solution : Dict
            Current solution dictionary
            
        Returns
        -------
        Tuple[List[Tuple], List[Tuple], List[float]]
            Row indices, column indices, and values for Jacobian matrix
        """
        rows = []
        cols = []
        values = []
        
        for i, (upper_idx, lower_idx) in enumerate(zip(self.upper_indices, self.lower_indices)):
            # Jacobian entries for density residual
            rows.append((i*4, 0))    # Row for density residual
            cols.append((upper_idx, 0))  # Column for upper density
            values.append(1.0)
            
            rows.append((i*4, 0))    # Row for density residual
            cols.append((lower_idx, 0))  # Column for lower density
            values.append(-1.0)
            
            # Jacobian entries for x-velocity residual
            rows.append((i*4+1, 0))    # Row for x-velocity residual
            cols.append((upper_idx, 1))  # Column for upper x-velocity
            values.append(1.0)
            
            rows.append((i*4+1, 0))    # Row for x-velocity residual
            cols.append((lower_idx, 1))  # Column for lower x-velocity
            values.append(-1.0)
            
            # Jacobian entries for y-velocity residual
            rows.append((i*4+2, 0))    # Row for y-velocity residual
            cols.append((upper_idx, 2))  # Column for upper y-velocity
            values.append(1.0)
            
            rows.append((i*4+2, 0))    # Row for y-velocity residual
            cols.append((lower_idx, 2))  # Column for lower y-velocity
            values.append(-1.0)
            
            # Jacobian entries for pressure residual
            rows.append((i*4+3, 0))    # Row for pressure residual
            cols.append((upper_idx, 3))  # Column for upper pressure
            values.append(1.0)
            
            rows.append((i*4+3, 0))    # Row for pressure residual
            cols.append((lower_idx, 3))  # Column for lower pressure
            values.append(-1.0)
        
        return rows, cols, values


class PhaseLagPeriodicityBC(PeriodicityBC):
    """
    Phase-lagged periodicity boundary condition for unsteady cascade simulations.
    
    This implements a time-shifted periodicity condition for unsteady simulations
    where a phase difference exists between adjacent blade passages.
    """
    
    def __init__(self, upper_indices: List[int], lower_indices: List[int], pitch: float,
                 interblade_phase_angle: float, time_period: float, time_history: Dict):
        """
        Initialize phase-lagged periodicity boundary condition.
        
        Parameters
        ----------
        upper_indices : List[int]
            Indices of grid points on the upper periodic boundary
        lower_indices : List[int]
            Indices of grid points on the lower periodic boundary
        pitch : float
            Blade-to-blade pitch (spacing) of the cascade
        interblade_phase_angle : float
            Interblade phase angle in radians
        time_period : float
            Time period of the unsteady flow
        time_history : Dict
            Dictionary containing time history of the solution at the boundaries
        """
        super().__init__(upper_indices, lower_indices, pitch)
        self.interblade_phase_angle = interblade_phase_angle
        self.time_period = time_period
        self.time_history = time_history
    
    def apply(self, solution: Dict, current_time: float) -> Dict:
        """
        Apply phase-lagged periodicity boundary condition to the solution.
        
        Parameters
        ----------
        solution : Dict
            Current solution dictionary
        current_time : float
            Current simulation time
            
        Returns
        -------
        Dict
            Updated solution dictionary
        """
        # Make a copy to avoid modifying the original solution
        updated_solution = solution.copy()
        
        # Get grid dimensions
        grid_size = len(solution['grid_x'])
        
        # Calculate phase time shift
        time_shift = self.interblade_phase_angle / (2 * np.pi) * self.time_period
        
        # Phase-shifted times for upper and lower boundaries
        upper_time = (current_time - time_shift) % self.time_period
        lower_time = (current_time + time_shift) % self.time_period
        
        # Update boundaries using time-shifted data from history
        # This is a simplified approach; a real implementation would need
        # to interpolate in time and handle storage of time history data
        
        # Store current boundary values in time history
        self._store_current_values(solution, current_time)
        
        # Apply time-shifted values from history
        self._apply_time_shifted_values(updated_solution, upper_time, lower_time)
        
        return updated_solution
    
    def _store_current_values(self, solution: Dict, current_time: float) -> None:
        """
        Store current boundary values in time history.
        
        Parameters
        ----------
        solution : Dict
            Current solution dictionary
        current_time : float
            Current simulation time
        """
        # Round time to nearest time step in the history
        time_index = int(round(current_time / self.time_period * len(self.time_history['times']))) % len(self.time_history['times'])
        
        # Get grid dimensions
        grid_size = len(solution['grid_x'])
        
        # Store current values for upper boundary
        for i, idx in enumerate(self.upper_indices):
            # Skip if index is out of bounds
            if idx >= grid_size:
                continue
                
            self.time_history['upper_density'][time_index][i] = solution['density'][idx]
            self.time_history['upper_velocity_x'][time_index][i] = solution['velocity_x'][idx]
            self.time_history['upper_velocity_y'][time_index][i] = solution['velocity_y'][idx]
            self.time_history['upper_pressure'][time_index][i] = solution['pressure'][idx]
        
        # Store current values for lower boundary
        for i, idx in enumerate(self.lower_indices):
            # Skip if index is out of bounds
            if idx >= grid_size:
                continue
                
            self.time_history['lower_density'][time_index][i] = solution['density'][idx]
            self.time_history['lower_velocity_x'][time_index][i] = solution['velocity_x'][idx]
            self.time_history['lower_velocity_y'][time_index][i] = solution['velocity_y'][idx]
            self.time_history['lower_pressure'][time_index][i] = solution['pressure'][idx]
    
    def _apply_time_shifted_values(self, solution: Dict, upper_time: float, lower_time: float) -> None:
        """
        Apply time-shifted values from history to the solution.
        
        Parameters
        ----------
        solution : Dict
            Solution dictionary to update
        upper_time : float
            Time for upper boundary data
        lower_time : float
            Time for lower boundary data
        """
        # Find nearest time indices in history
        upper_time_index = int(round(upper_time / self.time_period * len(self.time_history['times']))) % len(self.time_history['times'])
        lower_time_index = int(round(lower_time / self.time_period * len(self.time_history['times']))) % len(self.time_history['times'])
        
        # Get grid dimensions
        grid_size = len(solution['grid_x'])
        
        # Apply values to upper boundary from lower boundary at shifted time
        for i, idx in enumerate(self.upper_indices):
            # Skip if index is out of bounds
            if idx >= grid_size:
                continue
                
            solution['density'][idx] = self.time_history['lower_density'][lower_time_index][i]
            solution['velocity_x'][idx] = self.time_history['lower_velocity_x'][lower_time_index][i]
            solution['velocity_y'][idx] = self.time_history['lower_velocity_y'][lower_time_index][i]
            solution['pressure'][idx] = self.time_history['lower_pressure'][lower_time_index][i]
        
        # Apply values to lower boundary from upper boundary at shifted time
        for i, idx in enumerate(self.lower_indices):
            # Skip if index is out of bounds
            if idx >= grid_size:
                continue
                
            solution['density'][idx] = self.time_history['upper_density'][upper_time_index][i]
            solution['velocity_x'][idx] = self.time_history['upper_velocity_x'][upper_time_index][i]
            solution['velocity_y'][idx] = self.time_history['upper_velocity_y'][upper_time_index][i]
            solution['pressure'][idx] = self.time_history['upper_pressure'][upper_time_index][i]