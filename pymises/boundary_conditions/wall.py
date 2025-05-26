"""
PyMISES - Wall Boundary Conditions

This module implements solid wall boundary conditions for the Euler and coupled solver.
Includes both inviscid (slip) and viscous (with displacement thickness) wall conditions.
"""

import numpy as np
import abc
from typing import Dict, List, Tuple, Optional, Union, Callable

class WallBoundaryCondition(abc.ABC):
    """Base class for wall boundary conditions."""
    
    def __init__(self, grid_indices: List[int], normal_direction: str = "inner"):
        """
        Initialize wall boundary condition.
        
        Parameters
        ----------
        grid_indices : List[int]
            Indices of grid points on the wall boundary
        normal_direction : str
            Direction of the wall normal ('inner' or 'outer')
        """
        self.grid_indices = grid_indices
        self.normal_direction = normal_direction
        self._validate_inputs()
    
    def _validate_inputs(self) -> None:
        """Validate input parameters."""
        if not isinstance(self.grid_indices, list) or not all(isinstance(i, int) for i in self.grid_indices):
            raise ValueError("grid_indices must be a list of integers")
        
        if self.normal_direction not in ["inner", "outer"]:
            raise ValueError("normal_direction must be 'inner' or 'outer'")
    
    @abc.abstractmethod
    def apply(self, solution: Dict) -> Dict:
        """
        Apply boundary condition to the solution.
        
        Parameters
        ----------
        solution : Dict
            Current solution dictionary
            
        Returns
        -------
        Dict
            Updated solution dictionary
        """
        pass
    
    @abc.abstractmethod
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
            Residual contributions
        """
        pass
    
    @abc.abstractmethod
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
        pass


class InviscidWallBC(WallBoundaryCondition):
    """Inviscid (slip) wall boundary condition."""
    
    def __init__(self, grid_indices: List[int], normal_direction: str = "inner"):
        """
        Initialize inviscid wall boundary condition.
        
        Parameters
        ----------
        grid_indices : List[int]
            Indices of grid points on the wall boundary
        normal_direction : str
            Direction of the wall normal ('inner' or 'outer')
        """
        super().__init__(grid_indices, normal_direction)
    
    def apply(self, solution: Dict) -> Dict:
        """
        Apply inviscid wall boundary condition to the solution.
        Ensures that the velocity vector is tangent to the wall.
        
        Parameters
        ----------
        solution : Dict
            Current solution dictionary with keys:
            - 'density': Density field
            - 'velocity_x': X-component of velocity
            - 'velocity_y': Y-component of velocity
            - 'pressure': Pressure field
            - 'grid_x': X-coordinates of grid points
            - 'grid_y': Y-coordinates of grid points
            
        Returns
        -------
        Dict
            Updated solution dictionary
        """
        # Make a copy to avoid modifying the original solution
        updated_solution = solution.copy()
        
        # Extract grid point coordinates
        grid_x = solution['grid_x']
        grid_y = solution['grid_y']
        
        for i in self.grid_indices:
            # Calculate wall tangent vector
            # We need to use neighboring points to compute the tangent
            # For simplicity, we'll use central differencing if possible
            if i > 0 and i < len(grid_x) - 1:
                dx = grid_x[i+1] - grid_x[i-1]
                dy = grid_y[i+1] - grid_y[i-1]
            elif i == 0:
                dx = grid_x[i+1] - grid_x[i]
                dy = grid_y[i+1] - grid_y[i]
            else:  # i == len(grid_x) - 1
                dx = grid_x[i] - grid_x[i-1]
                dy = grid_y[i] - grid_y[i-1]
            
            # Normalize tangent vector
            magnitude = np.sqrt(dx**2 + dy**2)
            tx = dx / magnitude
            ty = dy / magnitude
            
            # Normal vector (perpendicular to tangent)
            if self.normal_direction == "inner":
                nx = -ty
                ny = tx
            else:  # "outer"
                nx = ty
                ny = -tx
            
            # Extract current velocity
            vx = solution['velocity_x'][i]
            vy = solution['velocity_y'][i]
            
            # Calculate dot product of velocity and normal vector
            vn = vx * nx + vy * ny
            
            # Remove normal component to ensure tangential flow
            updated_solution['velocity_x'][i] = vx - vn * nx
            updated_solution['velocity_y'][i] = vy - vn * ny
        
        return updated_solution
    
    def get_residual_contributions(self, solution: Dict) -> np.ndarray:
        """
        Get residual contributions from this boundary condition.
        For inviscid wall, we enforce zero normal velocity.
        
        Parameters
        ----------
        solution : Dict
            Current solution dictionary
            
        Returns
        -------
        np.ndarray
            Residual contributions, one for each boundary point
        """
        # Extract grid point coordinates
        grid_x = solution['grid_x']
        grid_y = solution['grid_y']
        
        # Initialize residuals
        residuals = np.zeros(len(self.grid_indices))
        
        for idx, i in enumerate(self.grid_indices):
            # Calculate wall normal vector using same method as in apply()
            if i > 0 and i < len(grid_x) - 1:
                dx = grid_x[i+1] - grid_x[i-1]
                dy = grid_y[i+1] - grid_y[i-1]
            elif i == 0:
                dx = grid_x[i+1] - grid_x[i]
                dy = grid_y[i+1] - grid_y[i]
            else:  # i == len(grid_x) - 1
                dx = grid_x[i] - grid_x[i-1]
                dy = grid_y[i] - grid_y[i-1]
            
            # Normalize tangent vector
            magnitude = np.sqrt(dx**2 + dy**2)
            tx = dx / magnitude
            ty = dy / magnitude
            
            # Normal vector (perpendicular to tangent)
            if self.normal_direction == "inner":
                nx = -ty
                ny = tx
            else:  # "outer"
                nx = ty
                ny = -tx
            
            # Extract current velocity
            vx = solution['velocity_x'][i]
            vy = solution['velocity_y'][i]
            
            # Calculate dot product of velocity and normal vector
            # This should be zero for a perfect wall BC
            vn = vx * nx + vy * ny
            
            # Residual is the normal velocity component
            residuals[idx] = vn
        
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
        # Extract grid point coordinates
        grid_x = solution['grid_x']
        grid_y = solution['grid_y']
        
        # Initialize lists for sparse matrix construction
        rows = []
        cols = []
        values = []
        
        for idx, i in enumerate(self.grid_indices):
            # Calculate wall normal vector (same as in other methods)
            if i > 0 and i < len(grid_x) - 1:
                dx = grid_x[i+1] - grid_x[i-1]
                dy = grid_y[i+1] - grid_y[i-1]
            elif i == 0:
                dx = grid_x[i+1] - grid_x[i]
                dy = grid_y[i+1] - grid_y[i]
            else:  # i == len(grid_x) - 1
                dx = grid_x[i] - grid_x[i-1]
                dy = grid_y[i] - grid_y[i-1]
            
            # Normalize tangent vector
            magnitude = np.sqrt(dx**2 + dy**2)
            tx = dx / magnitude
            ty = dy / magnitude
            
            # Normal vector (perpendicular to tangent)
            if self.normal_direction == "inner":
                nx = -ty
                ny = tx
            else:  # "outer"
                nx = ty
                ny = -tx
            
            # Derivatives of residual (normal velocity) with respect to velocities
            # δ(v·n)/δvx = nx
            # δ(v·n)/δvy = ny
            rows.append((i, 0))  # Row for this boundary condition residual
            cols.append((i, 0))  # Column for vx at this point
            values.append(nx)
            
            rows.append((i, 0))  # Row for this boundary condition residual
            cols.append((i, 1))  # Column for vy at this point
            values.append(ny)
            
            # Derivatives with respect to grid positions would also be included
            # for a fully general implementation, but omitted for simplicity
        
        return rows, cols, values


class ViscousWallBC(WallBoundaryCondition):
    """Viscous wall boundary condition with displacement thickness effect."""
    
    def __init__(self, grid_indices: List[int], normal_direction: str = "inner",
                 displacement_thickness_provider: Callable[[int], float] = None):
        """
        Initialize viscous wall boundary condition.
        
        Parameters
        ----------
        grid_indices : List[int]
            Indices of grid points on the wall boundary
        normal_direction : str
            Direction of the wall normal ('inner' or 'outer')
        displacement_thickness_provider : Callable[[int], float], optional
            Function that returns displacement thickness for a given grid index
        """
        super().__init__(grid_indices, normal_direction)
        self.displacement_thickness_provider = displacement_thickness_provider or (lambda i: 0.0)
    
    def apply(self, solution: Dict) -> Dict:
        """
        Apply viscous wall boundary condition to the solution.
        Accounts for displacement thickness by imposing transpiration velocity.
        
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
        
        # Extract grid point coordinates
        grid_x = solution['grid_x']
        grid_y = solution['grid_y']
        
        for i in self.grid_indices:
            # Calculate wall tangent and normal vectors (similar to InviscidWallBC)
            if i > 0 and i < len(grid_x) - 1:
                dx = grid_x[i+1] - grid_x[i-1]
                dy = grid_y[i+1] - grid_y[i-1]
            elif i == 0:
                dx = grid_x[i+1] - grid_x[i]
                dy = grid_y[i+1] - grid_y[i]
            else:  # i == len(grid_x) - 1
                dx = grid_x[i] - grid_x[i-1]
                dy = grid_y[i] - grid_y[i-1]
            
            # Normalize tangent vector
            magnitude = np.sqrt(dx**2 + dy**2)
            tx = dx / magnitude
            ty = dy / magnitude
            
            # Normal vector (perpendicular to tangent)
            if self.normal_direction == "inner":
                nx = -ty
                ny = tx
            else:  # "outer"
                nx = ty
                ny = -tx
            
            # Extract current velocity
            vx = solution['velocity_x'][i]
            vy = solution['velocity_y'][i]
            
            # Calculate edge velocity (tangential component)
            v_edge = vx * tx + vy * ty
            
            # Get displacement thickness at this point
            delta_star = self.displacement_thickness_provider(i)
            
            # Calculate transpiration velocity
            # Formula: v_n = d(δ* * v_edge)/ds = v_edge * dδ*/ds + δ* * dv_edge/ds
            # For simplicity, we'll approximate the derivatives
            if i > 0 and i < len(self.grid_indices) - 1:
                delta_star_prev = self.displacement_thickness_provider(i-1)
                delta_star_next = self.displacement_thickness_provider(i+1)
                
                vx_prev = solution['velocity_x'][i-1]
                vy_prev = solution['velocity_y'][i-1]
                v_edge_prev = vx_prev * tx + vy_prev * ty
                
                vx_next = solution['velocity_x'][i+1]
                vy_next = solution['velocity_y'][i+1]
                v_edge_next = vx_next * tx + vy_next * ty
                
                ddelta_star_ds = (delta_star_next - delta_star_prev) / (2 * magnitude)
                dv_edge_ds = (v_edge_next - v_edge_prev) / (2 * magnitude)
            elif i == 0:
                delta_star_next = self.displacement_thickness_provider(i+1)
                
                vx_next = solution['velocity_x'][i+1]
                vy_next = solution['velocity_y'][i+1]
                v_edge_next = vx_next * tx + vy_next * ty
                
                ddelta_star_ds = (delta_star_next - delta_star) / magnitude
                dv_edge_ds = (v_edge_next - v_edge) / magnitude
            else:  # i == len(self.grid_indices) - 1
                delta_star_prev = self.displacement_thickness_provider(i-1)
                
                vx_prev = solution['velocity_x'][i-1]
                vy_prev = solution['velocity_y'][i-1]
                v_edge_prev = vx_prev * tx + vy_prev * ty
                
                ddelta_star_ds = (delta_star - delta_star_prev) / magnitude
                dv_edge_ds = (v_edge - v_edge_prev) / magnitude
            
            # Calculate transpiration velocity
            v_transpiration = v_edge * ddelta_star_ds + delta_star * dv_edge_ds
            
            # Apply transpiration velocity to modify the normal component
            updated_solution['velocity_x'][i] = vx - (vx * nx + vy * ny - v_transpiration) * nx
            updated_solution['velocity_y'][i] = vy - (vx * nx + vy * ny - v_transpiration) * ny
        
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
            Residual contributions
        """
        # Extract grid point coordinates
        grid_x = solution['grid_x']
        grid_y = solution['grid_y']
        
        # Initialize residuals
        residuals = np.zeros(len(self.grid_indices))
        
        for idx, i in enumerate(self.grid_indices):
            # Calculate wall tangent and normal vectors
            if i > 0 and i < len(grid_x) - 1:
                dx = grid_x[i+1] - grid_x[i-1]
                dy = grid_y[i+1] - grid_y[i-1]
            elif i == 0:
                dx = grid_x[i+1] - grid_x[i]
                dy = grid_y[i+1] - grid_y[i]
            else:  # i == len(grid_x) - 1
                dx = grid_x[i] - grid_x[i-1]
                dy = grid_y[i] - grid_y[i-1]
            
            # Normalize tangent vector
            magnitude = np.sqrt(dx**2 + dy**2)
            tx = dx / magnitude
            ty = dy / magnitude
            
            # Normal vector (perpendicular to tangent)
            if self.normal_direction == "inner":
                nx = -ty
                ny = tx
            else:  # "outer"
                nx = ty
                ny = -tx
            
            # Extract current velocity
            vx = solution['velocity_x'][i]
            vy = solution['velocity_y'][i]
            
            # Calculate edge velocity (tangential component)
            v_edge = vx * tx + vy * ty
            
            # Get displacement thickness at this point
            delta_star = self.displacement_thickness_provider(i)
            
            # Calculate transpiration velocity (simplified from apply method)
            if i > 0 and i < len(self.grid_indices) - 1:
                delta_star_prev = self.displacement_thickness_provider(i-1)
                delta_star_next = self.displacement_thickness_provider(i+1)
                
                vx_prev = solution['velocity_x'][i-1]
                vy_prev = solution['velocity_y'][i-1]
                v_edge_prev = vx_prev * tx + vy_prev * ty
                
                vx_next = solution['velocity_x'][i+1]
                vy_next = solution['velocity_y'][i+1]
                v_edge_next = vx_next * tx + vy_next * ty
                
                ddelta_star_ds = (delta_star_next - delta_star_prev) / (2 * magnitude)
                dv_edge_ds = (v_edge_next - v_edge_prev) / (2 * magnitude)
            elif i == 0:
                delta_star_next = self.displacement_thickness_provider(i+1)
                
                vx_next = solution['velocity_x'][i+1]
                vy_next = solution['velocity_y'][i+1]
                v_edge_next = vx_next * tx + vy_next * ty
                
                ddelta_star_ds = (delta_star_next - delta_star) / magnitude
                dv_edge_ds = (v_edge_next - v_edge) / magnitude
            else:  # i == len(self.grid_indices) - 1
                delta_star_prev = self.displacement_thickness_provider(i-1)
                
                vx_prev = solution['velocity_x'][i-1]
                vy_prev = solution['velocity_y'][i-1]
                v_edge_prev = vx_prev * tx + vy_prev * ty
                
                ddelta_star_ds = (delta_star - delta_star_prev) / magnitude
                dv_edge_ds = (v_edge - v_edge_prev) / magnitude
            
            # Calculate transpiration velocity
            v_transpiration = v_edge * ddelta_star_ds + delta_star * dv_edge_ds
            
            # Normal velocity component
            v_normal = vx * nx + vy * ny
            
            # Residual is the difference between normal velocity and transpiration velocity
            residuals[idx] = v_normal - v_transpiration
        
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
        # Extract grid point coordinates
        grid_x = solution['grid_x']
        grid_y = solution['grid_y']
        
        # Initialize lists for sparse matrix construction
        rows = []
        cols = []
        values = []
        
        for idx, i in enumerate(self.grid_indices):
            # Calculate wall tangent and normal vectors
            if i > 0 and i < len(grid_x) - 1:
                dx = grid_x[i+1] - grid_x[i-1]
                dy = grid_y[i+1] - grid_y[i-1]
            elif i == 0:
                dx = grid_x[i+1] - grid_x[i]
                dy = grid_y[i+1] - grid_y[i]
            else:  # i == len(grid_x) - 1
                dx = grid_x[i] - grid_x[i-1]
                dy = grid_y[i] - grid_y[i-1]
            
            # Normalize tangent vector
            magnitude = np.sqrt(dx**2 + dy**2)
            tx = dx / magnitude
            ty = dy / magnitude
            
            # Normal vector (perpendicular to tangent)
            if self.normal_direction == "inner":
                nx = -ty
                ny = tx
            else:  # "outer"
                nx = ty
                ny = -tx
            
            # Direct contribution of velocity to residual
            # δ(v·n - v_transpiration)/δvx = nx - derivatives of transpiration term
            # δ(v·n - v_transpiration)/δvy = ny - derivatives of transpiration term
            # For simplicity, we'll include only the direct nx, ny contributions
            rows.append((i, 0))  # Row for this boundary condition residual
            cols.append((i, 0))  # Column for vx at this point
            values.append(nx)
            
            rows.append((i, 0))  # Row for this boundary condition residual
            cols.append((i, 1))  # Column for vy at this point
            values.append(ny)
            
            # For a complete implementation, we would also include the derivatives 
            # of transpiration velocity with respect to both local and neighboring
            # velocities, as well as displacement thickness dependencies
        
        return rows, cols, values