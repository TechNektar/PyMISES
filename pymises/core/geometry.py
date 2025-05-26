"""
Airfoil/cascade geometry handling module for PyMISES.

This module provides classes and functions for managing and manipulating
airfoil geometry and cascade configurations.
"""

import numpy as np
import os
from typing import Dict, Tuple, List, Union, Optional, Any, Callable

from pymises.utils.logger import get_logger
from pymises.utils.validation import validate_airfoil_coordinates

logger = get_logger(__name__)

class BladeGeometry:
    """
    Base class for airfoil/blade geometries.
    
    This class handles the representation and manipulation of airfoil
    or turbine blade geometries, providing methods for initialization,
    modification, and export.
    
    Attributes:
        x: x-coordinates of the geometry points.
        y: y-coordinates of the geometry points.
        name: Name of the geometry.
        n_points: Number of points defining the geometry.
        trailing_edge_idx: Index of the trailing edge point.
        is_closed: Boolean indicating if the geometry is closed.
    """
    
    @classmethod
    def create_flat_plate(cls, length: float = 1.0, n_points: int = 101) -> 'AirfoilGeometry':
        """
        Create a flat plate airfoil.
        
        Args:
            length: Length of the flat plate.
            n_points: Number of points to generate.
            
        Returns:
            AirfoilGeometry object representing a flat plate.
        """
        # Create equidistant points along the x-axis
        x = np.linspace(0, length, n_points)
        y = np.zeros_like(x)  # All y=0 for a flat plate
        
        # Create the airfoil geometry using the coordinates
        coords = np.column_stack((x, y))
        return cls(coords=coords, name="Flat Plate")
    
    @classmethod
    def create_naca(cls, code: str = '0012', n_points: int = 100, closed_te: bool = True) -> 'BladeGeometry':
        """
        Create a blade geometry based on NACA 4-digit airfoil.
        
        Args:
            code: NACA 4-digit code (e.g., '0012')
            n_points: Number of points to generate
            closed_te: Flag indicating if trailing edge should be closed
            
        Returns:
            BladeGeometry object with NACA airfoil shape
        """
        # Use the NACA4DigitAirfoil class to generate coordinates
        naca_airfoil = NACA4DigitAirfoil(code, n_points, closed_te)
        
        # Create a new BladeGeometry object with these coordinates
        return cls(coords=naca_airfoil.get_coordinates(), name=f"NACA{code}")
    
    def __init__(self, coords: Optional[Union[np.ndarray, List]] = None, 
                params: Optional[Dict[str, Any]] = None,
                name: str = "unnamed"):
        """
        Initialize blade geometry.
        
        Args:
            coords: Array of shape (n,2) containing [x,y] coordinates,
                   or separate x and y arrays if second param is an array,
                   or None to generate from parameters.
            params: Dictionary of parameters for generating the geometry,
                   or y-coordinates if coords is x-coordinates,
                   or None to use provided coordinates.
            name: Name of the geometry.
        """
        self.name = name
        
        # Handle different ways to provide coordinates
        if coords is not None and isinstance(params, (np.ndarray, list)):
            # Case: BladeGeometry(x, y)
            x = np.asarray(coords)
            y = np.asarray(params)
            self.x = x
            self.y = y
            self._params = {}
        elif coords is not None:
            # Case: BladeGeometry(coords) where coords is an array of shape (n,2)
            # Initialize from coordinates
            coords = np.asarray(coords)
            if coords.ndim != 2 or coords.shape[1] != 2:
                raise ValueError("Coordinates must be a Nx2 array")
            
            self.x = coords[:, 0]
            self.y = coords[:, 1]
            # Fix: Explicitly check if params is None instead of using boolean evaluation
            self._params = {} if params is None else params
        elif params is not None:
            # Generate geometry from parameters
            self._params = params
            self.generate_from_params()
        else:
            # No initialization data provided
            self.x = np.array([])
            self.y = np.array([])
            self._params = {}
        
        # Validate coordinates
        if len(self.x) > 0:
            try:
                validate_airfoil_coordinates(self.x, self.y)
            except ValueError as e:
                logger.warning(f"Geometry validation issue: {str(e)}")
            
            # Ensure the airfoil is closed if it has enough points
            if len(self.x) >= 3 and not self.is_closed:
                logger.info("Automatically closing airfoil geometry")
                self.ensure_closed()
    
    @property
    def n_points(self) -> int:
        """Get the number of points in the geometry."""
        return len(self.x)
    
    @property
    def is_closed(self) -> bool:
        """Check if the geometry is closed (first and last points are the same)."""
        if self.n_points < 2:
            return False
        return np.isclose(self.x[0], self.x[-1]) and np.isclose(self.y[0], self.y[-1])
    
    @property
    def trailing_edge_idx(self) -> int:
        """Get the index of the trailing edge point."""
        if self.n_points == 0:
            return -1
        
        # For closed airfoils, the trailing edge is typically at the end
        if self.is_closed:
            return self.n_points - 1
        
        # Otherwise, find the point with the maximum x value
        return np.argmax(self.x)
    
    def generate_from_params(self) -> None:
        """
        Generate geometry from parameters.
        
        This method should be implemented by subclasses to generate geometry
        points from parameter dictionaries.
        """
        raise NotImplementedError("Subclasses must implement generate_from_params")
    
    def get_coordinates(self) -> np.ndarray:
        """
        Get discretized (x,y) coordinates.
        
        Returns:
            Array of shape (n,2) containing [x,y] coordinates.
        """
        return np.column_stack((self.x, self.y))
    
    def get_arclength(self) -> np.ndarray:
        """
        Calculate arclength distribution along the geometry.
        
        Returns:
            Array of arclength values, starting from 0 at the first point.
        """
        if self.n_points < 2:
            return np.array([0.0] if self.n_points == 1 else [])
        
        # Calculate segment lengths
        dx = np.diff(self.x)
        dy = np.diff(self.y)
        segment_lengths = np.sqrt(dx**2 + dy**2)
        
        # Compute cumulative arclength
        arclength = np.zeros(self.n_points)
        arclength[1:] = np.cumsum(segment_lengths)
        
        return arclength
    
    def get_curvature(self) -> np.ndarray:
        """
        Calculate curvature distribution along the geometry.
        
        Returns:
            Array of curvature values at each point.
        """
        if self.n_points < 3:
            return np.zeros(self.n_points)
        
        # Get arclength distribution
        s = self.get_arclength()
        
        # Calculate derivatives with respect to arclength
        # Use central differences for interior points
        dx_ds = np.zeros(self.n_points)
        dy_ds = np.zeros(self.n_points)
        
        # Interior points (central difference)
        for i in range(1, self.n_points-1):
            if s[i+1] != s[i-1]:  # Avoid division by zero
                dx_ds[i] = (self.x[i+1] - self.x[i-1]) / (s[i+1] - s[i-1])
                dy_ds[i] = (self.y[i+1] - self.y[i-1]) / (s[i+1] - s[i-1])
        
        # End points (forward and backward differences)
        if s[1] != s[0]:
            dx_ds[0] = (self.x[1] - self.x[0]) / (s[1] - s[0])
            dy_ds[0] = (self.y[1] - self.y[0]) / (s[1] - s[0])
        
        if s[-1] != s[-2]:
            dx_ds[-1] = (self.x[-1] - self.x[-2]) / (s[-1] - s[-2])
            dy_ds[-1] = (self.y[-1] - self.y[-2]) / (s[-1] - s[-2])
        
        # Second derivatives
        d2x_ds2 = np.zeros(self.n_points)
        d2y_ds2 = np.zeros(self.n_points)
        
        # Interior points (central difference)
        for i in range(1, self.n_points-1):
            if s[i+1] != s[i-1]:  # Avoid division by zero
                d2x_ds2[i] = (dx_ds[i+1] - dx_ds[i-1]) / (s[i+1] - s[i-1])
                d2y_ds2[i] = (dy_ds[i+1] - dy_ds[i-1]) / (s[i+1] - s[i-1])
        
        # Curvature formula: κ = (x'y'' - y'x'') / (x'^2 + y'^2)^(3/2)
        denominator = (dx_ds**2 + dy_ds**2)**(1.5)
        curvature = np.zeros(self.n_points)
        
        # Avoid division by zero
        nonzero_denom = denominator > 1e-10
        if np.any(nonzero_denom):
            curvature[nonzero_denom] = (dx_ds[nonzero_denom] * d2y_ds2[nonzero_denom] - 
                                       dy_ds[nonzero_denom] * d2x_ds2[nonzero_denom]) / denominator[nonzero_denom]
        
        return curvature
    
    def get_normal_vectors(self) -> np.ndarray:
        """
        Calculate normal vectors at each point on the geometry.
        
        Returns:
            Array of shape (n,2) containing [nx,ny] normal vectors.
        """
        if self.n_points < 2:
            return np.zeros((self.n_points, 2))
        
        # Get arclength distribution
        s = self.get_arclength()
        
        # Calculate tangent vectors (derivatives with respect to arclength)
        dx_ds = np.zeros(self.n_points)
        dy_ds = np.zeros(self.n_points)
        
        # Interior points (central difference)
        for i in range(1, self.n_points-1):
            if s[i+1] != s[i-1]:  # Avoid division by zero
                dx_ds[i] = (self.x[i+1] - self.x[i-1]) / (s[i+1] - s[i-1])
                dy_ds[i] = (self.y[i+1] - self.y[i-1]) / (s[i+1] - s[i-1])
        
        # End points (forward and backward differences)
        if s[1] != s[0]:
            dx_ds[0] = (self.x[1] - self.x[0]) / (s[1] - s[0])
            dy_ds[0] = (self.y[1] - self.y[0]) / (s[1] - s[0])
        
        if s[-1] != s[-2]:
            dx_ds[-1] = (self.x[-1] - self.x[-2]) / (s[-1] - s[-2])
            dy_ds[-1] = (self.y[-1] - self.y[-2]) / (s[-1] - s[-2])
        
        # Normalize tangent vectors
        magnitude = np.sqrt(dx_ds**2 + dy_ds**2)
        nonzero_mag = magnitude > 1e-10
        
        if np.any(nonzero_mag):
            dx_ds[nonzero_mag] /= magnitude[nonzero_mag]
            dy_ds[nonzero_mag] /= magnitude[nonzero_mag]
        
        # Normal vectors are perpendicular to tangent vectors
        normal_vectors = np.zeros((self.n_points, 2))
        normal_vectors[:, 0] = -dy_ds  # nx = -dy/ds
        normal_vectors[:, 1] = dx_ds   # ny = dx/ds
        
        return normal_vectors
    
    def ensure_closed(self) -> None:
        """
        Ensure that the geometry is closed by making the first and last points identical.
        
        This method checks if the first and last points are the same and, if not,
        adds a duplicate of the first point at the end to close the geometry.
        """
        if self.n_points < 3:
            logger.warning("Cannot close geometry with less than 3 points")
            return
        
        # Check if already closed
        if self.is_closed:
            return
        
        # If not closed, add a copy of the first point at the end
        self.x = np.append(self.x, self.x[0])
        self.y = np.append(self.y, self.y[0])
        
        logger.debug(f"Closed geometry by adding duplicate point, now has {self.n_points} points")
    
    def redistribute_points(self, n_points: int = None, 
                          clustering: Optional[Dict[str, Any]] = None) -> None:
        """
        Redistribute points along the geometry with specified clustering.
        
        Args:
            n_points: New number of points (if None, keeps the current number).
            clustering: Dictionary with clustering parameters:
                - le_clustering: Leading edge clustering factor (0-1)
                - te_clustering: Trailing edge clustering factor (0-1)
                - method: Clustering method ('sine', 'tanh', 'exp')
        """
        if self.n_points < 2:
            logger.warning("Cannot redistribute points on geometry with less than 2 points")
            return
        
        # Default parameters
        n_new = n_points if n_points is not None else self.n_points
        cluster_params = clustering or {}
        
        le_clustering = cluster_params.get('le_clustering', 0.2)
        te_clustering = cluster_params.get('te_clustering', 0.3)
        method = cluster_params.get('method', 'sine')
        
        # Get current arclength distribution
        arclength = self.get_arclength()
        total_length = arclength[-1]
        
        # Create new arclength distribution with clustering
        if method == 'sine':
            # Sine-based clustering (more points at each end)
            t = np.linspace(0, np.pi, n_new)
            s_new = total_length * (1 - np.cos(t)) / 2
        
        elif method == 'tanh':
            # Hyperbolic tangent clustering
            beta = 5.0  # Clustering strength
            alpha_le = le_clustering * beta
            alpha_te = te_clustering * beta
            
            t = np.linspace(-1, 1, n_new)
            s_new = total_length * (1 + np.tanh(alpha_le * t) / np.tanh(alpha_le)) / 2
            
        elif method == 'exp':
            # Exponential clustering
            beta_le = -np.log(le_clustering) if le_clustering > 0 else 3.0
            beta_te = -np.log(te_clustering) if te_clustering > 0 else 3.0
            
            t = np.linspace(0, 1, n_new)
            scale_factor = 1 / (1 - np.exp(-beta_le) - np.exp(-beta_te) + np.exp(-beta_le-beta_te))
            s_new = total_length * scale_factor * (1 - np.exp(-beta_le * t) - np.exp(-beta_te * (1-t)) + 
                                                 np.exp(-beta_le * t - beta_te * (1-t)))
        
        else:
            # Linear distribution (no clustering)
            s_new = np.linspace(0, total_length, n_new)
        
        # Interpolate coordinates at new arclength positions
        x_new = np.interp(s_new, arclength, self.x)
        y_new = np.interp(s_new, arclength, self.y)
        
        # Update geometry
        self.x = x_new
        self.y = y_new
    
    def smooth(self, smoothing_factor: float = 0.5, n_iterations: int = 1) -> None:
        """
        Apply smoothing to the geometry.
        
        Args:
            smoothing_factor: Weight for averaging (0-1, higher means more smoothing).
            n_iterations: Number of smoothing iterations.
        """
        if self.n_points < 3:
            logger.warning("Cannot smooth geometry with less than 3 points")
            return
        
        # Save endpoints
        x_start, y_start = self.x[0], self.y[0]
        x_end, y_end = self.x[-1], self.y[-1]
        
        # Multiple smoothing iterations
        for _ in range(n_iterations):
            # Apply smoothing to interior points
            x_smooth = self.x.copy()
            y_smooth = self.y.copy()
            
            for i in range(1, self.n_points-1):
                x_smooth[i] = (1 - smoothing_factor) * self.x[i] + smoothing_factor * (self.x[i-1] + self.x[i+1]) / 2
                y_smooth[i] = (1 - smoothing_factor) * self.y[i] + smoothing_factor * (self.y[i-1] + self.y[i+1]) / 2
            
            self.x = x_smooth
            self.y = y_smooth
        
        # Restore endpoints if needed
        if self.is_closed:
            self.x[0] = self.x[-1] = x_start
            self.y[0] = self.y[-1] = y_start
        else:
            self.x[0], self.y[0] = x_start, y_start
            self.x[-1], self.y[-1] = x_end, y_end
    
    def scale(self, scale_factor: float = 1.0, reference_point: Optional[Tuple[float, float]] = None) -> None:
        """
        Scale the geometry.
        
        Args:
            scale_factor: Scaling factor to apply.
            reference_point: Point to scale around (default is the centroid).
        """
        if self.n_points == 0:
            return
        
        # Default reference point is the centroid
        if reference_point is None:
            ref_x = np.mean(self.x)
            ref_y = np.mean(self.y)
        else:
            ref_x, ref_y = reference_point
        
        # Apply scaling
        self.x = ref_x + scale_factor * (self.x - ref_x)
        self.y = ref_y + scale_factor * (self.y - ref_y)
    
    def get_surface_indices(self) -> np.ndarray:
        """
        Get the indices of points on the blade surface.
        
        For a simple blade geometry, this is all points.
        
        Returns:
            Array of indices representing surface points.
        """
        return np.arange(self.n_points)
    
    def rotate(self, angle_deg: float, reference_point: Optional[Tuple[float, float]] = None) -> None:
        """
        Rotate the geometry.
        
        Args:
            angle_deg: Rotation angle in degrees (positive is counterclockwise).
            reference_point: Point to rotate around (default is the centroid).
        """
        if self.n_points == 0:
            return
        
        # Default reference point is the centroid
        if reference_point is None:
            ref_x = np.mean(self.x)
            ref_y = np.mean(self.y)
        else:
            ref_x, ref_y = reference_point
        
        # Convert angle to radians
        angle_rad = np.radians(angle_deg)
        cos_theta = np.cos(angle_rad)
        sin_theta = np.sin(angle_rad)
        
        # Translate to reference point
        x_translated = self.x - ref_x
        y_translated = self.y - ref_y
        
        # Rotate
        x_rotated = x_translated * cos_theta - y_translated * sin_theta
        y_rotated = x_translated * sin_theta + y_translated * cos_theta
        
        # Translate back
        self.x = x_rotated + ref_x
        self.y = y_rotated + ref_y
    
    def translate(self, dx: float, dy: float) -> None:
        """
        Translate the geometry.
        
        Args:
            dx: Translation in x-direction.
            dy: Translation in y-direction.
        """
        if self.n_points == 0:
            return
        
        self.x += dx
        self.y += dy
    
    def deform(self, deformation_function: Callable[[np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray]]) -> None:
        """
        Apply a custom deformation function to the geometry.
        
        Args:
            deformation_function: Function that takes (x,y) arrays and returns new (x,y) arrays.
        """
        if self.n_points == 0:
            return
        
        # Apply the deformation function
        new_x, new_y = deformation_function(self.x, self.y)
        
        # Update the geometry
        if len(new_x) == len(self.x) and len(new_y) == len(self.y):
            self.x = new_x
            self.y = new_y
        else:
            logger.warning("Deformation function returned arrays of different length - ignoring")
    
    def export_to_file(self, filename: str, format: str = 'dat') -> bool:
        """
        Export geometry to a file.
        
        Args:
            filename: Output filename.
            format: File format ('dat', 'csv', etc.).
            
        Returns:
            True if export was successful, False otherwise.
        """
        if self.n_points == 0:
            logger.warning("Cannot export empty geometry")
            return False
        
        try:
            # Ensure the directory exists
            os.makedirs(os.path.dirname(os.path.abspath(filename)), exist_ok=True)
            
            # Write coordinates to file
            with open(filename, 'w') as f:
                # Write header based on format
                if format.lower() == 'csv':
                    f.write("x,y\n")
                elif format.lower() == 'dat':
                    f.write(f"{self.name}\n{self.n_points}\n")
                
                # Write coordinates
                for i in range(self.n_points):
                    if format.lower() == 'csv':
                        f.write(f"{self.x[i]},{self.y[i]}\n")
                    else:
                        f.write(f"{self.x[i]} {self.y[i]}\n")
            
            logger.info(f"Exported geometry to {filename}")
            return True
            
        except Exception as e:
            logger.error(f"Error exporting geometry: {str(e)}")
            return False
    
    @classmethod
    def import_from_file(cls, filename: str, format: str = None) -> 'BladeGeometry':
        """
        Import geometry from a file.
        
        Args:
            filename: Input filename.
            format: File format (if None, inferred from file extension).
            
        Returns:
            BladeGeometry object.
        """
        # Infer format from file extension if not provided
        if format is None:
            _, ext = os.path.splitext(filename)
            format = ext[1:]  # Remove leading dot
        
        try:
            # Read coordinates from file
            coords = []
            name = os.path.basename(filename)
            
            with open(filename, 'r') as f:
                lines = f.readlines()
                
                if format.lower() == 'csv':
                    # Skip header if present
                    start_line = 1 if ',' in lines[0].lower() and 'x' in lines[0].lower() else 0
                    
                    for line in lines[start_line:]:
                        if ',' in line:
                            parts = line.strip().split(',')
                            if len(parts) >= 2:
                                try:
                                    x = float(parts[0])
                                    y = float(parts[1])
                                    coords.append((x, y))
                                except ValueError:
                                    continue
                
                elif format.lower() in ['dat', 'txt']:
                    # First line might be name
                    name = lines[0].strip()
                    
                    # Second line might be number of points
                    try:
                        n_points = int(lines[1].strip())
                        start_line = 2
                    except ValueError:
                        start_line = 1
                    
                    for line in lines[start_line:]:
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            try:
                                x = float(parts[0])
                                y = float(parts[1])
                                coords.append((x, y))
                            except ValueError:
                                continue
            
            # Create geometry object
            if coords:
                return cls(np.array(coords), name=name)
            else:
                logger.warning(f"No valid coordinates found in {filename}")
                return cls(name=name)
            
        except Exception as e:
            logger.error(f"Error importing geometry: {str(e)}")
            return cls(name=os.path.basename(filename))


class NACA4DigitAirfoil(BladeGeometry):
    """
    NACA 4-digit series airfoil.
    
    This class implements the NACA 4-digit series airfoil generation method.
    The parameters dictionary should include:
    - 'code': NACA 4-digit code (e.g., '0012')
    - 'n_points': Number of points to generate
    - 'closed_te': Boolean indicating if trailing edge should be closed
    
    Attributes:
        code: NACA 4-digit code.
        closed_te: Flag indicating if trailing edge is closed.
    """
    
    def __init__(self, code: str = '0012', n_points: int = 100, closed_te: bool = True, 
                coords: Optional[np.ndarray] = None, name: str = None):
        """
        Initialize NACA 4-digit airfoil.
        
        Args:
            code: NACA 4-digit code (e.g., '0012').
            n_points: Number of points to generate.
            closed_te: Flag indicating if trailing edge should be closed.
            coords: Array of shape (n,2) containing [x,y] coordinates (if provided, other params are ignored).
            name: Airfoil name (default is "NACA" + code).
        """
        self.code = code
        self.closed_te = closed_te
        
        # Default name is NACA + code
        if name is None:
            name = f"NACA{code}"
        
        # Create parameter dictionary
        params = {
            'code': code,
            'n_points': n_points,
            'closed_te': closed_te
        }
        
        # Initialize base class
        super().__init__(coords=coords, params=params, name=name)
    
    def generate_from_params(self) -> None:
        """Generate NACA 4-digit airfoil geometry from parameters."""
        # Extract parameters
        code = self._params.get('code', '0012')
        n_points = self._params.get('n_points', 100)
        closed_te = self._params.get('closed_te', True)
        
        # Parse NACA code
        if len(code) != 4:
            logger.warning(f"Invalid NACA code: {code}, using default 0012")
            code = '0012'
        
        try:
            m = int(code[0]) / 100.0  # Maximum camber
            p = int(code[1]) / 10.0   # Location of maximum camber
            t = int(code[2:]) / 100.0  # Thickness
        except ValueError:
            logger.warning(f"Invalid NACA code: {code}, using default 0012")
            m = 0.0
            p = 0.0
            t = 0.12
        
        # Generate cosine-spaced points for better resolution near leading edge
        beta = np.linspace(0, np.pi, n_points)
        x = 0.5 * (1 - np.cos(beta))  # Cosine spacing from 0 to 1
        
        # Calculate thickness distribution
        y_t = t/0.2 * (0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)
        
        # Adjust trailing edge if closed
        if closed_te:
            y_t[-1] = 0.0
        
        # Calculate camber line and its slope
        y_c = np.zeros_like(x)
        dyc_dx = np.zeros_like(x)
        
        if m > 0:  # Cambered airfoil
            # Front part (x < p)
            mask = x <= p
            if p > 0:  # Avoid division by zero
                y_c[mask] = m/p**2 * (2*p*x[mask] - x[mask]**2)
                dyc_dx[mask] = 2*m/p**2 * (p - x[mask])
            
            # Rear part (x >= p)
            mask = x > p
            if p < 1:  # Avoid division by zero
                y_c[mask] = m/(1-p)**2 * (1 - 2*p + 2*p*x[mask] - x[mask]**2)
                dyc_dx[mask] = 2*m/(1-p)**2 * (p - x[mask])
        
        # Calculate the angle of the camber line
        theta = np.arctan(dyc_dx)
        
        # Calculate upper and lower surface coordinates
        x_upper = x - y_t * np.sin(theta)
        y_upper = y_c + y_t * np.cos(theta)
        
        x_lower = x + y_t * np.sin(theta)
        y_lower = y_c - y_t * np.cos(theta)
        
        # Combine to form the airfoil (counterclockwise from trailing edge)
        self.x = np.concatenate((x_upper[::-1], x_lower[1:]))
        self.y = np.concatenate((y_upper[::-1], y_lower[1:]))


class AirfoilGeometry(BladeGeometry):
    """
    Airfoil geometry class.
    
    This class represents an airfoil geometry and provides methods for
    analyzing and manipulating airfoil shapes. It inherits from BladeGeometry
    and adds airfoil-specific functionality.
    
    Attributes:
        All attributes from BladeGeometry plus:
        camber_line: Array of camber line y-coordinates.
        thickness_distribution: Array of thickness values along the chord.
    """
    
    def __init__(self, coords: Optional[np.ndarray] = None, 
                params: Optional[Dict[str, Any]] = None,
                name: str = "unnamed_airfoil"):
        """
        Initialize airfoil geometry.
        
        Args:
            coords: Array of shape (n,2) containing [x,y] coordinates,
                   or None to generate from parameters.
            params: Dictionary of parameters for generating the geometry,
                   or None to use provided coordinates.
            name: Name of the airfoil.
        """
        super().__init__(coords=coords, params=params, name=name)
        
        # Additional airfoil-specific attributes
        self.camber_line = None
        self.thickness_distribution = None
        
        # Calculate camber line and thickness if coordinates are provided
        if self.n_points > 0:
            self._calculate_camber_and_thickness()
    
    def _calculate_camber_and_thickness(self) -> None:
        """
        Calculate camber line and thickness distribution.
        
        This method divides the airfoil into upper and lower surfaces,
        interpolates to common x-locations, and then calculates the
        camber line and thickness distribution.
        """
        if self.n_points < 3:
            logger.warning("Not enough points to calculate camber and thickness")
            return
        
        # Find leading and trailing edge indices
        le_idx = np.argmin(self.x)
        te_idx = np.argmax(self.x)
        
        # Split into upper and lower surfaces
        if self.is_closed:
            # For closed airfoils, need to handle the loop carefully
            if le_idx < te_idx:
                # LE is before TE
                upper_indices = np.arange(le_idx, -1, -1)
                upper_indices = np.append(upper_indices, np.arange(self.n_points-1, te_idx, -1))
                lower_indices = np.arange(le_idx, te_idx + 1)
            else:
                # TE is before LE
                upper_indices = np.arange(le_idx, te_idx - 1, -1)
                lower_indices = np.arange(le_idx, self.n_points - 1)
                lower_indices = np.append(lower_indices, np.arange(0, te_idx + 1))
        else:
            # For open airfoils, simple split at LE
            upper_indices = np.arange(0, le_idx + 1)
            lower_indices = np.arange(le_idx, self.n_points)
        
        # Get upper and lower surface coordinates
        x_upper = self.x[upper_indices]
        y_upper = self.y[upper_indices]
        x_lower = self.x[lower_indices]
        y_lower = self.y[lower_indices]
        
        # Sort by x-coordinate
        upper_sort_idx = np.argsort(x_upper)
        lower_sort_idx = np.argsort(x_lower)
        
        x_upper = x_upper[upper_sort_idx]
        y_upper = y_upper[upper_sort_idx]
        x_lower = x_lower[lower_sort_idx]
        y_lower = y_lower[lower_sort_idx]
        
        # Create common x-coordinates for interpolation
        x_min = max(np.min(x_upper), np.min(x_lower))
        x_max = min(np.max(x_upper), np.max(x_lower))
        
        x_common = np.linspace(x_min, x_max, 100)
        
        # Interpolate upper and lower surfaces to common x-locations
        y_upper_interp = np.interp(x_common, x_upper, y_upper)
        y_lower_interp = np.interp(x_common, x_lower, y_lower)
        
        # Calculate camber line and thickness distribution
        self.camber_line = (y_upper_interp + y_lower_interp) / 2
        self.thickness_distribution = y_upper_interp - y_lower_interp
    
    def get_camber_line(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get the camber line coordinates.
        
        Returns:
            Tuple of (x, y) arrays for camber line.
        """
        if self.camber_line is None:
            self._calculate_camber_and_thickness()
        
        # Create normalized x-coordinates from LE to TE
        x_min, x_max = np.min(self.x), np.max(self.x)
        x_norm = np.linspace(x_min, x_max, len(self.camber_line))
        
        return x_norm, self.camber_line
    
    def get_thickness_distribution(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get the thickness distribution.
        
        Returns:
            Tuple of (x, thickness) arrays.
        """
        if self.thickness_distribution is None:
            self._calculate_camber_and_thickness()
        
        # Create normalized x-coordinates from LE to TE
        x_min, x_max = np.min(self.x), np.max(self.x)
        x_norm = np.linspace(x_min, x_max, len(self.thickness_distribution))
        
        return x_norm, self.thickness_distribution
    
    def get_max_thickness(self) -> Tuple[float, float]:
        """
        Get the maximum thickness and its location.
        
        Returns:
            Tuple of (max_thickness, x_location) where x_location 
            is normalized by chord length.
        """
        if self.thickness_distribution is None:
            self._calculate_camber_and_thickness()
        
        if self.thickness_distribution is not None and len(self.thickness_distribution) > 0:
            max_idx = np.argmax(self.thickness_distribution)
            max_thickness = self.thickness_distribution[max_idx]
            
            # Location normalized by chord
            x_min, x_max = np.min(self.x), np.max(self.x)
            chord = x_max - x_min
            x_norm = np.linspace(0, 1, len(self.thickness_distribution))
            x_location = x_norm[max_idx]
            
            return max_thickness, x_location
        else:
            return 0.0, 0.5
    
    def get_max_camber(self) -> Tuple[float, float]:
        """
        Get the maximum camber and its location.
        
        Returns:
            Tuple of (max_camber, x_location) where x_location 
            is normalized by chord length.
        """
        if self.camber_line is None:
            self._calculate_camber_and_thickness()
        
        if self.camber_line is not None and len(self.camber_line) > 0:
            max_idx = np.argmax(np.abs(self.camber_line))
            max_camber = self.camber_line[max_idx]
            
            # Location normalized by chord
            x_min, x_max = np.min(self.x), np.max(self.x)
            chord = x_max - x_min
            x_norm = np.linspace(0, 1, len(self.camber_line))
            x_location = x_norm[max_idx]
            
            return abs(max_camber), x_location
        else:
            return 0.0, 0.0
    
    @classmethod
    def from_naca4(cls, code: str = '0012', n_points: int = 100, closed_te: bool = True) -> 'AirfoilGeometry':
        """
        Create an airfoil from a NACA 4-digit code.
        
        Args:
            code: NACA 4-digit code (e.g., '0012').
            n_points: Number of points to generate.
            closed_te: Flag indicating if trailing edge should be closed.
            
        Returns:
            AirfoilGeometry object.
        """
        # Create a NACA 4-digit airfoil first
        naca_airfoil = NACA4DigitAirfoil(code, n_points, closed_te)
        
        # Convert to AirfoilGeometry
        return cls(coords=naca_airfoil.get_coordinates(), name=f"NACA{code}")
    
    @classmethod
    def create_naca(cls, code: str = '0012', n_points: int = 100, closed_te: bool = True) -> 'AirfoilGeometry':
        """
        Create an airfoil from a NACA 4-digit code.
        
        This is an alias for from_naca4 for API consistency with BladeGeometry.
        
        Args:
            code: NACA 4-digit code (e.g., '0012').
            n_points: Number of points to generate.
            closed_te: Flag indicating if trailing edge should be closed.
            
        Returns:
            AirfoilGeometry object.
        """
        return cls.from_naca4(code, n_points, closed_te)


class CascadeGeometry:
    """
    Cascade/blade row geometry configuration.
    
    This class represents a cascade of airfoils/blades with periodic
    boundary conditions, typical for turbomachinery applications.
    
    Attributes:
        blade: BladeGeometry object representing the blade profile.
        stagger_angle: Stagger angle in degrees.
        pitch: Blade-to-blade spacing.
        chord: Blade chord length.
        n_blades: Number of blades in the cascade.
    """
    
    def __init__(self, blade: BladeGeometry, stagger_angle: float = 0.0, 
                pitch: float = 1.0, chord: float = 1.0, n_blades: int = 3):
        """
        Initialize cascade geometry.
        
        Args:
            blade: BladeGeometry object representing the blade profile.
            stagger_angle: Stagger angle in degrees (positive clockwise).
            pitch: Blade-to-blade spacing.
            chord: Blade chord length.
            n_blades: Number of blades to include in the cascade.
        """
        self.blade = blade
        self.stagger_angle = stagger_angle
        self.pitch = pitch
        self.chord = chord
        self.n_blades = n_blades
        
        # Normalize and position the blade
        self._normalize_blade()
    
    def _normalize_blade(self) -> None:
        """
        Normalize and position the blade according to cascade parameters.
        
        This method scales the blade to the specified chord length,
        rotates it to the specified stagger angle, and positions it
        appropriately within the cascade.
        """
        if self.blade.n_points == 0:
            logger.warning("Cannot normalize empty blade geometry")
            return
        
        # Make a copy of the original blade
        original_x = self.blade.x.copy()
        original_y = self.blade.y.copy()
        
        # Find current chord length and midpoint
        x_min, x_max = np.min(original_x), np.max(original_x)
        current_chord = x_max - x_min
        mid_x = (x_min + x_max) / 2
        mid_y = np.mean(original_y)
        
        # Scale to desired chord length
        if current_chord > 0:
            scale_factor = self.chord / current_chord
            self.blade.x = mid_x + (original_x - mid_x) * scale_factor
            self.blade.y = mid_y + (original_y - mid_y) * scale_factor
        else:
            logger.warning("Cannot scale blade with zero chord length")
        
        # Rotate to stagger angle
        self.blade.rotate(-self.stagger_angle, (mid_x, mid_y))
        
        # Position at origin
        x_min, x_max = np.min(self.blade.x), np.max(self.blade.x)
        y_min, y_max = np.min(self.blade.y), np.max(self.blade.y)
        
        # Center the blade at the origin
        self.blade.translate(-0.5 * (x_min + x_max), -0.5 * (y_min + y_max))
    
    def get_blade_positions(self) -> List[Tuple[float, float]]:
        """
        Get the positions of all blades in the cascade.
        
        Returns:
            List of (x, y) tuples representing the position of each blade.
        """
        # Calculate the blade-to-blade vector based on pitch and stagger angle
        stagger_rad = np.radians(self.stagger_angle)
        dx = -self.pitch * np.sin(stagger_rad)
        dy = self.pitch * np.cos(stagger_rad)
        
        # Generate blade positions
        positions = []
        for i in range(self.n_blades):
            # Center the cascade around the origin
            offset = (i - (self.n_blades - 1) / 2)
            x = offset * dx
            y = offset * dy
            positions.append((x, y))
        
        return positions
    
    def get_blade_coordinates(self, blade_index: int) -> np.ndarray:
        """
        Get the coordinates of a specific blade in the cascade.
        
        Args:
            blade_index: Index of the blade (0 to n_blades-1).
            
        Returns:
            Array of shape (n,2) containing [x,y] coordinates of the blade.
        """
        if blade_index < 0 or blade_index >= self.n_blades:
            raise ValueError(f"Blade index out of range: {blade_index}")
        
        positions = self.get_blade_positions()
        dx, dy = positions[blade_index]
        
        # Translate the base blade to the specified position
        x = self.blade.x + dx
        y = self.blade.y + dy
        
        return np.column_stack((x, y))
    
    def get_all_blades_coordinates(self) -> List[np.ndarray]:
        """
        Get the coordinates of all blades in the cascade.
        
        Returns:
            List of arrays, each of shape (n,2) containing [x,y] coordinates of each blade.
        """
        return [self.get_blade_coordinates(i) for i in range(self.n_blades)]
    
    def get_tangential_pitch_line_coordinates(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get the coordinates of the tangential pitch lines.
        
        Returns:
            Tuple of (x_coords, y_coords) arrays defining the pitch lines.
        """
        # Calculate the extent of the cascade domain
        all_coords = self.get_all_blades_coordinates()
        min_x = min(np.min(coords[:, 0]) for coords in all_coords)
        max_x = max(np.max(coords[:, 0]) for coords in all_coords)
        min_y = min(np.min(coords[:, 1]) for coords in all_coords)
        max_y = max(np.max(coords[:, 1]) for coords in all_coords)
        
        # Add some margin
        margin = 0.5 * self.chord
        min_x -= margin
        max_x += margin
        min_y -= margin
        max_y += margin
        
        # Create pitch line coordinates
        positions = self.get_blade_positions()
        stagger_rad = np.radians(self.stagger_angle)
        
        # Direction vector parallel to stagger line
        dx = np.cos(stagger_rad)
        dy = np.sin(stagger_rad)
        
        # Normalize direction vector
        norm = np.sqrt(dx**2 + dy**2)
        dx, dy = dx / norm, dy / norm
        
        # Calculate pitch line length
        length = 2.0 * max(max_x - min_x, max_y - min_y)
        
        # Generate pitch line coordinates
        pitch_lines_x = []
        pitch_lines_y = []
        
        for x0, y0 in positions:
            # Create a line centered at the blade position
            t = np.linspace(-length/2, length/2, 100)
            x = x0 + t * dx
            y = y0 + t * dy
            
            pitch_lines_x.append(x)
            pitch_lines_y.append(y)
        
        return pitch_lines_x, pitch_lines_y