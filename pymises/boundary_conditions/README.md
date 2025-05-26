# PyMISES Boundary Conditions Module

This directory contains implementations of various boundary conditions used in the PyMISES solver. These boundary conditions are essential for properly defining the flow domain and enabling both direct analysis and inverse design capabilities.

## Overview

The boundary conditions module provides:

1. **Wall Boundary Conditions**: For solid surfaces like airfoils and cascade blades
2. **Farfield Boundary Conditions**: For domain outer boundaries
3. **Periodicity Boundary Conditions**: For cascade flow with repeating patterns
4. **Inverse Design Boundary Conditions**: For shape design based on pressure specifications

## Modules

### `wall.py`

This module implements boundary conditions for solid wall surfaces such as airfoils and cascade blades.

#### Key Classes

- **`WallBoundaryCondition`**: Abstract base class for all wall boundary conditions.
  - Defines the common interface for wall boundary conditions
  - Implements validation and utility methods
  - Manages grid indices along the wall

- **`InviscidWallBC`**: Implements inviscid (slip) wall boundary condition.
  - Enforces zero normal velocity at the wall
  - Preserves tangential velocity component
  - Computes residual and Jacobian contributions

- **`ViscousWallBC`**: Implements viscous wall boundary condition with displacement effects.
  - Accounts for boundary layer displacement thickness
  - Imposes transpiration velocity to represent viscous effects
  - Supports varying displacement thickness along the wall

#### Implementation Details

- For inviscid walls, the flow tangency condition is enforced by setting the normal velocity component to zero
- For viscous walls, the displacement thickness is translated into a transpiration velocity
- The transpiration velocity is calculated based on the derivative of the product of edge velocity and displacement thickness
- Jacobian contributions are computed analytically for use in the Newton solver
- Both conditions support inner and outer normal directions for different grid topologies

### `farfield.py`

This module implements boundary conditions for the far-field (external) boundaries of the domain.

#### Key Classes

- **`FarfieldBoundaryCondition`**: Abstract base class for all far-field boundary conditions.
  - Defines the common interface for far-field conditions
  - Manages grid indices along the far-field boundary
  - Provides methods for flow variable calculation

- **`SubsonicInflow`**: Implements boundary condition for subsonic inflow regions.
  - Specifies total pressure, total temperature, and flow direction
  - Extrapolates one characteristic from the interior

- **`SubsonicOutflow`**: Implements boundary condition for subsonic outflow regions.
  - Specifies static pressure (back pressure)
  - Extrapolates other variables from the interior

- **`VortexFarfieldBC`**: Implements far-field boundary condition with vortex correction.
  - Accounts for circulation around airfoil using a potential vortex
  - Provides accurate representation for external aerodynamics
  - Supports variable circulation for lift convergence

- **`CharacteristicBC`**: Implements characteristic-based boundary condition.
  - Uses characteristic theory to determine which variables to specify
  - Automatically handles inflow and outflow regions
  - Provides robust treatment for complex flow domains

#### Implementation Details

- Far-field conditions use characteristic theory to determine which variables to specify and which to extrapolate
- For subsonic inflows, total pressure, total temperature, and flow direction are specified
- For subsonic outflows, static pressure is specified while other variables are extrapolated
- Vortex correction accounts for the circulation around an airfoil using potential flow theory
- Characteristic-based approach automatically adapts to local flow direction and Mach number
- Implementation includes both strong (direct) and weak (via residuals) enforcement methods

### `periodicity.py`

This module implements boundary conditions for enforcing periodicity in cascade configurations.

#### Key Classes

- **`PeriodicityBC`**: Implements standard periodicity for steady flow.
  - Enforces equal flow conditions at corresponding points on periodic boundaries
  - Accounts for geometric transformation between periodic boundaries
  - Provides residual and Jacobian contributions for the Newton solver

- **`PhaseLagPeriodicityBC`**: Implements phase-lagged periodicity for unsteady flow.
  - Accounts for temporal phase shift between blade passages
  - Supports time-accurate and time-spectral formulations
  - Models blade row interactions with different blade counts

#### Implementation Details

- Periodicity conditions match flow variables between corresponding points on periodic boundaries
- For cascades, the matching accounts for the pitch vector (translation between adjacent blades)
- Transformation of velocity components accounts for the relative orientation of the periodic boundaries
- For phase-lagged conditions, temporal storage and interpolation are used to handle the phase shift
- Both conditions provide analytical Jacobian contributions for inclusion in the global system
- Implementation supports both direct replacement and residual formulations

### `inverse.py`

This module implements boundary conditions for inverse design capabilities.

#### Key Classes

- **`PressureSpecificationBC`**: Implements pressure-based inverse design.
  - Allows specification of target pressure distribution
  - Modifies geometry to achieve the target pressure
  - Preserves required geometrical constraints

- **`GeometricConstraintBC`**: Implements geometric constraints for inverse design.
  - Enforces fixed geometry at specific points (e.g., leading and trailing edges)
  - Maintains geometric continuity during the design process
  - Works in conjunction with pressure specification

- **`MixedInverseBC`**: Implements mixed inverse approach.
  - Combines pressure specification and fixed geometry regions
  - Allows for practical design with partial pressure specification
  - Maintains smooth transitions between design and analysis regions

- **`ModalInverseDesignBC`**: Implements modal shape function approach.
  - Uses shape functions (like Hicks-Henne bumps) for geometry modification
  - Reduces the design space to reasonable dimensionality
  - Ensures smooth geometries with physically realizable shapes

#### Implementation Details

- Pressure specification replaces the flow tangency condition with a pressure constraint
- Geometry modifications are solved as part of the global Newton system
- Geometric constraints fix specific points (like leading and trailing edges) to maintain closure
- Mixed inverse allows pressure specification on parts of the airfoil while keeping other parts fixed
- Modal approach uses shape functions to restrict geometry changes to smooth, well-behaved forms
- Jacobian contributions include sensitivity of pressure to geometry changes
- Implementation handles the additional degrees of freedom introduced by geometry movement

## Dependencies

- NumPy: For array operations
- SciPy: For sparse matrix operations

## Boundary Condition Interface

All boundary conditions implement a common interface:

```python
class BoundaryCondition:
    """Base class for all boundary conditions."""
    
    def __init__(self, grid_indices, **kwargs):
        """Initialize boundary condition with grid indices and parameters."""
        self.grid_indices = grid_indices
        # Additional parameters specific to boundary condition type
    
    def apply(self, solution):
        """Apply boundary condition to the solution."""
        # Modify solution in-place or return modified solution
    
    def get_residual_contributions(self, solution):
        """Get residual contributions from this boundary condition."""
        # Return residual terms for the Newton solver
    
    def get_jacobian_contributions(self, solution):
        """Get Jacobian contributions from this boundary condition."""
        # Return Jacobian terms for the Newton solver
```

## Example Usage

### Inviscid Wall Boundary Condition

```python
from pymises.boundary_conditions.wall import InviscidWallBC

# Grid indices for the airfoil surface
airfoil_indices = list(range(grid.ni))

# Create inviscid wall boundary condition
wall_bc = InviscidWallBC(airfoil_indices, normal_direction="inner")

# Add boundary condition to Euler solver
euler_solver.add_boundary_condition(wall_bc)
```

### Far-field Boundary Condition

```python
from pymises.boundary_conditions.farfield import VortexFarfieldBC

# Grid indices for the far-field boundary
farfield_indices = list(range(grid.ni, grid.ni * grid.nj))

# Create far-field boundary condition
farfield_bc = VortexFarfieldBC(
    farfield_indices, 
    mach=0.5, 
    alpha=np.radians(2.0),
    pressure=101325.0,
    temperature=288.15,
    circulation=0.0
)

# Add boundary condition to Euler solver
euler_solver.add_boundary_condition(farfield_bc)
```

### Periodicity for Cascade

```python
from pymises.boundary_conditions.periodicity import PeriodicityBC

# Grid indices for upper and lower periodic boundaries
upper_indices = list(range(grid.ni + 2*grid.nj, grid.ni + 3*grid.nj))
lower_indices = list(range(grid.ni + 3*grid.nj, grid.ni + 4*grid.nj))

# Pitch vector based on stagger angle
stagger_angle = np.radians(30.0)
pitch_vector = np.array([np.sin(stagger_angle), np.cos(stagger_angle)])

# Create periodicity boundary condition
periodicity_bc = PeriodicityBC(
    upper_indices,
    lower_indices,
    pitch_vector=pitch_vector
)

# Add boundary condition to Euler solver
euler_solver.add_boundary_condition(periodicity_bc)
```

### Inverse Design with Pressure Specification

```python
from pymises.boundary_conditions.inverse import PressureSpecificationBC, GeometricConstraintBC

# Grid indices for the airfoil surface
airfoil_indices = list(range(grid.ni))

# Create pressure specification boundary condition
pressure_bc = PressureSpecificationBC(
    grid_indices=airfoil_indices,
    target_pressure=target_pressure,
    normal_direction="inner"
)

# Add geometric constraints for leading and trailing edges
constraint_indices = [0, grid.ni//2, grid.ni-1]
constraint_bc = GeometricConstraintBC(
    grid_indices=[airfoil_indices[i] for i in constraint_indices],
    normal_direction="inner"
)

# Add boundary conditions to Euler solver
euler_solver.add_boundary_condition(pressure_bc)
euler_solver.add_boundary_condition(constraint_bc)
```

## Mathematical Background

### Wall Boundary Conditions

- **Inviscid Wall**: Flow tangency condition
  
  $\mathbf{v} \cdot \mathbf{n} = 0$
  
  where $\mathbf{v}$ is the velocity vector and $\mathbf{n}$ is the wall normal vector.

- **Viscous Wall**: Transpiration velocity to account for boundary layer displacement
  
  $\mathbf{v} \cdot \mathbf{n} = v_\text{transpiration}$
  
  where $v_\text{transpiration} = \frac{d}{ds}(u_e \delta^*)$ with $u_e$ being the edge velocity and $\delta^*$ the displacement thickness.

### Far-field Boundary Conditions

- **Subsonic Inflow**: Specify stagnation properties and flow direction
  
  $p_0$, $T_0$, and flow angle are specified, with one characteristic extrapolated from the interior.

- **Subsonic Outflow**: Specify back pressure
  
  $p$ is specified, with other variables extrapolated from the interior.

- **Vortex Correction**: Potential vortex influence
  
  $v_\theta = \frac{\Gamma}{2\pi r}$
  
  where $\Gamma$ is the circulation and $r$ is the distance from the airfoil quarter-chord.

### Periodicity Boundary Conditions

- **Standard Periodicity**: Matching between corresponding points
  
  $\mathbf{q}(\mathbf{x}) = \mathbf{q}(\mathbf{x} + \mathbf{\Delta})$
  
  where $\mathbf{q}$ is the vector of flow variables and $\mathbf{\Delta}$ is the pitch vector.

- **Phase-Lagged Periodicity**: Accounting for temporal phase shift
  
  $\mathbf{q}(\mathbf{x}, t) = \mathbf{q}(\mathbf{x} + \mathbf{\Delta}, t + \Delta t)$
  
  where $\Delta t$ is the temporal phase shift.

### Inverse Design Boundary Conditions

- **Pressure Specification**: Enforce target pressure instead of flow tangency
  
  $p(\mathbf{x}) = p_\text{target}(\mathbf{x})$
  
  with geometry degrees of freedom to satisfy the pressure constraint.

- **Modal Approach**: Restrict geometry modification using shape functions
  
  $\mathbf{x}_\text{new} = \mathbf{x}_\text{original} + \sum_{i=1}^N \alpha_i \mathbf{f}_i(\mathbf{x})$
  
  where $\mathbf{f}_i$ are shape functions (like Hicks-Henne bumps) and $\alpha_i$ are design variables.

For more detailed descriptions, refer to:
- Drela, M., "A Two-Dimensional Viscous Aerodynamic Design and Analysis Code", AIAA Paper 86-0424, 1986.
- Giles, M.B., "Non-reflecting Boundary Conditions for Euler Equation Calculations", AIAA Journal, Vol. 28, No. 12, 1990.
- Drela, M., "Design and Optimization Method for Multi-Element Airfoils", AIAA Paper 93-0969, 1993.
