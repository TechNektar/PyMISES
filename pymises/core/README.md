# PyMISES Core Module

This directory contains the core functionality of the PyMISES solver. These modules implement the fundamental algorithms and data structures required for airfoil and cascade flow analysis.

## Overview

The core module contains the following primary components:

1. **Grid Generation**: Implementation of structured grid generation for airfoils and cascades.
2. **Euler Solver**: Implementation of the streamline-based Euler equations solver.
3. **Boundary Layer Solver**: Implementation of the integral boundary layer equations.
4. **Coupling Framework**: Implementation of the viscous-inviscid coupling mechanism.
5. **Newton Solver**: Implementation of the Newton iteration framework for nonlinear systems.
6. **Geometry Handling**: Tools for defining and manipulating airfoil and cascade geometries.

## Modules

### `grid.py`

The grid module implements streamline-based structured grid generation for computational domains around airfoils and cascades.

#### Key Classes

- **`StreamlineGrid`**: Represents a structured grid with intrinsic streamlines.
  - Stores grid point coordinates (x, y)
  - Provides methods for grid quality assessment, smoothing, and redistribution
  - Computes grid metrics needed for the flow solver

- **`GridGenerator`**: Creates computational grids for different configurations.
  - Supports O-grid (closed loop around airfoil), C-grid (with wake region), and cascade configurations
  - Implements grid point clustering near leading/trailing edges and walls
  - Provides grid quality improvement methods

#### Grid Generation Workflow

1. Create a geometry (airfoil or cascade)
2. Initialize a `GridGenerator` with the geometry and configuration parameters
3. Call `generate_grid()` with the desired grid type
4. (Optional) Improve grid quality with `improve_grid_quality()`

#### Implementation Details

- Grid points are distributed around airfoils with clustering at leading and trailing edges
- Wall-normal distribution uses hyperbolic tangent or exponential clustering for boundary layer resolution
- For cascades, periodic boundaries are aligned with the stagger angle
- Grid quality metrics include orthogonality, aspect ratio, skewness, and smoothness
- Elliptic grid smoothing is available for grid improvement

### `euler.py`

The Euler solver module implements a streamline-based approach to solving the 2D Euler equations.

#### Key Classes

- **`EulerSolver`**: Solves the Euler equations on a streamline grid.
  - Manages conserved variables (density, momentum, energy)
  - Computes fluxes and residuals
  - Applies boundary conditions
  - Calculates forces and moments

#### Euler Solver Workflow

1. Create a grid
2. Initialize an `EulerSolver` with the grid
3. Add boundary conditions
4. Initialize the flow field with freestream conditions
5. Compute residuals and Jacobian for Newton's method
6. Solve the nonlinear system
7. Compute forces and other output quantities

#### Implementation Details

- Streamline-based discretization reduces the number of variables per node
- Uses the streamtube area-velocity formulation from the original MISES solver
- Artificial dissipation is added in transonic regions using a pressure-switch
- Conservation of mass, momentum (S and N directions), and energy equations
- Supersonic flow handling with shock capturing
- Efficient Jacobian calculation for Newton's method

### `boundary_layer.py`

The boundary layer module implements integral boundary layer equations for viscous effects.

#### Key Classes

- **`BoundaryLayerSolver`**: Solves the integral boundary layer equations.
  - Tracks boundary layer parameters (δ*, θ, H, cf)
  - Handles laminar and turbulent flow regimes
  - Includes transition prediction
  - Detects separation and reattachment

- **`BoundaryLayerFactory`**: Creates specialized boundary layer solvers.
  - Allows creation of solvers for specific configurations
  - Manages transition models and closure coefficients

#### Boundary Layer Solver Workflow

1. Define surface coordinates and edge velocity distribution
2. Create a boundary layer solver with appropriate Reynolds number and configuration
3. Solve the boundary layer equations to obtain δ*, θ, H, and cf
4. (Optional) Extract transition and separation points

#### Implementation Details

- Uses momentum integral and shape parameter equations
- Implements appropriate closure models for laminar and turbulent flow
- Uses modified Abu-Ghannam/Shaw transition model as described by Drela
- Includes lag equation for non-equilibrium effects in turbulent flow
- Special treatment for wake regions
- Robust numerical handling of laminar separation bubbles

### `coupling.py`

The coupling module implements the viscous-inviscid interaction method.

#### Key Classes

- **`CoupledSolver`**: Couples Euler and boundary layer solvers.
  - Manages the interaction between inviscid and viscous solutions
  - Updates displacement thickness and transpiration velocity
  - Integrates boundary layer variables into the global Newton system

#### Coupling Workflow

1. Obtain an initial inviscid solution
2. Initialize the coupled solver with Euler solver and BL factory
3. Set up the global solution vector including BL variables
4. Compute coupled residuals and Jacobian
5. Solve the coupled nonlinear system
6. Extract forces and boundary layer properties

#### Implementation Details

- Uses displacement body concept to represent viscous effects
- Calculates transpiration velocity from displacement thickness
- Fully integrates boundary layer equations into the global Newton system
- Handles wake displacement and circulation effects
- Special treatment for separation regions

### `newton.py`

The Newton solver module implements Newton's method for nonlinear systems.

#### Key Classes

- **`NewtonSolver`**: Solves nonlinear systems using Newton's method.
  - Manages the iterative solution process
  - Handles convergence monitoring
  - Provides relaxation strategies
  - Supports block matrix structure

#### Newton Solver Workflow

1. Define residual and Jacobian functions
2. Create a Newton solver with initial guess
3. Iteratively solve the system
4. Monitor convergence and adjust relaxation if needed

#### Implementation Details

- Supports various linear system solution methods
- Implements adaptive relaxation for improved robustness
- Handles block structure for coupled systems
- Monitors convergence using residual norms
- Implements line search and other globalization strategies

### `geometry.py`

The geometry module implements tools for defining and manipulating airfoil and cascade geometries.

#### Key Classes

- **`BladeGeometry`**: Base class for airfoil/blade geometries.
  - Stores and manipulates coordinates
  - Computes geometric properties
  - Supports parametric representations

- **`AirfoilGeometry`**: Specialized class for single airfoils.
  - Creates standard airfoil shapes (NACA, etc.)
  - Provides methods for airfoil manipulation
  - Extracts leading/trailing edge and surface properties

- **`CascadeGeometry`**: Represents cascade of multiple blades.
  - Manages stagger angle, pitch, and other cascade parameters
  - Transforms between blade and cascade coordinate systems
  - Provides methods for cascade geometry analysis

#### Geometry Workflow

1. Create or load an airfoil geometry
2. (Optional) Modify the geometry for design purposes
3. (Optional) Create a cascade geometry from the airfoil
4. Use the geometry to generate a computational grid

#### Implementation Details

- Supports loading standard airfoil formats (Selig, Lednicer)
- Built-in generation of NACA 4-digit airfoils
- Geometric analysis including curvature, thickness, camber
- Point redistribution for proper grid generation
- Coordinate transformations for cascade configurations
- Support for modal shape functions in design applications

## Dependencies

- NumPy: For array operations
- SciPy: For sparse matrix operations and optimization
- Matplotlib (optional): For visualization

## Example Usage

```python
# Create an airfoil geometry
airfoil = AirfoilGeometry.create_naca('0012', n_points=101)

# Generate a grid
grid_gen = GridGenerator(airfoil, {
    'ni': 101,
    'nj': 41,
    'far_field_distance': 20.0
})
grid = grid_gen.generate_grid(grid_type='o-grid')

# Create Euler solver
euler_solver = EulerSolver(grid)

# Add boundary conditions
from pymises.boundary_conditions.wall import InviscidWallBC
from pymises.boundary_conditions.farfield import VortexFarfieldBC

airfoil_indices = list(range(grid.ni))
wall_bc = InviscidWallBC(airfoil_indices, normal_direction="inner")

farfield_indices = list(range(grid.ni, grid.ni * grid.nj))
farfield_bc = VortexFarfieldBC(
    farfield_indices, 
    mach=0.5, 
    alpha=np.radians(2.0),
    pressure=101325.0,
    temperature=288.15,
    circulation=0.0
)

euler_solver.add_boundary_condition(wall_bc)
euler_solver.add_boundary_condition(farfield_bc)

# Initialize flow field
euler_solver.initialize(
    mach=0.5,
    alpha=np.radians(2.0),
    p0=101325.0,
    T0=288.15
)

# Create Newton solver and solve
newton_solver = NewtonSolver(
    residual_function=euler_solver.compute_residuals,
    jacobian_function=euler_solver.compute_jacobian,
    solution=euler_solver.get_solution_vector()
)

final_solution, convergence = newton_solver.solve(
    max_iter=50,
    tolerance=1e-6,
    relaxation=0.7
)

# Update solution and compute forces
euler_solver.set_solution_vector(final_solution)
forces = euler_solver.compute_forces()
print(f"Lift coefficient: {forces['cl']}")
print(f"Drag coefficient: {forces['cd']}")
```

## Mathematical Background

The PyMISES implementation follows the mathematical formulation of the original MISES code:

1. **Streamline-based Euler Discretization**:
   - Conservation of mass through streamtube area variation
   - S-momentum equation along streamlines
   - N-momentum equation normal to streamlines
   - Energy equation with constant stagnation enthalpy

2. **Integral Boundary Layer Equations**:
   - Momentum integral equation
   - Shape parameter equation
   - Lag equation for Reynolds stresses
   - Appropriate closure relations for laminar and turbulent flow

3. **Viscous-Inviscid Coupling**:
   - Displacement body concept
   - Transpiration velocity formulation
   - Fully-coupled Newton solution approach

4. **Modified Transition Model**:
   - H-based Abu-Ghannam/Shaw transition criterion
   - Combined amplification rate formulation
   - Edge velocity correction for pressure gradients

For a detailed discussion of the mathematical foundation, refer to:
- Drela, M., "A Two-Dimensional Viscous Aerodynamic Design and Analysis Code", AIAA Paper 86-0424, 1986.
- Drela, M., "An Integral Boundary Layer Formulation for Blunt Trailing Edges", AIAA Paper 89-2166, 1989.
- Drela, M., "MISES Implementation of Modified Abu-Ghannam/Shaw Transition Criterion", MIT ACDL Report, 1998.
