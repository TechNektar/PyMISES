# PyMISES Numerics Module

This directory contains implementations of numerical methods used in the PyMISES solver. These modules provide the essential numerical tools for solving the Euler and boundary layer equations.

## Overview

The numerics module provides:

1. **Jacobian calculation**: Methods for computing Jacobian matrices
2. **Linear solvers**: Algorithms for solving linear systems, particularly with block structure
3. **Stabilization techniques**: Methods for enhancing numerical stability

## Modules

### `jacobian.py`

This module provides tools for calculating and assembling Jacobian matrices for nonlinear systems.

#### Key Classes

- **`JacobianCalculator`**: Computes Jacobian matrices using various methods.
  - Supports finite difference approximation
  - Implements complex step differentiation for higher accuracy
  - Allows for analytical Jacobian specification
  - Provides options for sparsity pattern exploitation

- **`JacobianAssembler`**: Builds global Jacobian matrices from block components.
  - Manages block structure for efficient assembly
  - Supports sparse storage formats
  - Provides methods for adding component matrices
  - Creates the full Jacobian matrix needed for Newton's method

#### Implementation Details

- Multiple Jacobian calculation methods are supported:
  1. **Finite difference**: First-order approximation using forward differences
  2. **Complex step**: Higher accuracy method using complex perturbation
  3. **Analytical**: Direct calculation using user-provided function
  
- The finite difference method uses:
  
  $\frac{\partial F_i}{\partial x_j} \approx \frac{F_i(x+h e_j) - F_i(x)}{h}$
  
  where $h$ is a small step size and $e_j$ is the unit vector in the $j$-th direction.

- The complex step method uses:
  
  $\frac{\partial F_i}{\partial x_j} \approx \frac{\text{Im}[F_i(x+ih e_j)]}{h}$
  
  which avoids subtraction errors and provides higher accuracy.

- For large systems, sparsity patterns are exploited to reduce computation:
  1. Identify non-zero elements in the Jacobian
  2. Only compute derivatives for those elements
  3. Assemble the sparse matrix in CSR format

- The JacobianAssembler uses block structure to build the global Jacobian:
  1. Each block represents derivatives between specific variable groups
  2. Blocks are added to the global matrix at appropriate positions
  3. The resulting matrix preserves the block structure for efficient solution

### `linear_solver.py`

This module implements various algorithms for solving linear systems, particularly those arising in Newton's method.

#### Key Classes

- **`LinearSolver`**: Base class for linear solvers.
  - Defines common interface for all solvers
  - Provides functionality for preconditioner creation
  - Includes error checking and performance monitoring

- **`DirectSolver`**: Implements direct solution methods.
  - LU decomposition for dense matrices
  - Sparse direct solvers for large systems
  - Specialized methods for banded matrices

- **`IterativeSolver`**: Implements iterative solution methods.
  - GMRES for non-symmetric systems
  - Conjugate Gradient for symmetric systems
  - Includes convergence monitoring and control

- **`BlockEliminationSolver`**: Specialized solver for block-structured matrices.
  - Takes advantage of MISES-specific block structure
  - More efficient than general direct or iterative methods
  - Preserves block sparsity pattern during elimination

#### Implementation Details

- Direct solvers use LU decomposition for general matrices:
  
  $A = LU$
  
  followed by forward and backward substitution:
  
  $Ly = b$ and $Ux = y$

- For sparse matrices, specialized sparse direct solvers are used:
  1. Reordering to minimize fill-in
  2. Symbolic factorization to determine non-zero pattern
  3. Numerical factorization
  4. Solution via substitution

- Iterative solvers use Krylov subspace methods:
  - GMRES for non-symmetric systems
  - Preconditioned with incomplete LU factorization
  - Restart capability for memory efficiency
  - Convergence monitoring based on residual reduction

- The block elimination solver exploits the structure of the MISES system:
  1. The matrix is partitioned into blocks based on variable types
  2. Block-wise forward elimination reduces the system
  3. Global variables are solved first
  4. Local variables are found via back-substitution
  5. This approach is much more efficient for the specific structure

- Preconditioners are essential for iterative methods:
  1. Incomplete LU (ILU) for general systems
  2. Block Jacobi for systems with block structure
  3. Multigrid for systems with multiple scales

### `stabilization.py`

This module provides methods for enhancing the numerical stability of the flow solution.

#### Key Classes

- **`ArtificialDissipation`**: Implements artificial dissipation models.
  - Adds dissipation terms to stabilize shock capturing
  - Pressure-switch controls dissipation application
  - Parameters control dissipation strength and localization

- **`Relaxation`**: Implements solution relaxation techniques.
  - Under-relaxation for enhanced stability
  - Line relaxation for strongly coupled directions
  - Adaptive relaxation based on solution behavior

- **`Limiter`**: Implements flux and gradient limiters.
  - Prevents oscillations near discontinuities
  - Preserves monotonicity in the solution
  - Various formulations for different applications

#### Implementation Details

- Artificial dissipation for transonic flow follows Jameson's approach:
  
  $D_i = \epsilon^{(2)} (u_{i+1} - 2u_i + u_{i-1}) - \epsilon^{(4)} (u_{i+2} - 4u_{i+1} + 6u_i - 4u_{i-1} + u_{i-2})$
  
  where $\epsilon^{(2)}$ is controlled by a pressure switch:
  
  $\epsilon^{(2)}_i = \kappa^{(2)} \max(\nu_{i-1}, \nu_i, \nu_{i+1})$
  
  with $\nu_i = \left|\frac{p_{i+1} - 2p_i + p_{i-1}}{p_{i+1} + 2p_i + p_{i-1}}\right|$

- For MISES-specific implementation, the artificial dissipation takes the form:
  
  $\mu = \begin{cases}
  0, & M < M_c \\
  C_d(1-(M_c/M)^2)/(2M^2), & M > M_c
  \end{cases}$
  
  where $M_c$ is a threshold Mach number (typically 0.9-0.95) and $C_d$ is a dissipation coefficient.

- Under-relaxation for Newton iterations:
  
  $x^{n+1} = x^n + \omega \Delta x$
  
  where $\omega \in (0,1]$ is the relaxation factor.

- Adaptive relaxation strategies adjust $\omega$ based on solution behavior:
  1. Start with conservative (small) $\omega$
  2. Increase $\omega$ as the solution converges
  3. Decrease $\omega$ if divergence is detected
  4. Use line search to optimize $\omega$ for each step

- Limiters for gradient reconstruction:
  - Minmod: $\phi(r) = \max(0, \min(1, r))$
  - Van Leer: $\phi(r) = \frac{r + |r|}{1 + |r|}$
  - Superbee: $\phi(r) = \max(0, \min(2r, 1), \min(r, 2))$

## Dependencies

- NumPy: For array operations
- SciPy: For sparse matrix operations and linear algebra

## Example Usage

### Calculating Jacobian Matrix

```python
from pymises.numerics.jacobian import JacobianCalculator

# Define residual function
def residual_function(x):
    return np.array([
        x[0]**2 + x[1]**2 - 4,  # Circle equation
        x[0] - x[1]             # Line equation
    ])

# Create Jacobian calculator using finite differences
calculator = JacobianCalculator(method='finite_difference', step_size=1e-6)

# Calculate Jacobian at point [1.0, 1.0]
x0 = np.array([1.0, 1.0])
jacobian = calculator.calculate(residual_function, x0)

print(jacobian.toarray())
# Output:
# [[2.0, 2.0],
#  [1.0, -1.0]]
```

### Solving Linear System with Block Elimination

```python
from pymises.numerics.linear_solver import BlockEliminationSolver
import numpy as np
import scipy.sparse as sp

# Create a block-structured matrix
# [ A B ]
# [ C D ]
A = np.array([[2.0, 1.0], [1.0, 3.0]])
B = np.array([[0.5, 0.0], [0.0, 0.5]])
C = np.array([[0.2, 0.1], [0.1, 0.2]])
D = np.array([[3.0, 1.0], [1.0, 2.0]])

# Assemble full matrix
matrix = np.block([
    [A, B],
    [C, D]
])

# Convert to sparse format
A_sparse = sp.csr_matrix(A)
B_sparse = sp.csr_matrix(B)
C_sparse = sp.csr_matrix(C)
D_sparse = sp.csr_matrix(D)

# Create right-hand side
b = np.array([1.0, 2.0, 3.0, 4.0])

# Create block elimination solver
solver = BlockEliminationSolver(
    block_structure=[[0, 1], [2, 3]],  # Variable indices in each block
    global_indices=[0, 1]              # Global variables
)

# Provide block matrices
solver.set_blocks(A_sparse, B_sparse, C_sparse, D_sparse)

# Solve the system
x = solver.solve(b)

# Verify solution
print(np.allclose(matrix @ x, b))
# Output: True
```

### Using Artificial Dissipation

```python
from pymises.numerics.stabilization import ArtificialDissipation

# Create artificial dissipation model
dissipation = ArtificialDissipation(
    threshold_mach=0.95,
    dissipation_coeff=0.5
)

# Calculate artificial dissipation flux
# Grid and solution data
mach = np.array([0.8, 0.9, 1.1, 1.2, 1.0])
pressure = np.array([100.0, 95.0, 80.0, 75.0, 85.0])
density = np.array([1.2, 1.15, 1.0, 0.95, 1.05])

# Calculate dissipation term
diss_flux = dissipation.calculate(mach, pressure, density)

# Apply to residual calculation
residual = np.zeros_like(mach)
# ... (calculate physical fluxes)
residual += diss_flux
```

## Mathematical Background

### Jacobian Calculation

Finite difference approximation:

$\frac{\partial F_i}{\partial x_j} \approx \frac{F_i(x+h e_j) - F_i(x)}{h}$

Complex step method:

$\frac{\partial F_i}{\partial x_j} \approx \frac{\text{Im}[F_i(x+ih e_j)]}{h}$

Choice of step size $h$ involves a trade-off between truncation error and roundoff error.

### Linear System Solution

Direct methods solve $Ax = b$ through factorization:

$A = LU$ (LU decomposition)
$Ly = b$ (Forward substitution)
$Ux = y$ (Backward substitution)

Iterative methods use approximation sequence:

$x^{k+1} = x^k + M^{-1}(b - Ax^k)$

where $M$ is the preconditioner.

Block elimination for systems with structure:

$\begin{bmatrix} A & B \\ C & D \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} b_1 \\ b_2 \end{bmatrix}$

is solved by:

$S = D - CA^{-1}B$ (Schur complement)
$Sx_2 = b_2 - CA^{-1}b_1$ (Solve for $x_2$)
$Ax_1 = b_1 - Bx_2$ (Solve for $x_1$)

### Artificial Dissipation

The MISES artificial dissipation model adds a bulk viscosity term in regions of supersonic flow:

$\mu = \begin{cases}
0, & M < M_c \\
C_d(1-(M_c/M)^2)/(2M^2), & M > M_c
\end{cases}$

This is incorporated into the streamtube momentum equations to stabilize shock capturing.

For more detailed descriptions, refer to:
- Drela, M., "A Two-Dimensional Viscous Aerodynamic Design and Analysis Code", AIAA Paper 86-0424, 1986.
- Jameson, A., Schmidt, W., and Turkel, E., "Numerical Solutions of the Euler Equations by Finite Volume Methods Using Runge-Kutta Time-Stepping Schemes", AIAA Paper 81-1259, 1981.
- Saad, Y., "Iterative Methods for Sparse Linear Systems", SIAM, 2003.
