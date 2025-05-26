"""
Linear solver module for PyMISES.

This module provides classes and functions for solving linear systems
arising in the Newton method for the MISES solver, including specialized
block elimination techniques for the specific structure of the problem.
"""

import numpy as np
from typing import Dict, Tuple, List, Union, Optional, Any, Callable
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import time

from pymises.numerics.jacobian import BlockMatrixOperations, JacobianPreconditioner
from pymises.utils.logger import get_logger

logger = get_logger(__name__)

class LinearSolver:
    """
    Linear solver for systems arising in the Newton method.
    
    This class provides methods for solving the linear systems that
    arise in each iteration of Newton's method, including specialized
    block elimination techniques for the MISES solver.
    
    Attributes:
        config: Configuration dictionary with solver parameters.
        block_ops: Block matrix operations helper.
        preconditioner: Preconditioner for iterative solvers.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize linear solver.
        
        Args:
            config: Configuration dictionary with solver parameters:
                - solver_type: Type of linear solver ('direct', 'iterative', 'block_elimination')
                - iterative_method: Method for iterative solvers ('gmres', 'bicgstab')
                - tolerance: Convergence tolerance for iterative solvers
                - max_iterations: Maximum iterations for iterative solvers
                - use_preconditioner: Whether to use preconditioning with iterative solvers
                - matrix_format: Format for matrix storage ('dense', 'sparse', 'block')
        """
        self.config = config or {}
        
        # Solver parameters
        self.solver_type = self.config.get('solver_type', 'block_elimination')
        self.iterative_method = self.config.get('iterative_method', 'gmres')
        self.tolerance = self.config.get('tolerance', 1e-8)
        self.max_iterations = self.config.get('max_iterations', 100)
        self.use_preconditioner = self.config.get('use_preconditioner', True)
        self.matrix_format = self.config.get('matrix_format', 'sparse')
        
        # Initialize helpers
        block_config = {
            'matrix_format': self.matrix_format
        }
        self.block_ops = BlockMatrixOperations(block_config)
        
        precond_config = self.config.get('preconditioner_config', {})
        self.preconditioner = JacobianPreconditioner(precond_config)
        
        logger.info(f"Initialized linear solver: type={self.solver_type}, format={self.matrix_format}")
    
    def solve(self, matrix: Any, rhs: np.ndarray) -> Tuple[np.ndarray, Dict[str, Any]]:
        """
        Solve linear system Ax = b.
        
        Args:
            matrix: System matrix (can be dense, sparse, or block structure).
            rhs: Right-hand side vector.
            
        Returns:
            Tuple of (solution, info_dict).
        """
        # Start timer
        start_time = time.time()
        
        # Dispatch to appropriate solver
        if self.solver_type == 'direct':
            solution, info = self._solve_direct(matrix, rhs)
        elif self.solver_type == 'iterative':
            solution, info = self._solve_iterative(matrix, rhs)
        elif self.solver_type == 'block_elimination':
            solution, info = self._solve_block_elimination(matrix, rhs)
        else:
            logger.warning(f"Unknown solver type: {self.solver_type}, using direct method")
            solution, info = self._solve_direct(matrix, rhs)
        
        # Calculate solution time
        solution_time = time.time() - start_time
        
        # Add timing information
        info['solution_time'] = solution_time
        
        return solution, info
    
    def _solve_direct(self, matrix: Any, rhs: np.ndarray) -> Tuple[np.ndarray, Dict[str, Any]]:
        """
        Solve linear system using direct methods.
        
        Args:
            matrix: System matrix.
            rhs: Right-hand side vector.
            
        Returns:
            Tuple of (solution, info_dict).
        """
        # Different handling based on matrix type
        if isinstance(matrix, np.ndarray):
            # Dense matrix
            try:
                solution = np.linalg.solve(matrix, rhs)
                info = {'success': True, 'method': 'numpy.linalg.solve'}
            except np.linalg.LinAlgError as e:
                logger.warning(f"Failed to solve dense system: {str(e)}")
                # Try pseudo-inverse as fallback
                solution = np.linalg.lstsq(matrix, rhs, rcond=None)[0]
                info = {'success': True, 'method': 'numpy.linalg.lstsq', 'warning': str(e)}
        
        elif sp.issparse(matrix):
            # Sparse matrix
            try:
                solution = spla.spsolve(matrix, rhs)
                info = {'success': True, 'method': 'scipy.sparse.linalg.spsolve'}
            except Exception as e:
                logger.warning(f"Failed to solve sparse system: {str(e)}")
                # Try GMRES as fallback
                solution, exitCode = spla.gmres(matrix, rhs, tol=self.tolerance, atol=0)
                if exitCode == 0:
                    info = {'success': True, 'method': 'scipy.sparse.linalg.gmres (fallback)'}
                else:
                    # Last resort: convert to dense if small enough
                    if matrix.shape[0] <= 1000:
                        try:
                            dense_matrix = matrix.toarray()
                            solution = np.linalg.solve(dense_matrix, rhs)
                            info = {'success': True, 'method': 'dense conversion (fallback)'}
                        except Exception as e2:
                            logger.error(f"Failed all solution attempts: {str(e2)}")
                            solution = np.zeros_like(rhs)
                            info = {'success': False, 'error': str(e2)}
                    else:
                        logger.error(f"Matrix too large for dense conversion")
                        solution = np.zeros_like(rhs)
                        info = {'success': False, 'error': str(e)}
        
        elif isinstance(matrix, dict) and matrix.get('type') == 'block':
            # Block matrix structure
            # Convert to appropriate format and solve
            if self.matrix_format == 'sparse':
                # Convert blocks to sparse matrix
                blocks = matrix['blocks']
                block_pattern = matrix['block_pattern']
                n = matrix['n_total']
                
                sparse_matrix = sp.lil_matrix((n, n))
                
                # Fill in blocks
                for block_name, (i_range, j_range) in block_pattern.items():
                    if block_name in blocks:
                        i_start, i_end = i_range
                        j_start, j_end = j_range
                        block = blocks[block_name]
                        
                        # Convert block to appropriate format
                        if sp.issparse(block):
                            block_array = block.toarray()
                        else:
                            block_array = block
                        
                        # Fill in the sparse matrix
                        sparse_matrix[i_start:i_end, j_start:j_end] = block_array
                
                # Convert to CSR for efficient operations
                sparse_matrix = sparse_matrix.tocsr()
                
                # Solve the sparse system
                solution, info = self._solve_direct(sparse_matrix, rhs)
                info['method'] = 'block_to_sparse conversion'
            
            else:
                # Use block elimination method
                solution, info = self._solve_block_elimination(matrix, rhs)
        
        else:
            # Unsupported matrix type
            logger.error(f"Unsupported matrix type: {type(matrix)}")
            solution = np.zeros_like(rhs)
            info = {'success': False, 'error': f"Unsupported matrix type: {type(matrix)}"}
        
        return solution, info
    
    def _solve_iterative(self, matrix: Any, rhs: np.ndarray) -> Tuple[np.ndarray, Dict[str, Any]]:
        """
        Solve linear system using iterative methods.
        
        Args:
            matrix: System matrix.
            rhs: Right-hand side vector.
            
        Returns:
            Tuple of (solution, info_dict).
        """
        # Ensure matrix is in sparse format
        if not sp.issparse(matrix):
            if isinstance(matrix, np.ndarray):
                matrix = sp.csr_matrix(matrix)
            elif isinstance(matrix, dict) and matrix.get('type') == 'block':
                # Convert block structure to sparse matrix
                blocks = matrix['blocks']
                block_pattern = matrix['block_pattern']
                n = matrix['n_total']
                
                sparse_matrix = sp.lil_matrix((n, n))
                
                # Fill in blocks
                for block_name, (i_range, j_range) in block_pattern.items():
                    if block_name in blocks:
                        i_start, i_end = i_range
                        j_start, j_end = j_range
                        block = blocks[block_name]
                        
                        # Convert block to appropriate format
                        if sp.issparse(block):
                            block_array = block.toarray()
                        else:
                            block_array = block
                        
                        # Fill in the sparse matrix
                        sparse_matrix[i_start:i_end, j_start:j_end] = block_array
                
                # Convert to CSR for efficient operations
                matrix = sparse_matrix.tocsr()
            else:
                logger.error(f"Cannot convert to sparse format: {type(matrix)}")
                return np.zeros_like(rhs), {'success': False, 'error': "Cannot convert to sparse format"}
        
        # Create preconditioner if needed
        M = None
        if self.use_preconditioner:
            M = self.preconditioner.create_preconditioner(matrix)
        
        # Choose iterative method
        if self.iterative_method == 'gmres':
            # GMRES method
            solution, exitCode = spla.gmres(
                matrix, rhs,
                tol=self.tolerance,
                atol=0,
                maxiter=self.max_iterations,
                M=M
            )
            
            if exitCode == 0:
                info = {
                    'success': True,
                    'method': 'gmres',
                    'preconditioned': self.use_preconditioner
                }
            else:
                logger.warning(f"GMRES failed to converge in {self.max_iterations} iterations")
                info = {
                    'success': False,
                    'method': 'gmres',
                    'preconditioned': self.use_preconditioner,
                    'error': f"Failed to converge in {self.max_iterations} iterations"
                }
        
        elif self.iterative_method == 'bicgstab':
            # BiCGSTAB method
            solution, exitCode = spla.bicgstab(
                matrix, rhs,
                tol=self.tolerance,
                atol=0,
                maxiter=self.max_iterations,
                M=M
            )
            
            if exitCode == 0:
                info = {
                    'success': True,
                    'method': 'bicgstab',
                    'preconditioned': self.use_preconditioner
                }
            else:
                logger.warning(f"BiCGSTAB failed to converge in {self.max_iterations} iterations")
                info = {
                    'success': False,
                    'method': 'bicgstab',
                    'preconditioned': self.use_preconditioner,
                    'error': f"Failed to converge in {self.max_iterations} iterations"
                }
        
        else:
            # Unknown method, fallback to GMRES
            logger.warning(f"Unknown iterative method: {self.iterative_method}, using GMRES")
            solution, exitCode = spla.gmres(
                matrix, rhs,
                tol=self.tolerance,
                atol=0,
                maxiter=self.max_iterations,
                M=M
            )
            
            if exitCode == 0:
                info = {
                    'success': True,
                    'method': 'gmres (fallback)',
                    'preconditioned': self.use_preconditioner
                }
            else:
                logger.warning(f"GMRES failed to converge in {self.max_iterations} iterations")
                info = {
                    'success': False,
                    'method': 'gmres (fallback)',
                    'preconditioned': self.use_preconditioner,
                    'error': f"Failed to converge in {self.max_iterations} iterations"
                }
        
        return solution, info
    
    def _solve_block_elimination(self, matrix: Any, rhs: np.ndarray) -> Tuple[np.ndarray, Dict[str, Any]]:
        """
        Solve linear system using block elimination.
        
        Args:
            matrix: System matrix (block structure).
            rhs: Right-hand side vector.
            
        Returns:
            Tuple of (solution, info_dict).
        """
        # Check if matrix is in block format
        if not (isinstance(matrix, dict) and matrix.get('type') == 'block'):
            logger.warning("Matrix is not in block format, cannot use block elimination")
            return self._solve_direct(matrix, rhs)
        
        # Extract block information
        blocks = matrix['blocks']
        block_pattern = matrix['block_pattern']
        block_sizes = matrix['block_sizes']
        
        # Solve using block operations
        try:
            solution = self.block_ops.solve_block_system(blocks, block_pattern, block_sizes, rhs)
            
            info = {
                'success': True,
                'method': 'block_elimination'
            }
        except Exception as e:
            logger.warning(f"Block elimination failed: {str(e)}")
            
            # Fallback to direct solver
            solution, info = self._solve_direct(matrix, rhs)
            info['method'] = 'direct (block elimination fallback)'
        
        return solution, info


class BlockEliminationSolver:
    """
    Specialized solver using block elimination for the MISES structure.
    
    This class implements a tailored block elimination method for the
    specific structure of the MISES equations, which have both local
    variables (defined on the grid) and global variables (such as
    mass fluxes).
    
    Attributes:
        config: Configuration dictionary with solver parameters.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize block elimination solver.
        
        Args:
            config: Configuration dictionary with solver parameters:
                - matrix_format: Format for matrix storage ('dense', 'sparse')
                - direct_solver: Type of direct solver for diagonal blocks
        """
        self.config = config or {}
        
        # Solver parameters
        self.matrix_format = self.config.get('matrix_format', 'sparse')
        self.direct_solver = self.config.get('direct_solver', 'lu')
        
        logger.info(f"Initialized block elimination solver: format={self.matrix_format}")
    
    def solve(self, blocks: Dict[str, Any], rhs: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """
        Solve the MISES system using specialized block elimination.
        
        This method implements the block elimination approach described in
        Drela's MISES papers, which handles the specific structure of the
        coupled Euler and boundary layer equations.
        
        Args:
            blocks: Dictionary with matrix blocks.
            rhs: Dictionary with right-hand side vectors.
            
        Returns:
            Dictionary with solution components.
        """
        # Extract blocks and RHS components
        A = blocks.get('A')  # Euler equation Jacobian
        B = blocks.get('B')  # Coupling: Euler -> BL
        C = blocks.get('C')  # Coupling: BL -> Euler
        Z = blocks.get('Z')  # BL equation Jacobian
        
        G = blocks.get('G')  # Euler -> Global
        H = blocks.get('H')  # Global -> Euler
        K = blocks.get('K')  # Global -> BL
        L = blocks.get('L')  # BL -> Global
        M = blocks.get('M')  # Global equation Jacobian
        
        # Extract RHS components
        r_euler = rhs.get('euler')
        r_bl = rhs.get('bl')
        r_global = rhs.get('global')
        
        # Check if all required blocks are present
        required_blocks = ['A', 'Z']
        missing_blocks = [block for block in required_blocks if block not in blocks]
        
        if missing_blocks:
            logger.error(f"Missing required blocks: {missing_blocks}")
            return {}
        
        # Check if global variables are present
        has_global = all(block in blocks for block in ['G', 'H', 'M'])
        
        # Initialize solution components
        solution = {}
        
        # Step 1: Eliminate Z block to get Schur complement for A
        # This is done by solving Z * delta_z = r_bl - B * delta_a for delta_z
        # and substituting into the Euler equation to get
        # (A - C * Z^-1 * B) * delta_a = r_euler - C * Z^-1 * r_bl
        
        # Initialize modified system
        if self.matrix_format == 'sparse':
            # Use sparse matrix operations
            A_modified = A.copy()
            r_euler_modified = r_euler.copy()
            
            # Compute Z^-1 * r_bl
            if sp.issparse(Z):
                Z_inv_r_bl = spla.spsolve(Z, r_bl)
            else:
                Z_inv_r_bl = np.linalg.solve(Z, r_bl)
            
            # Compute Z^-1 * B
            if B is not None:
                if sp.issparse(Z) and sp.issparse(B):
                    Z_inv_B = spla.spsolve(Z, B.tocsc()).tocsr()
                else:
                    Z_dense = Z.toarray() if sp.issparse(Z) else Z
                    B_dense = B.toarray() if sp.issparse(B) else B
                    Z_inv_B = np.linalg.solve(Z_dense, B_dense)
                
                # Compute Schur complement A_modified = A - C * Z^-1 * B
                if C is not None:
                    if sp.issparse(C) and sp.issparse(Z_inv_B):
                        A_modified = A - C * Z_inv_B
                    else:
                        C_dense = C.toarray() if sp.issparse(C) else C
                        Z_inv_B_dense = Z_inv_B.toarray() if sp.issparse(Z_inv_B) else Z_inv_B
                        A_modified = A - sp.csr_matrix(C_dense @ Z_inv_B_dense)
            
            # Compute modified RHS r_euler_modified = r_euler - C * Z^-1 * r_bl
            if C is not None:
                if sp.issparse(C):
                    r_euler_modified = r_euler - C * Z_inv_r_bl
                else:
                    r_euler_modified = r_euler - C @ Z_inv_r_bl
        
        else:
            # Use dense matrix operations
            A_dense = A.toarray() if sp.issparse(A) else A
            Z_dense = Z.toarray() if sp.issparse(Z) else Z
            
            A_modified = np.copy(A_dense)
            r_euler_modified = np.copy(r_euler)
            
            # Compute Z^-1 * r_bl
            Z_inv_r_bl = np.linalg.solve(Z_dense, r_bl)
            
            # Compute Schur complement and modified RHS
            if B is not None and C is not None:
                B_dense = B.toarray() if sp.issparse(B) else B
                C_dense = C.toarray() if sp.issparse(C) else C
                
                Z_inv_B = np.linalg.solve(Z_dense, B_dense)
                A_modified = A_dense - C_dense @ Z_inv_B
                r_euler_modified = r_euler - C_dense @ Z_inv_r_bl
        
        # Step 2: Handle global variables if present
        if has_global:
            # Eliminate G and L blocks to get Schur complement for M
            # (M - L * Z^-1 * K - H * A_modified^-1 * G) * delta_g = ...
            if self.matrix_format == 'sparse':
                # Global variable handling with sparse matrices
                # This is a simplified implementation
                # Full implementation would be more complex
                
                # Compute A_modified^-1 * G
                if sp.issparse(A_modified) and sp.issparse(G):
                    A_mod_inv_G = spla.spsolve(A_modified, G.tocsc()).tocsr()
                else:
                    A_mod_dense = A_modified.toarray() if sp.issparse(A_modified) else A_modified
                    G_dense = G.toarray() if sp.issparse(G) else G
                    A_mod_inv_G = np.linalg.solve(A_mod_dense, G_dense)
                
                # Compute Schur complement for M
                M_modified = M
                
                if H is not None:
                    if sp.issparse(H) and sp.issparse(A_mod_inv_G):
                        M_modified = M - H * A_mod_inv_G
                    else:
                        H_dense = H.toarray() if sp.issparse(H) else H
                        A_mod_inv_G_dense = A_mod_inv_G.toarray() if sp.issparse(A_mod_inv_G) else A_mod_inv_G
                        M_modified = M - sp.csr_matrix(H_dense @ A_mod_inv_G_dense)
                
                if L is not None and K is not None:
                    # Add contribution from BL variables
                    if sp.issparse(Z) and sp.issparse(K):
                        Z_inv_K = spla.spsolve(Z, K.tocsc()).tocsr()
                    else:
                        Z_dense = Z.toarray() if sp.issparse(Z) else Z
                        K_dense = K.toarray() if sp.issparse(K) else K
                        Z_inv_K = np.linalg.solve(Z_dense, K_dense)
                    
                    if sp.issparse(L) and sp.issparse(Z_inv_K):
                        M_modified = M_modified - L * Z_inv_K
                    else:
                        L_dense = L.toarray() if sp.issparse(L) else L
                        Z_inv_K_dense = Z_inv_K.toarray() if sp.issparse(Z_inv_K) else Z_inv_K
                        M_modified = M_modified - sp.csr_matrix(L_dense @ Z_inv_K_dense)
                
                # Compute modified RHS for global equations
                r_global_modified = r_global
                
                # Contribution from Euler variables
                if H is not None:
                    # Compute A_modified^-1 * r_euler_modified
                    if sp.issparse(A_modified):
                        A_mod_inv_r_euler = spla.spsolve(A_modified, r_euler_modified)
                    else:
                        A_mod_inv_r_euler = np.linalg.solve(A_modified, r_euler_modified)
                    
                    if sp.issparse(H):
                        r_global_modified = r_global - H * A_mod_inv_r_euler
                    else:
                        r_global_modified = r_global - H @ A_mod_inv_r_euler
                
                # Contribution from BL variables
                if L is not None:
                    if sp.issparse(L):
                        r_global_modified = r_global_modified - L * Z_inv_r_bl
                    else:
                        r_global_modified = r_global_modified - L @ Z_inv_r_bl
                
                # Solve for global variables
                if sp.issparse(M_modified):
                    delta_g = spla.spsolve(M_modified, r_global_modified)
                else:
                    delta_g = np.linalg.solve(M_modified, r_global_modified)
                
                # Store global solution
                solution['global'] = delta_g
                
                # Update Euler RHS with global solution
                if G is not None:
                    if sp.issparse(G):
                        r_euler_modified = r_euler_modified - G * delta_g
                    else:
                        r_euler_modified = r_euler_modified - G @ delta_g
                
                # Update BL RHS with global solution
                if K is not None:
                    if sp.issparse(K):
                        r_bl_modified = r_bl - K * delta_g
                    else:
                        r_bl_modified = r_bl - K @ delta_g
                else:
                    r_bl_modified = r_bl
            
            else:
                # Global variable handling with dense matrices
                # Similar to above but using dense operations
                # This is a simplified implementation
                
                # Convert to dense if needed
                A_mod_dense = A_modified.toarray() if sp.issparse(A_modified) else A_modified
                G_dense = G.toarray() if sp.issparse(G) else G
                H_dense = H.toarray() if sp.issparse(H) else H
                M_dense = M.toarray() if sp.issparse(M) else M
                
                # Compute A_modified^-1 * G
                A_mod_inv_G = np.linalg.solve(A_mod_dense, G_dense)
                
                # Compute Schur complement for M
                M_modified = M_dense - H_dense @ A_mod_inv_G
                
                if L is not None and K is not None:
                    # Add contribution from BL variables
                    Z_dense = Z.toarray() if sp.issparse(Z) else Z
                    K_dense = K.toarray() if sp.issparse(K) else K
                    L_dense = L.toarray() if sp.issparse(L) else L
                    
                    Z_inv_K = np.linalg.solve(Z_dense, K_dense)
                    M_modified = M_modified - L_dense @ Z_inv_K
                
                # Compute modified RHS for global equations
                r_global_modified = r_global
                
                # Contribution from Euler variables
                A_mod_inv_r_euler = np.linalg.solve(A_mod_dense, r_euler_modified)
                r_global_modified = r_global - H_dense @ A_mod_inv_r_euler
                
                # Contribution from BL variables
                if L is not None:
                    r_global_modified = r_global_modified - L_dense @ Z_inv_r_bl
                
                # Solve for global variables
                delta_g = np.linalg.solve(M_modified, r_global_modified)
                
                # Store global solution
                solution['global'] = delta_g
                
                # Update Euler RHS with global solution
                r_euler_modified = r_euler_modified - G_dense @ delta_g
                
                # Update BL RHS with global solution
                if K is not None:
                    r_bl_modified = r_bl - K_dense @ delta_g
                else:
                    r_bl_modified = r_bl
        
        else:
            # No global variables, use original RHS
            r_bl_modified = r_bl
        
        # Step 3: Solve for Euler variables
        if self.matrix_format == 'sparse':
            if sp.issparse(A_modified):
                delta_a = spla.spsolve(A_modified, r_euler_modified)
            else:
                delta_a = np.linalg.solve(A_modified, r_euler_modified)
        else:
            delta_a = np.linalg.solve(A_modified, r_euler_modified)
        
        # Store Euler solution
        solution['euler'] = delta_a
        
        # Step 4: Solve for BL variables
        # Back-substitute into BL equations
        if B is not None:
            if self.matrix_format == 'sparse':
                if sp.issparse(B):
                    r_bl_modified = r_bl_modified - B * delta_a
                else:
                    r_bl_modified = r_bl_modified - B @ delta_a
            else:
                B_dense = B.toarray() if sp.issparse(B) else B
                r_bl_modified = r_bl_modified - B_dense @ delta_a
        
        if self.matrix_format == 'sparse':
            if sp.issparse(Z):
                delta_z = spla.spsolve(Z, r_bl_modified)
            else:
                delta_z = np.linalg.solve(Z, r_bl_modified)
        else:
            Z_dense = Z.toarray() if sp.issparse(Z) else Z
            delta_z = np.linalg.solve(Z_dense, r_bl_modified)
        
        # Store BL solution
        solution['bl'] = delta_z
        
        return solution


def solve_tridiagonal_system(a: np.ndarray, b: np.ndarray, c: np.ndarray,
                           d: np.ndarray) -> np.ndarray:
    """
    Solve a tridiagonal system using the Thomas algorithm.
    
    This function solves a system of the form:
    a_i * x_{i-1} + b_i * x_i + c_i * x_{i+1} = d_i
    
    Args:
        a: Lower diagonal (first element is not used).
        b: Main diagonal.
        c: Upper diagonal (last element is not used).
        d: Right-hand side.
        
    Returns:
        Solution vector.
    """
    n = len(d)
    
    # Initialize solution vector
    x = np.zeros(n)
    
    # Forward elimination
    c_prime = np.zeros(n)
    d_prime = np.zeros(n)
    
    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]
    
    for i in range(1, n):
        denominator = b[i] - a[i] * c_prime[i-1]
        if abs(denominator) < 1e-12:
            # Avoid division by zero
            denominator = 1e-12 if denominator >= 0 else -1e-12
        
        c_prime[i] = c[i] / denominator if i < n-1 else 0
        d_prime[i] = (d[i] - a[i] * d_prime[i-1]) / denominator
    
    # Back substitution
    x[n-1] = d_prime[n-1]
    
    for i in range(n-2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i+1]
    
    return x


def solve_block_tridiagonal_system(A: List[np.ndarray], B: List[np.ndarray],
                                 C: List[np.ndarray], D: List[np.ndarray]) -> List[np.ndarray]:
    """
    Solve a block tridiagonal system using block Thomas algorithm.
    
    This function solves a system of the form:
    A_i * X_{i-1} + B_i * X_i + C_i * X_{i+1} = D_i
    where A_i, B_i, C_i are matrices and X_i, D_i are vectors.
    
    Args:
        A: List of lower diagonal block matrices (first element is not used).
        B: List of main diagonal block matrices.
        C: List of upper diagonal block matrices (last element is not used).
        D: List of right-hand side vectors.
        
    Returns:
        List of solution vectors.
    """
    n = len(D)
    
    # Initialize solution vectors
    X = [np.zeros_like(D[i]) for i in range(n)]
    
    # Forward elimination
    C_prime = [None] * n
    D_prime = [None] * n
    
    # First block - direct inverse
    if isinstance(B[0], np.ndarray):
        B_inv = np.linalg.inv(B[0])
        C_prime[0] = B_inv @ C[0]
        D_prime[0] = B_inv @ D[0]
    else:
        # Assume B[0] is a scalar
        C_prime[0] = C[0] / B[0]
        D_prime[0] = D[0] / B[0]
    
    # Forward sweep
    for i in range(1, n):
        # Calculate B_i - A_i * C'_{i-1}
        if isinstance(B[i], np.ndarray):
            # Matrix operations
            modified_B = B[i] - A[i] @ C_prime[i-1]
            B_inv = np.linalg.inv(modified_B)
            
            # Calculate C'_i = B_inv * C_i
            C_prime[i] = B_inv @ C[i] if i < n-1 else np.zeros_like(B[i])
            
            # Calculate D'_i = B_inv * (D_i - A_i * D'_{i-1})
            D_prime[i] = B_inv @ (D[i] - A[i] @ D_prime[i-1])
        else:
            # Scalar operations
            modified_B = B[i] - A[i] * C_prime[i-1]
            
            # Avoid division by zero
            if abs(modified_B) < 1e-12:
                modified_B = 1e-12 if modified_B >= 0 else -1e-12
            
            C_prime[i] = C[i] / modified_B if i < n-1 else 0
            D_prime[i] = (D[i] - A[i] * D_prime[i-1]) / modified_B
    
    # Back substitution
    X[n-1] = D_prime[n-1]
    
    for i in range(n-2, -1, -1):
        if isinstance(C_prime[i], np.ndarray):
            X[i] = D_prime[i] - C_prime[i] @ X[i+1]
        else:
            X[i] = D_prime[i] - C_prime[i] * X[i+1]
    
    return X


class IterativeSolver:
    """
    Iterative solver for sparse linear systems.
    
    This class provides iterative methods for solving sparse linear systems,
    including GMRES, BiCGSTAB, and other Krylov subspace methods.
    
    Attributes:
        config: Configuration dictionary with solver parameters.
        method: Iterative method to use ('gmres', 'bicgstab', etc.).
        tolerance: Convergence tolerance.
        max_iterations: Maximum number of iterations.
        use_preconditioner: Whether to use preconditioning.
        preconditioner: Preconditioner object or callable.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize iterative solver.
        
        Args:
            config: Configuration dictionary with solver parameters:
                - method: Iterative method ('gmres', 'bicgstab', 'cg')
                - tolerance: Convergence tolerance
                - max_iterations: Maximum iterations
                - use_preconditioner: Whether to use preconditioning
                - preconditioner_type: Type of preconditioner to use
        """
        self.config = config or {}
        
        # Solver parameters
        self.method = self.config.get('method', 'gmres')
        self.tolerance = self.config.get('tolerance', 1e-8)
        self.max_iterations = self.config.get('max_iterations', 100)
        self.use_preconditioner = self.config.get('use_preconditioner', True)
        
        # Initialize preconditioner
        self.preconditioner = None
        if self.use_preconditioner:
            self.preconditioner = self.setup_preconditioner()
        
        logger.info(f"Initialized iterative solver: method={self.method}, "
                   f"tolerance={self.tolerance}, max_iterations={self.max_iterations}")
    
    def setup_preconditioner(self) -> Callable:
        """
        Set up preconditioner based on configuration.
        
        Returns:
            Preconditioner function or object.
        """
        precond_type = self.config.get('preconditioner_type', 'ilu')
        
        if precond_type == 'ilu':
            def ilu_preconditioner(A, b=None):
                try:
                    ilu = spla.spilu(A.tocsc(), drop_tol=1e-4)
                    M = spla.LinearOperator(A.shape, lambda x: ilu.solve(x))
                    return M
                except Exception as e:
                    logger.warning(f"Failed to create ILU preconditioner: {str(e)}")
                    return None
            
            return ilu_preconditioner
        
        elif precond_type == 'diagonal':
            def diag_preconditioner(A, b=None):
                try:
                    diag = A.diagonal()
                    diag_inv = 1.0 / np.maximum(np.abs(diag), 1e-10)
                    M = spla.LinearOperator(A.shape, lambda x: diag_inv * x)
                    return M
                except Exception as e:
                    logger.warning(f"Failed to create diagonal preconditioner: {str(e)}")
                    return None
            
            return diag_preconditioner
        
        elif precond_type == 'custom':
            # Custom preconditioner provided in config
            custom_precond = self.config.get('custom_preconditioner')
            if callable(custom_precond):
                return custom_precond
        
        # Default to None if no valid preconditioner
        logger.warning(f"Unknown preconditioner type: {precond_type}, using no preconditioning")
        return None
    
    def solve(self, A: Union[np.ndarray, sp.spmatrix], b: np.ndarray) -> Tuple[np.ndarray, Dict[str, Any]]:
        """
        Solve linear system Ax = b using iterative methods.
        
        Args:
            A: System matrix (dense or sparse).
            b: Right-hand side vector.
            
        Returns:
            Tuple of (solution, info_dict).
        """
        # Ensure matrix is in sparse format
        if not sp.issparse(A):
            A = sp.csr_matrix(A)
        
        # Start timer
        start_time = time.time()
        
        # Apply preconditioner if requested
        M = None
        if self.use_preconditioner and self.preconditioner is not None:
            try:
                M = self.preconditioner(A, b)
            except Exception as e:
                logger.warning(f"Failed to apply preconditioner: {str(e)}")
        
        # Select and apply iterative method
        if self.method == 'gmres':
            try:
                solution, info = spla.gmres(A, b, M=M, tol=self.tolerance, atol=1e-10, 
                                           maxiter=self.max_iterations, restart=min(50, len(b)))
                success = (info == 0)
            except Exception as e:
                logger.error(f"GMRES solver failed: {str(e)}")
                solution = np.zeros_like(b)
                info = -1
                success = False
        
        elif self.method == 'bicgstab':
            try:
                solution, info = spla.bicgstab(A, b, M=M, tol=self.tolerance, atol=1e-10, 
                                              maxiter=self.max_iterations)
                success = (info == 0)
            except Exception as e:
                logger.error(f"BiCGSTAB solver failed: {str(e)}")
                solution = np.zeros_like(b)
                info = -1
                success = False
        
        elif self.method == 'cg':
            # Check if matrix is approximately symmetric
            if not self.is_approximately_symmetric(A):
                logger.warning("CG method used with non-symmetric matrix")
            
            try:
                solution, info = spla.cg(A, b, M=M, tol=self.tolerance, atol=1e-10, 
                                        maxiter=self.max_iterations)
                success = (info == 0)
            except Exception as e:
                logger.error(f"CG solver failed: {str(e)}")
                solution = np.zeros_like(b)
                info = -1
                success = False
        
        else:
            logger.error(f"Unknown iterative method: {self.method}")
            solution = np.zeros_like(b)
            info = -1
            success = False
        
        # Calculate solution time
        solution_time = time.time() - start_time
        
        # Prepare info dictionary
        info_dict = {
            'success': success,
            'iterations': info if success else self.max_iterations,
            'method': self.method,
            'solution_time': solution_time,
            'error': None if success else f"Solver failed to converge"
        }
        
        return solution, info_dict
    
    @staticmethod
    def is_approximately_symmetric(A: sp.spmatrix, tol: float = 1e-8) -> bool:
        """
        Check if a matrix is approximately symmetric.
        
        Args:
            A: Input matrix.
            tol: Tolerance for symmetry check.
            
        Returns:
            True if matrix is approximately symmetric, False otherwise.
        """
        # For large matrices, check only a sample of entries
        if A.shape[0] > 1000:
            # Take a sample of rows/columns
            sample_size = min(100, A.shape[0])
            indices = np.random.choice(A.shape[0], size=sample_size, replace=False)
            
            # Check symmetry on the sampled submatrix
            for i in indices:
                for j in indices:
                    if abs(A[i, j] - A[j, i]) > tol * (abs(A[i, j]) + abs(A[j, i]) + 1e-10):
                        return False
            
            return True
        else:
            # For smaller matrices, check all entries
            diff = A - A.T
            return np.max(np.abs(diff.data)) <= tol * (np.max(np.abs(A.data)) + 1e-10)