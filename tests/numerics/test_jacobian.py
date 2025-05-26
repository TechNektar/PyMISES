"""
Tests for the Jacobian module in pymises.numerics.jacobian.
"""

import os
import sys
import numpy as np
import unittest
import scipy.sparse as sp

from pymises.numerics.jacobian import JacobianAssembler, JacobianCalculator
from pymises.utils.validation import validate_jacobian

class TestJacobianCalculator(unittest.TestCase):
    """Test the JacobianCalculator class."""
    
    def setUp(self):
        """Set up test cases for the Jacobian calculator."""
        # Define a simple nonlinear function: F(x) = [x1^2 + x2, x1*x2]
        def function(x):
            return np.array([x[0]**2 + x[1], x[0] * x[1]])
        
        self.function = function
        self.x0 = np.array([2.0, 3.0])  # Test point
        
        # Create a Jacobian calculator
        self.calculator = JacobianCalculator(method='finite_difference')
    
    def test_finite_difference_jacobian(self):
        """Test finite difference Jacobian calculation."""
        # Calculate Jacobian using finite differences
        J_fd = self.calculator.calculate(self.function, self.x0)
        
        # Expected Jacobian at x0 = [2, 3]:
        # J = [2*x1, 1]
        #     [x2,   x1]
        # J = [4, 1]
        #     [3, 2]
        J_expected = np.array([
            [4.0, 1.0],
            [3.0, 2.0]
        ])
        
        # Check Jacobian values
        np.testing.assert_allclose(J_fd.toarray(), J_expected, rtol=1e-5)
    
    def test_complex_step_jacobian(self):
        """Test complex step Jacobian calculation."""
        # Create a calculator using complex step method
        calculator = JacobianCalculator(method='complex_step')
        
        # Define a function that can handle complex inputs
        def complex_function(x):
            return np.array([x[0]**2 + x[1], x[0] * x[1]])
        
        # Calculate Jacobian using complex step
        J_cs = calculator.calculate(complex_function, self.x0)
        
        # Expected Jacobian at x0 = [2, 3]:
        J_expected = np.array([
            [4.0, 1.0],
            [3.0, 2.0]
        ])
        
        # Check Jacobian values
        np.testing.assert_allclose(J_cs.toarray(), J_expected, rtol=1e-10)
    
    def test_analytical_jacobian(self):
        """Test analytical Jacobian calculation."""
        # Define analytical Jacobian function
        def jacobian_function(x):
            J = np.array([
                [2*x[0], 1.0],
                [x[1], x[0]]
            ])
            return sp.csr_matrix(J)
        
        # Create a calculator using analytical method
        calculator = JacobianCalculator(method='analytical', jacobian_function=jacobian_function)
        
        # Calculate Jacobian
        J_analytical = calculator.calculate(self.function, self.x0)
        
        # Expected Jacobian at x0 = [2, 3]:
        J_expected = np.array([
            [4.0, 1.0],
            [3.0, 2.0]
        ])
        
        # Check Jacobian values
        np.testing.assert_allclose(J_analytical.toarray(), J_expected, rtol=1e-10)
    
    def test_step_size_effect(self):
        """Test the effect of step size on finite difference accuracy."""
        # Create calculators with different step sizes
        calc1 = JacobianCalculator(method='finite_difference', step_size=1e-4)
        calc2 = JacobianCalculator(method='finite_difference', step_size=1e-8)
        
        # Calculate Jacobians
        J1 = calc1.calculate(self.function, self.x0)
        J2 = calc2.calculate(self.function, self.x0)
        
        # Expected Jacobian
        J_expected = np.array([
            [4.0, 1.0],
            [3.0, 2.0]
        ])
        
        # Calculate errors
        error1 = np.max(np.abs(J1.toarray() - J_expected))
        error2 = np.max(np.abs(J2.toarray() - J_expected))
        
        # Smaller step size should give better accuracy up to a point
        # But very small step sizes might cause numerical issues
        # For this simple function, we expect reasonable accuracy with either step size
        self.assertLess(error1, 1e-3)
        self.assertLess(error2, 1e-3)

class TestJacobianAssembler(unittest.TestCase):
    """Test the JacobianAssembler class."""
    
    def setUp(self):
        """Set up test cases for the Jacobian assembler."""
        # Create a JacobianAssembler for a simple 2D system
        self.n_variables = 2
        self.n_equations = 2
        self.assembler = JacobianAssembler(self.n_variables, self.n_equations)
    
    def test_assembly_from_components(self):
        """Test assembling a Jacobian from component matrices."""
        # Create component matrices
        A = np.array([[1.0, 2.0], [3.0, 4.0]])
        B = np.array([[5.0, 6.0], [7.0, 8.0]])
        C = np.array([[9.0, 10.0], [11.0, 12.0]])
        D = np.array([[13.0, 14.0], [15.0, 16.0]])
        
        # Add components to the assembler
        self.assembler.add_component(0, 0, sp.csr_matrix(A))
        self.assembler.add_component(0, 1, sp.csr_matrix(B))
        self.assembler.add_component(1, 0, sp.csr_matrix(C))
        self.assembler.add_component(1, 1, sp.csr_matrix(D))
        
        # Assemble the full Jacobian
        J = self.assembler.assemble()
        
        # Expected full matrix
        J_expected = np.array([
            [1.0, 2.0, 5.0, 6.0],
            [3.0, 4.0, 7.0, 8.0],
            [9.0, 10.0, 13.0, 14.0],
            [11.0, 12.0, 15.0, 16.0]
        ])
        
        # Check assembled Jacobian
        np.testing.assert_allclose(J.toarray(), J_expected)
    
    def test_empty_components(self):
        """Test assembly with some empty components."""
        # Create component matrices
        A = np.array([[1.0, 2.0], [3.0, 4.0]])
        D = np.array([[13.0, 14.0], [15.0, 16.0]])
        
        # Add only some components (B and C are zero)
        self.assembler.add_component(0, 0, sp.csr_matrix(A))
        self.assembler.add_component(1, 1, sp.csr_matrix(D))
        
        # Assemble the full Jacobian
        J = self.assembler.assemble()
        
        # Expected full matrix
        J_expected = np.array([
            [1.0, 2.0, 0.0, 0.0],
            [3.0, 4.0, 0.0, 0.0],
            [0.0, 0.0, 13.0, 14.0],
            [0.0, 0.0, 15.0, 16.0]
        ])
        
        # Check assembled Jacobian
        np.testing.assert_allclose(J.toarray(), J_expected)
    
    def test_sparse_component_addition(self):
        """Test adding sparse components."""
        # Create sparse component matrices
        rows = [0, 1]
        cols = [0, 1]
        data = [1.0, 4.0]
        A = sp.csr_matrix((data, (rows, cols)), shape=(2, 2))
        
        rows = [0, 1]
        cols = [0, 1]
        data = [13.0, 16.0]
        D = sp.csr_matrix((data, (rows, cols)), shape=(2, 2))
        
        # Add sparse components
        self.assembler.add_component(0, 0, A)
        self.assembler.add_component(1, 1, D)
        
        # Assemble the full Jacobian
        J = self.assembler.assemble()
        
        # Expected full matrix
        J_expected = np.array([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 4.0, 0.0, 0.0],
            [0.0, 0.0, 13.0, 0.0],
            [0.0, 0.0, 0.0, 16.0]
        ])
        
        # Check assembled Jacobian
        np.testing.assert_allclose(J.toarray(), J_expected)
    
    def test_incorrect_component_shape(self):
        """Test error handling for components with incorrect shape."""
        # Create component with wrong shape
        A = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])  # 2x3, not 2x2
        
        # Adding component with wrong shape should raise ValueError
        with self.assertRaises(ValueError):
            self.assembler.add_component(0, 0, sp.csr_matrix(A))

class TestJacobianValidation(unittest.TestCase):
    """Test validation of Jacobian calculations."""
    
    def test_jacobian_validation(self):
        """Test the validate_jacobian function."""
        # Define a simple function and its Jacobian
        def function(x):
            return np.array([x[0]**2 + x[1], x[0] * x[1]])
        
        def jacobian_function(x):
            J = np.array([
                [2*x[0], 1.0],
                [x[1], x[0]]
            ])
            return J
        
        # Test point
        x0 = np.array([2.0, 3.0])
        
        # Validate the Jacobian
        is_valid, max_error = validate_jacobian(
            function, 
            jacobian_function,
            x0,
            h=1e-6,
            tolerance=1e-4
        )
        
        # The analytical Jacobian should match finite differences
        self.assertTrue(is_valid)
        self.assertLess(max_error, 1e-4)
        
        # Test with an incorrect Jacobian
        def incorrect_jacobian(x):
            J = np.array([
                [2*x[0], 2.0],  # Should be [2*x[0], 1.0]
                [x[1], x[0]]
            ])
            return J
        
        is_valid, max_error = validate_jacobian(
            function,
            incorrect_jacobian,
            x0,
            h=1e-6,
            tolerance=1e-4
        )
        
        # Should detect the error in the Jacobian
        self.assertFalse(is_valid)
        self.assertGreater(max_error, 1e-4)

if __name__ == '__main__':
    unittest.main()
