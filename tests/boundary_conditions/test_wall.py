"""
Tests for the wall boundary conditions in pymises.boundary_conditions.wall.
"""

import os
import sys
import numpy as np
import unittest

from pymises.boundary_conditions.wall import WallBoundaryCondition, InviscidWallBC, ViscousWallBC

class TestInviscidWallBC(unittest.TestCase):
    """Test the InviscidWallBC class."""
    
    def setUp(self):
        """Set up test case for inviscid wall boundary condition."""
        # Simple 2D test grid coordinates
        n_points = 10
        self.grid_indices = list(range(n_points))
        
        # Create a test solution dictionary
        # For simplicity, we'll use a flow along a flat wall at y=0
        self.solution = {
            'grid_x': np.linspace(0, 1, n_points),
            'grid_y': np.zeros(n_points),
            'density': np.ones(n_points),
            'velocity_x': np.ones(n_points),  # Flow along x-direction
            'velocity_y': 0.1 * np.ones(n_points),  # Small non-zero y-component to test BC
            'pressure': np.ones(n_points),
            'energy': 2.5 * np.ones(n_points)
        }
        
        # Create boundary condition
        self.wall_bc = InviscidWallBC(self.grid_indices, normal_direction="inner")
    
    def test_initialization(self):
        """Test initialization of wall boundary condition."""
        self.assertEqual(self.wall_bc.grid_indices, self.grid_indices)
        self.assertEqual(self.wall_bc.normal_direction, "inner")
    
    def test_input_validation(self):
        """Test input validation for wall boundary condition."""
        # Test invalid grid indices
        with self.assertRaises(ValueError):
            InviscidWallBC("not_a_list", normal_direction="inner")
        
        # Test invalid normal direction
        with self.assertRaises(ValueError):
            InviscidWallBC(self.grid_indices, normal_direction="invalid")
    
    def test_apply_flat_wall(self):
        """Test application of BC to a flat wall."""
        # Apply boundary condition
        updated_solution = self.wall_bc.apply(self.solution)
        
        # Check that velocity_y becomes zero at the wall (flat wall case)
        for i in self.grid_indices:
            self.assertAlmostEqual(updated_solution['velocity_y'][i], 0.0)
        
        # Velocity_x should remain unchanged
        for i in self.grid_indices:
            self.assertAlmostEqual(updated_solution['velocity_x'][i], self.solution['velocity_x'][i])
    
    def test_apply_curved_wall(self):
        """Test application of BC to a curved wall."""
        # Create a curved wall with specified tangent directions
        n_points = 10
        curved_solution = {
            'grid_x': np.linspace(0, 1, n_points),
            'grid_y': 0.1 * np.sin(np.pi * np.linspace(0, 1, n_points)),  # Curved wall
            'density': np.ones(n_points),
            'velocity_x': np.ones(n_points),
            'velocity_y': np.zeros(n_points),
            'pressure': np.ones(n_points),
            'energy': 2.5 * np.ones(n_points)
        }
        
        # Apply boundary condition
        updated_solution = self.wall_bc.apply(curved_solution)
        
        # For a curved wall, velocity should be tangent to the wall
        # Check by computing the normal vector and ensuring v·n = 0
        for i in range(1, n_points-1):
            # Tangent vector along the wall (using central differencing)
            tx = curved_solution['grid_x'][i+1] - curved_solution['grid_x'][i-1]
            ty = curved_solution['grid_y'][i+1] - curved_solution['grid_y'][i-1]
            
            # Normalize
            t_norm = np.sqrt(tx**2 + ty**2)
            tx /= t_norm
            ty /= t_norm
            
            # Normal vector (perpendicular to tangent, pointing inward)
            nx = -ty
            ny = tx
            
            # Velocity components
            vx = updated_solution['velocity_x'][i]
            vy = updated_solution['velocity_y'][i]
            
            # Dot product of velocity and normal should be close to zero
            v_dot_n = vx * nx + vy * ny
            self.assertAlmostEqual(v_dot_n, 0.0, places=5)
    
    def test_residual_contributions(self):
        """Test the residual contributions from the wall BC."""
        # Apply boundary condition
        updated_solution = self.wall_bc.apply(self.solution)
        
        # Get residual contributions
        residuals = self.wall_bc.get_residual_contributions(updated_solution)
        
        # There should be one residual per boundary point
        self.assertEqual(len(residuals), len(self.grid_indices))
        
        # For a perfect application of the BC, residuals should be near zero
        for res in residuals:
            self.assertAlmostEqual(res, 0.0, places=5)
    
    def test_jacobian_contributions(self):
        """Test the Jacobian contributions from the wall BC."""
        # Get Jacobian contributions
        rows, cols, values = self.wall_bc.get_jacobian_contributions(self.solution)
        
        # Should have non-empty Jacobian terms
        self.assertGreater(len(rows), 0)
        self.assertEqual(len(rows), len(cols))
        self.assertEqual(len(rows), len(values))

class TestViscousWallBC(unittest.TestCase):
    """Test the ViscousWallBC class."""
    
    def setUp(self):
        """Set up test case for viscous wall boundary condition."""
        # Simple 2D test grid coordinates
        n_points = 10
        self.grid_indices = list(range(n_points))
        
        # Create a test solution dictionary
        self.solution = {
            'grid_x': np.linspace(0, 1, n_points),
            'grid_y': np.zeros(n_points),
            'density': np.ones(n_points),
            'velocity_x': np.ones(n_points),
            'velocity_y': 0.1 * np.ones(n_points),
            'pressure': np.ones(n_points),
            'energy': 2.5 * np.ones(n_points)
        }
        
        # Create a displacement thickness provider function
        def displacement_thickness_provider(i):
            return 0.01 * np.sqrt(self.solution['grid_x'][i])  # Simple sqrt growth
        
        # Create boundary condition
        self.wall_bc = ViscousWallBC(
            self.grid_indices, 
            normal_direction="inner",
            displacement_thickness_provider=displacement_thickness_provider
        )
    
    def test_initialization(self):
        """Test initialization of viscous wall boundary condition."""
        self.assertEqual(self.wall_bc.grid_indices, self.grid_indices)
        self.assertEqual(self.wall_bc.normal_direction, "inner")
        self.assertIsNotNone(self.wall_bc.displacement_thickness_provider)
    
    def test_apply_with_displacement(self):
        """Test application of BC with displacement thickness."""
        # Apply boundary condition
        updated_solution = self.wall_bc.apply(self.solution)
        
        # For viscous wall, we should have a transpiration velocity
        # The value should be non-zero due to displacement thickness
        has_nonzero_v_normal = False
        for i in range(1, len(self.grid_indices)-1):
            v_normal = updated_solution['velocity_y'][i]  # For a flat wall, normal is in y-direction
            if abs(v_normal) > 1e-8:
                has_nonzero_v_normal = True
                break
                
        self.assertTrue(has_nonzero_v_normal, "No transpiration velocity found with displacement thickness")
    
    def test_zero_displacement(self):
        """Test that BC reduces to inviscid case with zero displacement."""
        # Create boundary condition with zero displacement thickness
        zero_disp_bc = ViscousWallBC(
            self.grid_indices, 
            normal_direction="inner",
            displacement_thickness_provider=lambda i: 0.0
        )
        
        # Apply boundary condition
        updated_solution = zero_disp_bc.apply(self.solution)
        
        # With zero displacement, all normal velocities should be zero
        for i in self.grid_indices:
            self.assertAlmostEqual(updated_solution['velocity_y'][i], 0.0, places=5)

class TestWallBoundaryConditionSubclassing(unittest.TestCase):
    """Test the proper subclassing of WallBoundaryCondition."""
    
    def test_abstract_methods(self):
        """Test that abstract methods must be implemented by subclasses."""
        try:
            # Create a subclass without implementing abstract methods
            class IncompleteBC(WallBoundaryCondition):
                pass
            
            # Creating an instance should fail
            IncompleteBC([0, 1, 2], normal_direction="inner")
            self.fail("Should have raised TypeError")
        except TypeError:
            # This is the expected behavior
            pass
        
        # Create a proper subclass
        class ProperBC(WallBoundaryCondition):
            def apply(self, solution):
                return solution
            
            def get_residual_contributions(self, solution):
                return np.zeros(len(self.grid_indices))
            
            def get_jacobian_contributions(self, solution):
                return [], [], []
        
        # Should be able to create an instance
        bc = ProperBC([0, 1, 2], normal_direction="inner")
        self.assertIsInstance(bc, WallBoundaryCondition)

if __name__ == '__main__':
    unittest.main()
