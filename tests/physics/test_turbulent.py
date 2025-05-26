"""
Tests for the turbulent closure models in pymises.physics.turbulent.
"""

import os
import sys
import numpy as np
import unittest

from pymises.physics.turbulent import TurbulentClosure

class TestTurbulentClosure(unittest.TestCase):
    """Test the turbulent closure relations."""
    
    def setUp(self):
        """Set up the test case."""
        self.closure = TurbulentClosure({
            'lag_constant': 5.6
        })
    
    def test_kinematic_shape_parameter(self):
        """Test kinematic shape parameter calculation."""
        # For incompressible flow, Hk = H
        H_incomp = 1.5
        Hk_incomp = self.closure.kinematic_shape_parameter(H=H_incomp, Mach=0.0)
        self.assertAlmostEqual(Hk_incomp, H_incomp, places=6)
        
        # Compressibility effect should reduce Hk at higher Mach
        H_comp = 1.5
        Hk_comp = self.closure.kinematic_shape_parameter(H=H_comp, Mach=0.8)
        self.assertLess(Hk_comp, H_comp, 
                       "Compressibility should reduce kinematic shape parameter")
    
    def test_energy_shape_parameter(self):
        """Test energy shape parameter calculation."""
        # Test at different Hk values
        H_star_1 = self.closure.energy_shape_parameter(Hk=1.5, Re_theta=1000.0)
        H_star_2 = self.closure.energy_shape_parameter(Hk=2.0, Re_theta=1000.0)
        
        # H* should decrease as Hk increases
        self.assertGreater(H_star_1, H_star_2, 
                          "Energy shape parameter should decrease with increasing Hk")
        
        # Reynolds number correction should be more significant at low Re
        H_star_low_re = self.closure.energy_shape_parameter(Hk=1.5, Re_theta=100.0)
        self.assertGreater(H_star_low_re, H_star_1, 
                          "Energy shape parameter should be higher at low Reynolds number")
    
    def test_skin_friction_coefficient(self):
        """Test skin friction coefficient calculation."""
        # Test at different Hk values and Reynolds numbers
        Cf_1 = self.closure.skin_friction_coefficient(Hk=1.5, Re_theta=1000.0, Mach=0.0)
        Cf_2 = self.closure.skin_friction_coefficient(Hk=2.0, Re_theta=1000.0, Mach=0.0)
        
        # Cf should decrease as Hk increases (more adverse pressure gradient)
        self.assertGreater(Cf_1, Cf_2, 
                          "Skin friction should decrease with increasing Hk")
        
        # Cf should decrease as Reynolds number increases
        Cf_high_re = self.closure.skin_friction_coefficient(Hk=1.5, Re_theta=10000.0, Mach=0.0)
        self.assertGreater(Cf_1, Cf_high_re, 
                          "Skin friction should decrease with increasing Reynolds number")
        
        # Compressibility effect should reduce Cf at higher Mach
        Cf_comp = self.closure.skin_friction_coefficient(Hk=1.5, Re_theta=1000.0, Mach=0.8)
        self.assertLess(Cf_comp, Cf_1, 
                       "Compressibility should reduce skin friction")
    
    def test_equilibrium_dissipation_coefficient(self):
        """Test equilibrium dissipation coefficient calculation."""
        # Test at different Hk values
        Cf = 0.005  # Reference skin friction
        CD_eq_1 = self.closure.equilibrium_dissipation_coefficient(Hk=1.5, Cf=Cf)
        CD_eq_2 = self.closure.equilibrium_dissipation_coefficient(Hk=2.0, Cf=Cf)
        
        # CD_eq should increase as Hk increases (further from equilibrium)
        self.assertLess(CD_eq_1, CD_eq_2, 
                        "Equilibrium dissipation coefficient should increase with Hk")
        
        # CD_eq should be proportional to Cf
        CD_eq_low_cf = self.closure.equilibrium_dissipation_coefficient(Hk=1.5, Cf=0.001)
        self.assertAlmostEqual(CD_eq_low_cf/0.001, CD_eq_1/0.005, places=2, msg="Equilibrium dissipation should scale approximately with Cf")
    
    def test_shear_stress_coefficient(self):
        """Test initial shear stress coefficient at transition."""
        H = 2.5  # Shape parameter
        Cf = 0.003  # Skin friction
        CT_eq = 0.01  # Equilibrium shear stress
        
        CT = self.closure.shear_stress_coefficient(H, Cf, CT_eq)
        
        # CT should be positive
        self.assertGreater(CT, 0.0, "Shear stress coefficient should be positive")
        
        # For higher H (stronger adverse gradient), initial CT should be closer to equilibrium
        H_high = 3.0
        CT_high_H = self.closure.shear_stress_coefficient(H_high, Cf, CT_eq)
        
        # The ratio to equilibrium should be higher for higher H
        self.assertGreater(CT_high_H/CT_eq, CT/CT_eq, 
                          "Initial CT ratio should be higher for adverse pressure gradient")
    
    def test_update_shear_stress(self):
        """Test Reynolds stress lag equation."""
        CT = 0.005  # Current shear stress
        CT_eq = 0.010  # Equilibrium shear stress (higher than current)
        dx = 0.01  # Streamwise step
        delta = 0.05  # Boundary layer thickness
        
        # Update CT using lag equation
        CT_new = self.closure.update_shear_stress(CT, CT_eq, dx, delta)
        
        # If CT_eq > CT, then CT_new should increase
        self.assertGreater(CT_new, CT, "CT should increase when below equilibrium")
        
        # If CT_eq < CT, then CT_new should decrease
        CT_high = 0.015
        CT_new_high = self.closure.update_shear_stress(CT_high, CT_eq, dx, delta)
        self.assertLess(CT_new_high, CT_high, "CT should decrease when above equilibrium")
    
    def test_compute_closure_relations(self):
        """Test computing all closure relations in one call."""
        # Test at typical boundary layer conditions
        H = 1.6
        Re_theta = 1000.0
        Mach = 0.3
        
        # Compute all relations
        results = self.closure.compute_closure_relations(H, Re_theta, Mach)
        
        # Check if all expected keys are present
        expected_keys = ['Hk', 'H_star', 'Cf', 'CD_eq', 'CT_eq', 'Us', 'H_star_star']
        for key in expected_keys:
            self.assertIn(key, results, f"Missing key: {key}")
        
        # Check if values are in reasonable ranges
        self.assertGreater(results['Cf'], 0.0, "Skin friction should be positive")
        self.assertGreater(results['CD_eq'], 0.0, "Dissipation coefficient should be positive")
        self.assertGreater(results['H_star'], 1.0, "Energy shape parameter should be > 1.0")
        
        # Test with CT provided
        CT = 0.005
        results_with_CT = self.closure.compute_closure_relations(H, Re_theta, Mach, CT)
        
        # Should include CD with CT provided
        self.assertIn('CD', results_with_CT, "CD should be included when CT is provided")
        self.assertIn('CT', results_with_CT, "CT should be included in results")
        
        # CD is a combination of wall term and wake term
        self.assertGreater(results_with_CT['CD'], 0.0, "Dissipation coefficient should be positive")

if __name__ == '__main__':
    unittest.main()
