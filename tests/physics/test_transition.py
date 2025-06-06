"""
Tests for the transition prediction models in pymises.physics.transition.
"""

import os
import sys
import numpy as np
import unittest

from pymises.physics.transition import (
    TransitionModel,
    ModifiedAGSTransition,
    EnvelopeEnMethod,
    TransitionPredictor,
    create_transition_model
)

class TestModifiedAGSTransition(unittest.TestCase):
    """Test the Modified Abu-Ghannam/Shaw transition model."""
    
    def setUp(self):
        """Set up the test case."""
        self.transition_model = ModifiedAGSTransition({
            'amplification_constant': 0.10,
            'ramp_width': 0.30,
            'turbulence_level': 0.01
        })
    
    def test_critical_reynolds_number(self):
        """Test critical Reynolds number calculation."""
        # Test with different shape parameters and turbulence levels
        # H=2.3 is close to flat plate Blasius value
        Re_theta_crit_1 = self.transition_model.critical_reynolds_number(H=2.3, turbulence_level=0.01)
        self.assertGreater(Re_theta_crit_1, 0, "Critical Reynolds number should be positive")
        
        # H=2.8 is in adverse pressure gradient region
        Re_theta_crit_2 = self.transition_model.critical_reynolds_number(H=2.8, turbulence_level=0.01)
        self.assertGreater(Re_theta_crit_2, 0, "Critical Reynolds number should be positive")
        
        # Higher H should generally result in lower critical Reynolds number
        self.assertLess(Re_theta_crit_2, Re_theta_crit_1, 
                        "Critical Reynolds number should decrease with increasing H")
        
        # Higher turbulence level should decrease critical Reynolds number
        Re_theta_crit_high_tu = self.transition_model.critical_reynolds_number(H=2.3, turbulence_level=0.05)
        self.assertLess(Re_theta_crit_high_tu, Re_theta_crit_1, 
                        "Critical Reynolds number should decrease with higher turbulence")
    
    def test_amplification_rate(self):
        """Test amplification rate calculation."""
        # Calculate for below, at, and above critical Reynolds
        H = 2.3
        Re_theta_crit = self.transition_model.critical_reynolds_number(H)
        
        # Below critical - should have minimal bypass contribution
        rate_below = self.transition_model.amplification_rate(H, 0.8 * Re_theta_crit, Re_theta_crit)
        
        # At critical - should have modest bypass contribution
        rate_at = self.transition_model.amplification_rate(H, Re_theta_crit, Re_theta_crit)
        
        # Above critical - should have strong bypass contribution
        rate_above = self.transition_model.amplification_rate(H, 1.5 * Re_theta_crit, Re_theta_crit)
        
        # Amplification rate should increase as we go above critical Reynolds
        self.assertLess(rate_below, rate_at, 
                        "Amplification rate should increase at critical Reynolds")
        self.assertLess(rate_at, rate_above, 
                        "Amplification rate should increase above critical Reynolds")
    
    def test_update_amplification_factor(self):
        """Test amplification factor update with edge velocity correction."""
        # Initial n-factor
        n = 3.0
        
        # Case 1: Constant edge velocity (ratio = 1.0)
        n_updated1 = self.transition_model.update_amplification_factor(
            n, amplification_rate=0.1, dx=0.01, edge_velocity_ratio=1.0
        )
        # Expected: n + 0.1*0.01 = 3.001
        self.assertAlmostEqual(n_updated1, 3.001, places=6)
        
        # Case 2: Accelerating flow (ratio > 1.0) - should decrease n-factor
        n_updated2 = self.transition_model.update_amplification_factor(
            n, amplification_rate=0.1, dx=0.01, edge_velocity_ratio=1.1
        )
        # Logarithm term should produce negative contribution
        self.assertLess(n_updated2, 3.001)
        
        # Case 3: Decelerating flow (ratio < 1.0) - should increase n-factor
        n_updated3 = self.transition_model.update_amplification_factor(
            n, amplification_rate=0.1, dx=0.01, edge_velocity_ratio=0.9
        )
        # Logarithm term should produce positive contribution
        self.assertGreater(n_updated3, 3.001)

class TestTransitionPredictor(unittest.TestCase):
    """Test the TransitionPredictor class."""
    
    def setUp(self):
        """Set up the test case."""
        self.x = np.linspace(0.0, 1.0, 101)
        self.H = np.ones_like(self.x) * 2.3  # Constant shape parameter for testing
        self.Re_theta = np.linspace(100, 1000, 101)  # Increasing Reynolds number
        self.edge_velocity = np.ones_like(self.x)  # Constant edge velocity for simplicity
        
        # Create transition predictor
        self.predictor = TransitionPredictor(
            model_type='modified_ags',
            config={'turbulence_level': 0.01}
        )
    
    def test_initialization(self):
        """Test initialization of the predictor."""
        self.predictor.initialize(self.x, self.H)
        
        # Check if arrays are properly initialized
        self.assertEqual(len(self.predictor.x), len(self.x))
        self.assertEqual(len(self.predictor.n_factor), len(self.x))
        self.assertEqual(self.predictor.n_factor[0], 0.0)  # Initial n-factor should be zero
        self.assertFalse(self.predictor.transition_occurred)
        self.assertIsNone(self.predictor.transition_index)
    
    def test_predict_transition(self):
        """Test transition prediction for a Blasius boundary layer."""
        result = self.predictor.predict_transition(
            self.x, self.H, self.Re_theta, self.edge_velocity
        )
        
        # Check if transition is detected for this case
        # Note: whether transition occurs depends on the test conditions
        if result['transition_occurred']:
            self.assertGreater(result['transition_index'], 0)
            self.assertLess(result['transition_index'], len(self.x))
            self.assertIsNotNone(result['transition_x'])
        
        # Check if n_factor is calculated for all points
        self.assertEqual(len(result['n_factor']), len(self.x))
        self.assertEqual(result['n_factor'][0], 0.0)  # Initial n-factor should be zero
        
        # n-factor should be non-decreasing for this setup
        for i in range(1, len(self.x)):
            self.assertGreaterEqual(result['n_factor'][i], result['n_factor'][i-1])

class TestFactoryFunction(unittest.TestCase):
    """Test the model factory function."""
    
    def test_create_modified_ags(self):
        """Test creating a ModifiedAGSTransition model."""
        model = create_transition_model('modified_ags')
        self.assertIsInstance(model, ModifiedAGSTransition)
    
    def test_create_envelope_en(self):
        """Test creating an EnvelopeEnMethod model."""
        model = create_transition_model('envelope_en')
        self.assertIsInstance(model, EnvelopeEnMethod)
    
    def test_invalid_model_type(self):
        """Test error handling for invalid model type."""
        with self.assertRaises(ValueError):
            create_transition_model('invalid_model_type')

if __name__ == '__main__':
    unittest.main()
