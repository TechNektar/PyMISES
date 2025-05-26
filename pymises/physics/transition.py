"""
Transition prediction models for PyMISES.

This module provides implementations of laminar-turbulent transition
prediction methods, including the modified Abu-Ghannam/Shaw criteria.
"""

import numpy as np
from typing import Dict, Tuple, Union, Optional, Any

from pymises.utils.logger import get_logger

logger = get_logger(__name__)

class TransitionModel:
    """
    Base class for transition prediction models.

    This class provides the interface for transition prediction models
    used in the boundary layer solver to determine the laminar-turbulent
    transition location.

    Attributes:
        config: Configuration dictionary with model parameters.
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the transition model.

        Args:
            config: Configuration dictionary with model parameters (optional).
        """
        self.config = config or {}

    def critical_reynolds_number(self, H: float, turbulence_level: float) -> float:
        """
        Calculate critical Reynolds number for transition onset.

        Args:
            H: Shape parameter.
            turbulence_level: Freestream turbulence level (0-1).

        Returns:
            Critical Reynolds number based on momentum thickness.
        """
        raise NotImplementedError("Subclasses must implement critical_reynolds_number")

    def amplification_rate(self, H: float, Re_theta: float, Re_theta_crit: float) -> float:
        """
        Calculate amplification rate for transition prediction.

        Args:
            H: Shape parameter.
            Re_theta: Reynolds number based on momentum thickness.
            Re_theta_crit: Critical Reynolds number for transition onset.

        Returns:
            Amplification rate.
        """
        raise NotImplementedError("Subclasses must implement amplification_rate")

    def update_amplification_factor(self, n: float, amplification_rate: float,
                                   dx: float, edge_velocity_ratio: float) -> float:
        """
        Update amplification factor based on local conditions.

        Args:
            n: Current amplification factor.
            amplification_rate: Local amplification rate.
            dx: Streamwise distance step.
            edge_velocity_ratio: Ratio of new to old edge velocity.

        Returns:
            Updated amplification factor.
        """
        raise NotImplementedError("Subclasses must implement update_amplification_factor")

    def is_transition(self, n: float, n_crit: float) -> bool:
        """
        Determine if transition occurs based on amplification factor.

        Args:
            n: Current amplification factor.
            n_crit: Critical amplification factor for transition.

        Returns:
            True if transition occurs, False otherwise.
        """
        return n >= n_crit


class ModifiedAGSTransition(TransitionModel):
    """
    Modified Abu-Ghannam/Shaw transition model.

    This model implements Drela's modified version of the Abu-Ghannam/Shaw
    transition criterion, which addresses ill-posedness issues in the original
    formulation by using shape parameter H instead of pressure gradient parameter lambda.

    Attributes:
        A: Amplification rate constant for bypass transition.
        B: Ramp width parameter for bypass transition.
        turbulence_level: Freestream turbulence level.
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the modified AGS transition model.

        Args:
            config: Configuration dictionary with model parameters:
                - A: Amplification rate constant (default: 0.10)
                - B: Ramp width parameter (default: 0.30)
                - turbulence_level: Freestream turbulence level (default: 0.01)
        """
        super().__init__(config)

        self.A = self.config.get('amplification_constant', 0.10)
        self.B = self.config.get('ramp_width', 0.30)
        self.turbulence_level = self.config.get('turbulence_level', 0.01)

        logger.info(f"Initialized modified AGS transition model with A={self.A}, B={self.B}, "
                   f"turbulence_level={self.turbulence_level}")

    def critical_reynolds_number(self, H: float, turbulence_level: Optional[float] = None) -> float:
        """
        Calculate critical Reynolds number for transition onset using modified AGS criterion.

        Args:
            H: Shape parameter.
            turbulence_level: Freestream turbulence level (0-1, optional).
                              If not provided, uses the value from initialization.

        Returns:
            Critical Reynolds number based on momentum thickness.
        """
        # Use provided turbulence level or default
        tu = turbulence_level if turbulence_level is not None else self.turbulence_level

        # Limit turbulence level to avoid numerical issues
        tu = max(min(tu, 1.0), 1e-6)

        # Modified turbulence level term from Drela's paper
        tau_prime = 2.7 * np.tanh(tu / 2.7)

        # Critical amplification factor based on turbulence level
        n_crit = -8.43 - 2.4 * np.log(tau_prime / 100.0)

        # H-based parameterization that works even in separation regions
        tanh_term = np.tanh(10.0 / (H - 1.0) - 5.5)

        # Critical Reynolds number correlation
        Re_theta_crit = 155.0 + 89.0 * (0.25 * tanh_term + 1.0) * (n_crit**1.25)

        return max(Re_theta_crit, 155.0)  # Ensure minimum value

    def ts_wave_growth_rate(self, H: float, Re_theta: float) -> float:
        """
        Calculate Tollmien-Schlichting wave growth rate based on Orr-Sommerfeld results.

        Args:
            H: Shape parameter.
            Re_theta: Reynolds number based on momentum thickness.

        Returns:
            TS wave growth rate.
        """
        # Simple correlation based on Orr-Sommerfeld results
        # Higher growth rates for adverse pressure gradients (higher H)
        H_term = max(0.0, H - 2.3)
        Re_term = min(1.0, Re_theta / 300.0)

        # Base growth rate
        return 0.01 * Re_term * (1.0 + 5.0 * H_term**2)

    def bypass_growth_rate(self, r: float) -> float:
        """
        Calculate bypass transition growth rate using cubic ramp function.

        Args:
            r: Normalized Reynolds number parameter.

        Returns:
            Bypass growth rate.
        """
        if r < 0:
            return 0.0
        elif r < 1:
            # Cubic ramp function for smooth transition
            return self.A * (3.0 * r**2 - 2.0 * r**3)
        else:
            return self.A

    def amplification_rate(self, H: float, Re_theta: float, Re_theta_crit: float) -> float:
        """
        Calculate combined amplification rate for transition prediction.

        Args:
            H: Shape parameter.
            Re_theta: Reynolds number based on momentum thickness.
            Re_theta_crit: Critical Reynolds number for transition onset.

        Returns:
            Combined TS-wave and bypass amplification rate.
        """
        # TS wave component from Orr-Sommerfeld
        f_rate = self.ts_wave_growth_rate(H, Re_theta)

        # Compute normalized Reynolds number parameter for bypass transition
        r = (1.0 / self.B) * ((Re_theta / Re_theta_crit) - 1.0) + 0.5

        # Bypass transition component using cubic ramp function
        g_rate = self.bypass_growth_rate(r)

        # Combined growth rate
        return f_rate + g_rate

    def update_amplification_factor(self, n: float, amplification_rate: float,
                                   dx: float, edge_velocity_ratio: float = 1.0) -> float:
        """
        Update amplification factor with edge velocity correction.

        Args:
            n: Current amplification factor.
            amplification_rate: Local amplification rate.
            dx: Streamwise distance step.
            edge_velocity_ratio: Ratio of new to old edge velocity.

        Returns:
            Updated amplification factor.
        """
        # Edge velocity correction term (accounts for acceleration/deceleration effects)
        dln_ue = np.log(edge_velocity_ratio)

        # Modified amplification equation (eq. 31 in Drela's paper)
        dn = amplification_rate * dx - dln_ue

        # Update amplification factor
        n_new = n + dn

        # Ensure n doesn't go negative
        return max(0.0, n_new)

    def initial_condition(self, x: float, H: float) -> Tuple[float, float]:
        """
        Set initial condition for transition prediction at start of boundary layer.

        Args:
            x: Initial x-coordinate.
            H: Initial shape parameter.

        Returns:
            Tuple of (x_initial, n_initial) for transition prediction.
        """
        # Initially zero amplification factor at the start of the boundary layer
        return x, 0.0


class EnvelopeEnMethod(TransitionModel):
    """
    Envelope e^n transition prediction method.

    This class implements the envelope method for e^n transition prediction,
    which uses a simplified approach based on maximum growth rates across
    all frequencies.

    Attributes:
        turbulence_level: Freestream turbulence level.
        n_crit: Critical amplification factor.
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the envelope e^n transition model.

        Args:
            config: Configuration dictionary with model parameters:
                - turbulence_level: Freestream turbulence level (default: 0.01)
        """
        super().__init__(config)

        self.turbulence_level = self.config.get('turbulence_level', 0.01)

        # Calculate n_crit based on turbulence level
        tu = max(min(self.turbulence_level, 1.0), 1e-6)
        self.n_crit = -8.43 - 2.4 * np.log(tu / 100.0)

        logger.info(f"Initialized envelope e^n transition model with turbulence_level={self.turbulence_level}, "
                   f"n_crit={self.n_crit}")

    def critical_reynolds_number(self, H: float, turbulence_level: Optional[float] = None) -> float:
        """
        Calculate critical Reynolds number for neutral stability.

        Args:
            H: Shape parameter.
            turbulence_level: Freestream turbulence level (not used in this method).

        Returns:
            Critical Reynolds number for neutral stability.
        """
        # Reynolds number where instability first appears (neutral stability point)
        # Correlation based on Orr-Sommerfeld results for Falkner-Skan profiles

        # Convert H to Pohlhausen parameter lambda
        if H < 2.5:
            pohlhausen = 0.22 - 1.57 * (2.5 - H)  # Favorable gradient (lambda > 0)
        else:
            pohlhausen = 0.22 - 0.8 * (H - 2.5)   # Adverse gradient (lambda < 0)

        # Limit lambda to reasonable range
        pohlhausen = max(-0.18, min(0.1, pohlhausen))

        # Critical Reynolds number correlation
        if pohlhausen >= 0:
            return 125.0 + 2000.0 * pohlhausen
        else:
            return 125.0 + 1150.0 * pohlhausen

    def amplification_rate(self, H: float, Re_theta: float, Re_theta_crit: float) -> float:
        """
        Calculate amplification rate using envelope method.

        Args:
            H: Shape parameter.
            Re_theta: Reynolds number based on momentum thickness.
            Re_theta_crit: Critical Reynolds number for neutral stability.

        Returns:
            Amplification rate.
        """
        # No amplification below critical Reynolds number
        if Re_theta < Re_theta_crit:
            return 0.0

        # Envelope method approximation for maximum growth rate
        # Based on Drela's approximation of envelope method results

        # Normalized distance from critical point
        normalized_distance = (Re_theta - Re_theta_crit) / Re_theta_crit

        # Base growth rate
        growth_rate = 0.01 * normalized_distance

        # Adjust for pressure gradient effect through H
        if H > 2.5:  # Adverse gradient
            adjustment = 1.0 + 3.0 * (H - 2.5)**2
            growth_rate *= adjustment

        return growth_rate

    def update_amplification_factor(self, n: float, amplification_rate: float,
                                   dx: float, edge_velocity_ratio: float = 1.0) -> float:
        """
        Update amplification factor with edge velocity correction.

        Args:
            n: Current amplification factor.
            amplification_rate: Local amplification rate.
            dx: Streamwise distance step.
            edge_velocity_ratio: Ratio of new to old edge velocity.

        Returns:
            Updated amplification factor.
        """
        # Edge velocity ratio effect
        dln_ue = np.log(edge_velocity_ratio)

        # Update amplification factor (spatial amplification rate)
        dn = amplification_rate * dx - dln_ue

        # Ensure n doesn't go negative
        return max(0.0, n + dn)

    def initial_condition(self, x: float, H: float) -> Tuple[float, float]:
        """
        Set initial condition for transition prediction at neutral stability point.

        Args:
            x: Initial x-coordinate.
            H: Initial shape parameter.

        Returns:
            Tuple of (x_initial, n_initial) for transition prediction.
        """
        # Start with zero amplification factor
        return x, 0.0


class TransitionPredictor:
    """
    Predicts boundary layer transition using various models.

    Attributes:
        model: Transition prediction model instance.
        config: Configuration dictionary with predictor parameters.
        n_factor: Array of amplification factors.
        transition_occurred: Flag indicating if transition has occurred.
    """

    def __init__(self, model=None, model_type: str = 'modified_ags', config: Optional[Dict[str, Any]] = None):
        """
        Initialize the transition predictor.

        Args:
            model: Direct instance of a transition model (takes precedence if provided)
            model_type: Type of transition model to use if model not provided.
            config: Configuration dictionary.
        """
        self.config = config or {}

        # Create transition model
        if model is not None:
            # Use the provided model directly
            self.model = model
        elif model_type == 'modified_ags':
            self.model = ModifiedAGSTransition(self.config)
        else:
            raise ValueError(f"Unknown transition model type: {model_type}")

        # Initialize arrays
        self.n_points = 0
        self.x = None
        self.edge_velocity = None
        self.theta = None
        self.reynolds_number = None
        self.n_factor = None
        self.amplification_rate = None
        self.transition_occurred = False
        self.transition_index = None

        model_name = self.model.__class__.__name__
        logger.info(f"Initialized transition predictor with {model_name} model")

    def initialize(self, x: np.ndarray, edge_velocity: np.ndarray) -> None:
        """Initialize transition predictor with grid and flow properties.

        Args:
            x: Streamwise coordinate array.
            edge_velocity: Edge velocity array.
        """
        self.x = x
        self.edge_velocity = edge_velocity
        self.n_points = len(x)

        # Initialize arrays for transition prediction
        self.n_factor = np.zeros(self.n_points)
        self.amplification_rate = np.zeros(self.n_points)
        self.transition_occurred = False
        self.transition_index = None

        # Initialize model parameters
        config = {
            'amplification_constant': self.config.get('A', 0.1),
            'ramp_width': self.config.get('B', 0.3),
            'turbulence_level': self.config.get('turbulence_level', 0.01)
        }
        self.model = create_transition_model('modified_ags', config)

        logger.info(f"Initialized transition predictor with modified_ags model")

    def update(self, i: int, H: float, Ue: float) -> None:
        """
        Update the transition prediction.

        Args:
            i: Current point index.
            H: Current shape parameter.
            Ue: Current edge velocity.
        """
        try:
            # Ensure H and Ue are in reasonable ranges
            H = max(1.05, min(4.0, H))  # Constrain H to reasonable values
            Ue = max(1e-8, Ue)  # Ensure Ue is positive

            # Update n-factor
            if i > 0 and i < len(self.n_factor):
                # Compute amplification rate
                if H > 2.2:  # Only amplify if H is above critical value
                    try:
                        # Use a simple approximation for dH/dx based on H
                        dH_dx = 0.1 * (H - 2.2)  # Positive in adverse pressure gradient

                        # Compute Reynolds number based on momentum thickness
                        if hasattr(self, 'theta') and i < len(self.theta):
                            theta_i = self.theta[i]
                        else:
                            # Use a default value if theta is not available
                            theta_i = 1e-4

                        Re_theta = self.reynolds_number * Ue * theta_i
                        Re_theta = max(1.0, Re_theta)  # Ensure positive Reynolds number

                        Re_theta_crit = self.model.critical_reynolds_number(H)
                        Re_theta_crit = max(1.0, Re_theta_crit)  # Ensure positive critical Reynolds number

                        # Compute amplification rate
                        self.amplification_rate[i] = self.model.amplification_rate(H, Re_theta, Re_theta_crit)
                    except Exception as e:
                        logger.warning(f"Error computing amplification rate: {str(e)}")
                        self.amplification_rate[i] = 0.0
                else:
                    self.amplification_rate[i] = 0.0

                # Update n-factor with error handling
                try:
                    dx = self.x[i] - self.x[i-1]
                    dx = max(1e-8, dx)  # Ensure positive dx
                    self.n_factor[i] = self.n_factor[i-1] + self.amplification_rate[i] * dx
                except Exception as e:
                    logger.warning(f"Error updating n-factor: {str(e)}")
                    self.n_factor[i] = self.n_factor[i-1]  # Use previous value if update fails
            else:
                # Initialize first point or handle out-of-bounds index
                if i < len(self.amplification_rate):
                    self.amplification_rate[i] = 0.0
                if i < len(self.n_factor):
                    self.n_factor[i] = 0.0
        except Exception as e:
            logger.warning(f"Transition prediction update failed: {str(e)}")

    def check_transition(self, x: float, Re_theta: float, H: float, dUe_dx: float) -> bool:
        """Check if transition occurs at the current location.

        Args:
            x: Streamwise coordinate.
            Re_theta: Reynolds number based on momentum thickness.
            H: Shape parameter.
            dUe_dx: Edge velocity gradient.

        Returns:
            True if transition occurs, False otherwise.
        """
        try:
            # Ensure inputs are in reasonable ranges
            x = max(0.0, x)  # Ensure x is non-negative
            Re_theta = max(1.0, Re_theta)  # Ensure Re_theta is positive
            H = max(1.05, min(4.0, H))  # Constrain H to reasonable values

            # Skip transition check near leading edge
            if x < 0.05:
                return False

            try:
                # Compute critical Reynolds number
                Re_theta_crit = self.model.critical_reynolds_number(H, dUe_dx)
                Re_theta_crit = max(1.0, Re_theta_crit)  # Ensure positive critical Reynolds number

                # Check for bypass transition
                if Re_theta > Re_theta_crit:
                    logger.info(f"Bypass transition detected at x = {x:.3f}")
                    self.transition_occurred = True
                    return True
            except Exception as e:
                logger.warning(f"Error in bypass transition check: {str(e)}")

            # Check for separation-induced transition
            if H > 3.5 and Re_theta > 200:
                logger.info(f"Separation-induced transition detected at x = {x:.3f}")
                self.transition_occurred = True
                return True
        except Exception as e:
            logger.warning(f"Transition check failed: {str(e)}")

        return False

    def is_transition(self) -> bool:
        """
        Check if transition occurs at the current point.

        Returns:
            True if transition occurs, False otherwise.
        """
        try:
            # Get current n-factor and calculate critical n-factor
            if len(self.n_factor) > 0:
                n = self.n_factor[-1]
            else:
                return False  # No n-factor data available

            # Ensure turbulence level is in a reasonable range
            turbulence_level = max(1e-6, min(0.2, self.model.turbulence_level))  # Between 0.0001% and 20%

            # Calculate critical n-factor with safety checks
            try:
                n_crit = -8.43 - 2.4 * np.log(turbulence_level / 100.0)
                n_crit = max(0.1, min(20.0, n_crit))  # Reasonable bounds for n_crit
            except Exception as e:
                logger.warning(f"Error calculating critical n-factor: {str(e)}")
                n_crit = 9.0  # Default value for 0.1% turbulence

            # Check if transition occurs
            try:
                return self.model.is_transition(n, n_crit)
            except Exception as e:
                logger.warning(f"Error in transition check: {str(e)}")
                # Fallback to simple comparison
                return n > n_crit
        except Exception as e:
            logger.warning(f"Transition detection failed: {str(e)}")
            return False

    def get_state(self) -> Dict[str, Any]:
        """
        Get current state of transition prediction.

        Returns:
            Dictionary with current state:
            - n_factor: Array of amplification factors
            - amplification_rate: Array of amplification rates
            - transition_occurred: Flag indicating if transition has occurred
            - transition_index: Index where transition occurred (if any)
        """
        return {
            'n_factor': self.n_factor,
            'amplification_rate': self.amplification_rate,
            'transition_occurred': self.transition_occurred,
            'transition_index': self.transition_index
        }

    def predict_transition(self, x: np.ndarray, H: np.ndarray, Re_theta: np.ndarray,
                         edge_velocity: np.ndarray) -> Dict[str, Any]:
        """
        Run transition prediction on a given boundary layer solution.

        Args:
            x: Streamwise coordinate array.
            H: Shape parameter array.
            Re_theta: Reynolds number based on momentum thickness array.
            edge_velocity: Edge velocity array.

        Returns:
            Dictionary with transition prediction results:
            - n_factor: Array of amplification factors
            - transition_occurred: Flag indicating if transition has occurred
            - transition_index: Index where transition occurred (if any)
            - transition_x: x-coordinate where transition occurred (if any)
        """
        # Initialize arrays
        n_points = len(x)
        n_factor = np.zeros(n_points)
        transition_occurred = False
        transition_index = None
        transition_x = None

        # Check array lengths
        if len(H) != n_points or len(Re_theta) != n_points or len(edge_velocity) != n_points:
            raise ValueError("All input arrays must have the same length")

        # Calculate critical n-factor
        n_crit = -8.43 - 2.4 * np.log(self.model.turbulence_level / 100.0)

        # Loop over all points
        for i in range(1, n_points):
            # Skip if transition has already occurred
            if transition_occurred:
                n_factor[i] = n_factor[i-1]
                continue

            # Calculate critical Reynolds number
            Re_theta_crit = self.model.critical_reynolds_number(H[i])

            # Calculate edge velocity ratio
            edge_vel_ratio = edge_velocity[i] / max(edge_velocity[i-1], 1e-6)

            # Calculate amplification rate
            rate = self.model.amplification_rate(H[i], Re_theta[i], Re_theta_crit)

            # Update n-factor
            n_factor[i] = self.model.update_amplification_factor(
                n_factor[i-1], rate, x[i] - x[i-1], edge_vel_ratio
            )

            # Check for transition
            if n_factor[i] >= n_crit or Re_theta[i] > 2*Re_theta_crit or H[i] > 3.5:
                transition_occurred = True
                transition_index = i
                transition_x = x[i]
                logger.info(f"Transition predicted at x = {x[i]:.4f}, Reynolds = {Re_theta[i]:.1f}")

        # Return results
        return {
            'n_factor': n_factor,
            'transition_occurred': transition_occurred,
            'transition_index': transition_index,
            'transition_x': transition_x
        }


def create_transition_model(model_type: str, config: Optional[Dict[str, Any]] = None) -> TransitionModel:
    """
    Factory function to create a transition model.

    Args:
        model_type: Type of transition model ('modified_ags' or 'envelope_en').
        config: Configuration dictionary with model parameters.

    Returns:
        Instantiated transition model.

    Raises:
        ValueError: If the specified model type is not supported.
    """
    if model_type.lower() == 'modified_ags':
        return ModifiedAGSTransition(config)
    elif model_type.lower() == 'envelope_en':
        return EnvelopeEnMethod(config)
    else:
        raise ValueError(f"Unsupported transition model type: {model_type}. "
                        f"Supported types are: 'modified_ags', 'envelope_en'")
