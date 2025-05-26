"""
Thermodynamic relations module for PyMISES.

This module provides functions for calculating thermodynamic
properties for the compressible flow solver in PyMISES.
"""

import numpy as np
from typing import Dict, Tuple, Union, Optional, Any

from pymises.utils.logger import get_logger

logger = get_logger(__name__)

# Constants for air (SI units)
GAMMA = 1.4  # Specific heat ratio for air
R_GAS = 287.05  # Gas constant for air [J/kg/K]
CP = 1004.5  # Specific heat at constant pressure [J/kg/K]
CV = 717.5  # Specific heat at constant volume [J/kg/K]

class Thermodynamics:
    """
    Class for thermodynamic calculations in compressible flow analysis.

    This class provides methods for computing various thermodynamic
    properties for use in the Euler and boundary layer solvers.

    Attributes:
        gamma (float): Specific heat ratio.
        R (float): Gas constant.
        cp (float): Specific heat at constant pressure.
        cv (float): Specific heat at constant volume.
    """

    def __init__(self, gamma: float = GAMMA, R: float = R_GAS,
                 cp: Optional[float] = None, cv: Optional[float] = None,
                 p0: float = 101325.0, T0: float = 288.15):
        """
        Initialize thermodynamics with specified gas constants.

        Args:
            gamma: Specific heat ratio cp/cv.
            R: Gas constant [J/kg/K].
            cp: Specific heat at constant pressure [J/kg/K] (optional, calculated from gamma and R if not provided).
            cv: Specific heat at constant volume [J/kg/K] (optional, calculated from gamma and R if not provided).
            p0: Reference stagnation pressure [Pa].
            T0: Reference stagnation temperature [K].
        """
        self.gamma = gamma
        self.R = R
        self.p0 = p0
        self.T0 = T0

        # If cp and cv are not provided, calculate them from gamma and R
        if cp is None:
            self.cp = self.R * self.gamma / (self.gamma - 1.0)
        else:
            self.cp = cp

        if cv is None:
            self.cv = self.R / (self.gamma - 1.0)
        else:
            self.cv = cv

        # Calculate reference speed of sound
        self.a0 = np.sqrt(self.gamma * self.R * self.T0)

        # Verify consistency
        if not np.isclose(self.gamma, self.cp / self.cv):
            logger.warning(f"Inconsistent thermodynamic parameters: gamma={self.gamma}, cp/cv={self.cp/self.cv}")

    def pressure_from_density_temp(self, density: Union[float, np.ndarray],
                                  temperature: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate pressure from density and temperature using the ideal gas law.

        Args:
            density: Density [kg/m^3].
            temperature: Temperature [K].

        Returns:
            Pressure [Pa].
        """
        return density * self.R * temperature

    def density_from_pressure_temp(self, pressure: Union[float, np.ndarray],
                                  temperature: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate density from pressure and temperature using the ideal gas law.

        Args:
            pressure: Pressure [Pa].
            temperature: Temperature [K].

        Returns:
            Density [kg/m^3].
        """
        return pressure / (self.R * temperature)

    def temperature_from_pressure_density(self, pressure: Union[float, np.ndarray],
                                        density: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate temperature from pressure and density using the ideal gas law.

        Args:
            pressure: Pressure [Pa].
            density: Density [kg/m^3].

        Returns:
            Temperature [K].
        """
        return pressure / (density * self.R)

    def internal_energy(self, temperature: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate specific internal energy from temperature.

        Args:
            temperature: Temperature [K].

        Returns:
            Specific internal energy [J/kg].
        """
        return self.cv * temperature

    def enthalpy(self, temperature: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate specific enthalpy from temperature.

        Args:
            temperature: Temperature [K].

        Returns:
            Specific enthalpy [J/kg].
        """
        return self.cp * temperature

    def temperature_from_enthalpy(self, h: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate temperature from specific enthalpy.

        Args:
            h: Specific enthalpy [J/kg].

        Returns:
            Temperature [K].
        """
        return h / self.cp

    def speed_of_sound(self, temperature: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate speed of sound from temperature.

        Args:
            temperature: Temperature [K].

        Returns:
            Speed of sound [m/s].
        """
        return np.sqrt(self.gamma * self.R * temperature)

    def total_temperature(self, temperature: Union[float, np.ndarray],
                         mach: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate total (stagnation) temperature from static temperature and Mach number.

        Args:
            temperature: Static temperature [K].
            mach: Mach number.

        Returns:
            Total temperature [K].
        """
        return temperature * (1.0 + 0.5 * (self.gamma - 1.0) * mach**2)

    def static_temperature(self, total_temperature: Union[float, np.ndarray],
                          mach: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate static temperature from total temperature and Mach number.

        Args:
            total_temperature: Total temperature [K].
            mach: Mach number.

        Returns:
            Static temperature [K].
        """
        return total_temperature / (1.0 + 0.5 * (self.gamma - 1.0) * mach**2)

    def total_pressure(self, pressure: Union[float, np.ndarray],
                      mach: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate total (stagnation) pressure from static pressure and Mach number.

        Args:
            pressure: Static pressure [Pa].
            mach: Mach number.

        Returns:
            Total pressure [Pa].
        """
        return pressure * (1.0 + 0.5 * (self.gamma - 1.0) * mach**2)**(self.gamma / (self.gamma - 1.0))

    def static_pressure(self, total_pressure: Union[float, np.ndarray],
                       mach: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate static pressure from total pressure and Mach number.

        Args:
            total_pressure: Total pressure [Pa].
            mach: Mach number.

        Returns:
            Static pressure [Pa].
        """
        return total_pressure / (1.0 + 0.5 * (self.gamma - 1.0) * mach**2)**(self.gamma / (self.gamma - 1.0))

    def total_density(self, density: Union[float, np.ndarray],
                     mach: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate total (stagnation) density from static density and Mach number.

        Args:
            density: Static density [kg/m^3].
            mach: Mach number.

        Returns:
            Total density [kg/m^3].
        """
        return density * (1.0 + 0.5 * (self.gamma - 1.0) * mach**2)**(1.0 / (self.gamma - 1.0))

    def static_density(self, total_density: Union[float, np.ndarray],
                      mach: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate static density from total density and Mach number.

        Args:
            total_density: Total density [kg/m^3].
            mach: Mach number.

        Returns:
            Static density [kg/m^3].
        """
        return total_density / (1.0 + 0.5 * (self.gamma - 1.0) * mach**2)**(1.0 / (self.gamma - 1.0))

    def mach_from_pressure_ratio(self, p_ratio: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate Mach number from pressure ratio (p/p0).

        Args:
            p_ratio: Pressure ratio (static/total).

        Returns:
            Mach number.
        """
        # p_ratio = (1.0 + 0.5 * (gamma - 1.0) * M^2)^(-gamma/(gamma-1))
        # Solve for M
        exp = -1.0 * (self.gamma - 1.0) / self.gamma
        term = p_ratio**exp - 1.0
        return np.sqrt(2.0 * term / (self.gamma - 1.0))

    def mach_from_temperature_ratio(self, t_ratio: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate Mach number from temperature ratio (T/T0).

        Args:
            t_ratio: Temperature ratio (static/total).

        Returns:
            Mach number.
        """
        # t_ratio = 1.0 / (1.0 + 0.5 * (gamma - 1.0) * M^2)
        # Solve for M
        return np.sqrt(2.0 * (1.0 / t_ratio - 1.0) / (self.gamma - 1.0))

    def critical_pressure_ratio(self) -> float:
        """
        Calculate the critical pressure ratio for choked flow.

        Returns:
            Critical pressure ratio p*/p0 where p* is the critical pressure.
        """
        return (2.0 / (self.gamma + 1.0))**(self.gamma / (self.gamma - 1.0))

    def critical_temperature_ratio(self) -> float:
        """
        Calculate the critical temperature ratio for choked flow.

        Returns:
            Critical temperature ratio T*/T0 where T* is the critical temperature.
        """
        return 2.0 / (self.gamma + 1.0)

    def critical_density_ratio(self) -> float:
        """
        Calculate the critical density ratio for choked flow.

        Returns:
            Critical density ratio rho*/rho0 where rho* is the critical density.
        """
        return (2.0 / (self.gamma + 1.0))**(1.0 / (self.gamma - 1.0))

    def total_enthalpy(self, velocity: Union[float, np.ndarray],
                      enthalpy: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate total enthalpy from velocity and static enthalpy.

        Args:
            velocity: Velocity magnitude [m/s].
            enthalpy: Static enthalpy [J/kg].

        Returns:
            Total enthalpy [J/kg].
        """
        return enthalpy + 0.5 * velocity**2

    def entropy_change(self, p1: Union[float, np.ndarray], rho1: Union[float, np.ndarray],
                      p2: Union[float, np.ndarray], rho2: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate entropy change between two states.

        Args:
            p1: Initial pressure [Pa].
            rho1: Initial density [kg/m^3].
            p2: Final pressure [Pa].
            rho2: Final density [kg/m^3].

        Returns:
            Entropy change [J/kg/K].
        """
        # ds = cv * ln(p2/p1 * (rho1/rho2)^gamma)
        return self.cv * np.log((p2/p1) * (rho1/rho2)**self.gamma)

    def isentropic_relations(self, mach: Union[float, np.ndarray]) -> Dict[str, Union[float, np.ndarray]]:
        """
        Calculate various isentropic flow ratios for a given Mach number.

        Args:
            mach: Mach number.

        Returns:
            Dictionary containing pressure, temperature, density, and area ratios.
        """
        # Pressure ratio (p/p0)
        p_ratio = (1.0 + 0.5 * (self.gamma - 1.0) * mach**2)**(-self.gamma / (self.gamma - 1.0))

        # Temperature ratio (T/T0)
        t_ratio = 1.0 / (1.0 + 0.5 * (self.gamma - 1.0) * mach**2)

        # Density ratio (rho/rho0)
        rho_ratio = (1.0 + 0.5 * (self.gamma - 1.0) * mach**2)**(-1.0 / (self.gamma - 1.0))

        # Area ratio (A/A*) where A* is the throat area
        area_ratio = (1.0/mach) * ((2.0 + (self.gamma - 1.0) * mach**2) / (self.gamma + 1.0))**((self.gamma + 1.0) / (2.0 * (self.gamma - 1.0)))

        return {
            'pressure_ratio': p_ratio,
            'temperature_ratio': t_ratio,
            'density_ratio': rho_ratio,
            'area_ratio': area_ratio
        }

    def normal_shock_relations(self, mach1: Union[float, np.ndarray]) -> Dict[str, Union[float, np.ndarray]]:
        """
        Calculate flow properties across a normal shock wave.

        Args:
            mach1: Upstream Mach number.

        Returns:
            Dictionary containing downstream Mach number, pressure ratio, temperature ratio,
            density ratio, total pressure ratio, and entropy change.
        """
        # Ensure Mach number is supersonic
        if np.any(mach1 <= 1.0):
            logger.warning("Normal shock relations require supersonic upstream Mach number")
            if isinstance(mach1, np.ndarray):
                mach1 = np.where(mach1 <= 1.0, 1.001, mach1)
            else:
                mach1 = 1.001 if mach1 <= 1.0 else mach1

        # Downstream Mach number
        numerator = 1.0 + 0.5 * (self.gamma - 1.0) * mach1**2
        denominator = self.gamma * mach1**2 - 0.5 * (self.gamma - 1.0)
        mach2 = np.sqrt(numerator / denominator)

        # Pressure ratio (p2/p1)
        p_ratio = 1.0 + 2.0 * self.gamma / (self.gamma + 1.0) * (mach1**2 - 1.0)

        # Temperature ratio (T2/T1)
        t_ratio = (1.0 + 0.5 * (self.gamma - 1.0) * mach1**2) * \
                 (2.0 * self.gamma * mach1**2 - (self.gamma - 1.0)) / \
                 ((self.gamma + 1.0) * mach1**2)

        # Density ratio (rho2/rho1)
        rho_ratio = (self.gamma + 1.0) * mach1**2 / (2.0 + (self.gamma - 1.0) * mach1**2)

        # Total pressure ratio (p02/p01)
        p0_ratio = ((self.gamma + 1.0) * mach1**2 / (2.0 + (self.gamma - 1.0) * mach1**2))**(
            self.gamma / (self.gamma - 1.0)) * \
            (1.0 / ((self.gamma + 1.0) / (2.0 * self.gamma * mach1**2 - (self.gamma - 1.0))))**(1.0 / (self.gamma - 1.0))

        # Entropy change
        ds = self.cv * np.log(p0_ratio**(-1.0))

        return {
            'mach2': mach2,
            'pressure_ratio': p_ratio,
            'temperature_ratio': t_ratio,
            'density_ratio': rho_ratio,
            'total_pressure_ratio': p0_ratio,
            'entropy_change': ds
        }

    def prandtl_meyer_function(self, mach: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate the Prandtl-Meyer function for a given Mach number.

        Args:
            mach: Mach number (must be >= 1).

        Returns:
            Prandtl-Meyer angle [radians].
        """
        # Ensure Mach number is supersonic
        if np.any(mach < 1.0):
            logger.warning("Prandtl-Meyer function requires Mach number >= 1")
            if isinstance(mach, np.ndarray):
                mach = np.where(mach < 1.0, 1.0, mach)
            else:
                mach = 1.0 if mach < 1.0 else mach

        # Calculate Prandtl-Meyer function
        k = np.sqrt((self.gamma + 1.0) / (self.gamma - 1.0))
        nu = k * np.arctan(np.sqrt((mach**2 - 1.0) / (k**2))) - np.arctan(np.sqrt(mach**2 - 1.0))

        return nu

    def mach_from_prandtl_meyer(self, nu: Union[float, np.ndarray],
                               mach_guess: float = 2.0,
                               tolerance: float = 1e-8,
                               max_iter: int = 100) -> Union[float, np.ndarray]:
        """
        Calculate the Mach number from the Prandtl-Meyer angle using an iterative approach.

        Args:
            nu: Prandtl-Meyer angle [radians].
            mach_guess: Initial guess for Mach number.
            tolerance: Convergence tolerance.
            max_iter: Maximum number of iterations.

        Returns:
            Mach number corresponding to the given Prandtl-Meyer angle.
        """
        # If nu is a scalar, convert to array for processing, then convert back at the end
        scalar_input = False
        if np.isscalar(nu):
            nu = np.array([nu])
            scalar_input = True

        # Initialize result array
        mach = np.full_like(nu, mach_guess, dtype=float)

        # Handle special case: nu = 0 => M = 1
        mach[nu == 0] = 1.0

        # For each nu > 0, find the corresponding Mach number
        for i, nu_val in enumerate(nu):
            if nu_val <= 0:
                continue

            # Initial guess
            m = mach_guess
            for _ in range(max_iter):
                pm = self.prandtl_meyer_function(m)
                if np.abs(pm - nu_val) < tolerance:
                    break

                # Use Newton-Raphson method for faster convergence
                # Derivative of Prandtl-Meyer function with respect to Mach number
                k = np.sqrt((self.gamma + 1.0) / (self.gamma - 1.0))
                dpm_dm = np.sqrt(m**2 - 1.0) / m * (
                    k / (1.0 + (m**2 - 1.0) / k**2) - 1.0 / np.sqrt(m**2 - 1.0)
                )

                # Update estimate
                m = m - (pm - nu_val) / dpm_dm

                # Ensure Mach number remains valid
                if m < 1.0:
                    m = 1.001

            mach[i] = m

        return mach[0] if scalar_input else mach

    def oblique_shock_relations(self, mach1: float, shock_angle: float) -> Dict[str, float]:
        """
        Calculate flow properties across an oblique shock wave.

        Args:
            mach1: Upstream Mach number.
            shock_angle: Shock angle with respect to the upstream flow [radians].

        Returns:
            Dictionary containing downstream Mach number, pressure ratio, temperature ratio,
            density ratio, total pressure ratio, and entropy change.
        """
        # Ensure Mach number is supersonic
        if mach1 <= 1.0:
            logger.warning("Oblique shock relations require supersonic upstream Mach number")
            mach1 = 1.001

        # Calculate normal component of Mach number
        mach1n = mach1 * np.sin(shock_angle)

        # Use normal shock relations for the normal component
        normal_shock = self.normal_shock_relations(mach1n)

        # Calculate downstream Mach number
        mach2n = normal_shock['mach2']
        mach2 = mach2n / np.sin(shock_angle - np.arctan(np.tan(shock_angle) / normal_shock['density_ratio']))

        # Add deflection angle and Mach angles to the results
        deflection = np.arctan(2.0 * np.cot(shock_angle) *
                              (mach1**2 * np.sin(shock_angle)**2 - 1.0) /
                              (mach1**2 * (self.gamma + np.cos(2.0 * shock_angle)) + 2.0))

        mach_angle1 = np.arcsin(1.0 / mach1)
        mach_angle2 = np.arcsin(1.0 / mach2)

        result = normal_shock.copy()
        result.update({
            'mach2': mach2,
            'deflection_angle': deflection,
            'mach_angle1': mach_angle1,
            'mach_angle2': mach_angle2
        })

        return result

    def max_deflection_angle(self, mach: float) -> Tuple[float, float]:
        """
        Calculate the maximum deflection angle and corresponding shock angle for a given Mach number.

        Args:
            mach: Upstream Mach number.

        Returns:
            Tuple containing (max_deflection_angle, shock_angle) in radians.
        """
        if mach <= 1.0:
            logger.warning("Max deflection angle calculation requires supersonic Mach number")
            return 0.0, np.pi/2.0

        # Use a numerical approach to find the maximum deflection angle
        # Start with a range of shock angles and narrow down
        angles = np.linspace(np.arcsin(1.0/mach), np.pi/2.0, 100)
        deflections = []

        for angle in angles:
            try:
                result = self.oblique_shock_relations(mach, angle)
                deflections.append(result['deflection_angle'])
            except:
                # If calculation fails, append a negative value
                deflections.append(-1.0)

        # Find the maximum deflection angle
        deflections = np.array(deflections)
        valid_idx = np.where(deflections >= 0.0)[0]
        if len(valid_idx) == 0:
            return 0.0, np.pi/2.0

        max_idx = valid_idx[np.argmax(deflections[valid_idx])]
        max_deflection = deflections[max_idx]
        corresponding_shock_angle = angles[max_idx]

        return max_deflection, corresponding_shock_angle

    def shock_angle_from_deflection(self, mach: float, deflection: float,
                                  weak_shock: bool = True) -> float:
        """
        Calculate the shock angle for a given Mach number and deflection angle.

        Args:
            mach: Upstream Mach number.
            deflection: Flow deflection angle [radians].
            weak_shock: If True, return the weak shock solution, else the strong shock solution.

        Returns:
            Shock angle [radians].
        """
        if mach <= 1.0:
            logger.warning("Shock angle calculation requires supersonic Mach number")
            return np.pi/2.0

        # Ensure deflection is within valid range
        max_deflection, _ = self.max_deflection_angle(mach)
        if deflection > max_deflection:
            logger.warning(f"Deflection angle {deflection*180/np.pi:.2f}° exceeds maximum {max_deflection*180/np.pi:.2f}° for Mach {mach:.2f}")
            return np.pi/2.0

        # Define the function to find roots of
        def f(beta):
            return np.tan(deflection) - 2.0 * np.cot(beta) * (mach**2 * np.sin(beta)**2 - 1.0) / (mach**2 * (self.gamma + np.cos(2.0*beta)) + 2.0)

        # Bisection method to find the shock angle
        # For weak shock: beta between arcsin(1/M) and beta_max
        # For strong shock: beta between beta_max and pi/2

        # First find an estimate of beta_max
        _, beta_max_est = self.max_deflection_angle(mach)

        if weak_shock:
            a = np.arcsin(1.0/mach)
            b = beta_max_est
        else:
            a = beta_max_est
            b = np.pi/2.0

        # Ensure we have opposite signs at the endpoints
        if f(a) * f(b) > 0:
            # If not, adjust the range
            logger.warning("Adjusting shock angle search range")
            if weak_shock:
                a = np.arcsin(1.0/mach)
                b = np.pi/2.0
            else:
                logger.warning("Could not find strong shock solution")
                return np.pi/2.0

        # Bisection method
        tolerance = 1e-10
        max_iter = 100
        for _ in range(max_iter):
            c = (a + b) / 2.0
            fc = f(c)

            if abs(fc) < tolerance or (b - a) < tolerance:
                return c

            if f(a) * fc < 0:
                b = c
            else:
                a = c

        logger.warning("Shock angle calculation did not converge")
        return (a + b) / 2.0

    def artificial_dissipation_coef(self, mach: Union[float, np.ndarray],
                                   threshold_mach: float = 0.95,
                                   dissipation_coeff: float = 0.5) -> Union[float, np.ndarray]:
        """
        Calculate artificial dissipation coefficient for Euler solver.

        Args:
            mach: Local Mach number.
            threshold_mach: Threshold Mach number for dissipation activation.
            dissipation_coeff: Dissipation coefficient.

        Returns:
            Artificial dissipation coefficient.
        """
        # Only apply dissipation in supersonic regions
        if np.isscalar(mach):
            if mach < threshold_mach:
                return 0.0
            else:
                return dissipation_coeff * (1.0 - (threshold_mach/mach)**2) / (2.0 * mach**2)
        else:
            result = np.zeros_like(mach)
            supersonic = mach >= threshold_mach
            result[supersonic] = dissipation_coeff * (1.0 - (threshold_mach/mach[supersonic])**2) / (2.0 * mach[supersonic]**2)
            return result
