# PyMISES Physics Module

This directory contains implementations of physical models used in the PyMISES solver. These modules provide the fundamental physical relationships and closure models that govern fluid flow behavior.

## Overview

The physics module provides:

1. **Thermodynamic Models**: Equations of state and thermodynamic relations
2. **Transition Prediction**: Models for laminar-turbulent transition
3. **Turbulence Models**: Turbulent closure relations for boundary layer equations
4. **Dissipation Models**: Physical and artificial dissipation models

## Modules

### `thermo.py`

This module implements thermodynamic relations and equations of state for compressible flow.

#### Key Functions and Classes

- **`ThermodynamicModel`**: Base class for thermodynamic models.
  - Provides interface for computing thermodynamic properties
  - Supports various equations of state
  - Handles unit conversions and normalization

- **`IdealGas`**: Implements ideal gas thermodynamic relations.
  - Perfect gas law (p = ρRT)
  - Specific heat models (constant, temperature-dependent)
  - Enthalpy and entropy calculations
  - Isentropic flow relations

- **`CompressibleFlowFunctions`**: Implements compressible flow relationships.
  - Isentropic relations between Mach number and flow properties
  - Normal shock relations
  - Oblique shock and expansion relations
  - Fanno and Rayleigh flow relations

#### Implementation Details

- Ideal gas law implemented as:
  
  $p = \rho R T$
  
  where $p$ is pressure, $\rho$ is density, $R$ is the gas constant, and $T$ is temperature.

- Isentropic relations for compressible flow:
  
  $\frac{T_0}{T} = 1 + \frac{\gamma-1}{2}M^2$
  
  $\frac{p_0}{p} = \left(1 + \frac{\gamma-1}{2}M^2\right)^{\frac{\gamma}{\gamma-1}}$
  
  $\frac{\rho_0}{\rho} = \left(1 + \frac{\gamma-1}{2}M^2\right)^{\frac{1}{\gamma-1}}$
  
  where $M$ is the Mach number, $\gamma$ is the ratio of specific heats, and the subscript 0 denotes stagnation conditions.

- Speed of sound calculation:
  
  $a = \sqrt{\gamma R T} = \sqrt{\gamma p / \rho}$

- Entropy calculation:
  
  $s - s_{ref} = c_p \ln\left(\frac{T}{T_{ref}}\right) - R \ln\left(\frac{p}{p_{ref}}\right)$

- Tools for converting between different variable sets:
  - Primitive variables (ρ, u, v, p)
  - Conserved variables (ρ, ρu, ρv, E)
  - Characteristic variables

### `transition.py`

This module implements models for predicting laminar-turbulent transition in boundary layers.

#### Key Classes

- **`TransitionModel`**: Base class for transition prediction models.
  - Defines common interface for all transition models
  - Tracks amplification factor and other transition metrics
  - Provides methods to detect transition onset

- **`ModifiedAGSTransition`**: Implements the modified Abu-Ghannam/Shaw model.
  - H-based formulation for improved numerical stability
  - Combined formulation for TS-wave and bypass mechanisms
  - Edge velocity correction for non-equilibrium flows
  - Based on Drela's 1998 paper on transition modeling

- **`EnvelopeEnMethod`**: Implements the envelope $e^n$ method.
  - Tracks amplification of the most unstable frequency
  - Uses database of growth rates from stability analysis
  - Identifies transition when n-factor exceeds critical value
  - Incorporates curvature and crossflow effects

- **`TransitionPredictor`**: Manages transition prediction process.
  - Orchestrates data flow between flow solver and transition model
  - Handles transition model initialization and state tracking
  - Processes transition events and model switching
  - Provides boundary layer parameters to the transition model

#### Implementation Details

- The modified Abu-Ghannam/Shaw model calculates critical Reynolds number as:
  
  $Re_{\theta,crit}(H, Tu) = 155 + 89.0 \left[0.25 \tanh\left(\frac{10}{H-1} - 5.5\right) + 1\right] (n_{crit})^{1.25}$
  
  where $H$ is the shape parameter, $Tu$ is the turbulence level, and $n_{crit}$ is related to $Tu$ by:
  
  $n_{crit} = -8.43 - 2.4 \ln\left(\frac{Tu'}{100}\right)$
  
  with $Tu' = 2.7 \tanh\left(\frac{Tu}{2.7}\right)$

- Combined amplification rate formulation:
  
  $\frac{dn}{dx} = f(H, Re_\theta) + g(H, Re_\theta)$
  
  where $f$ is the TS-wave component and $g$ is the bypass component:
  
  $g(H, Re_\theta) = \begin{cases}
  0, & r < 0 \\
  A(3r^2 - 2r^3), & 0 \leq r \leq 1 \\
  A, & r > 1
  \end{cases}$
  
  with $r = \frac{1}{B}\left(\frac{Re_\theta}{Re_{\theta,crit}} - 1\right) + \frac{1}{2}$

- Edge velocity correction for varying pressure gradients:
  
  $\frac{dn}{dx} + \frac{d \ln u_e}{dx} = \frac{1}{\theta}\left[f(H, Re_\theta) + g(H, Re_\theta)\right]$

- Transition is detected when the amplification factor $n$ exceeds a critical value (typically 9-11 for low turbulence environments, lower for higher turbulence).

### `laminar.py`

This module implements closure models for laminar boundary layer flows.

#### Key Functions

- **`laminar_skin_friction`**: Computes skin friction for laminar flow.
  - Based on the Falkner-Skan family of similarity solutions
  - Accounts for pressure gradient effects through shape parameter
  - Includes compressibility corrections

- **`laminar_dissipation`**: Computes dissipation coefficient for laminar flow.
  - Energy dissipation due to viscous effects
  - Relates to boundary layer profile shape
  - Used in the shape parameter equation

- **`laminar_heat_transfer`**: Computes heat transfer for laminar flow.
  - Conduction through boundary layer
  - Reynolds analogy relations
  - Temperature recovery effects

#### Implementation Details

- Laminar skin friction coefficient based on shape parameter:
  
  $C_f Re_\theta = -0.067 + 0.01977 \frac{(7.4-H_k)^2}{H_k-1}$
  
  where $H_k$ is the kinematic shape parameter.

- Dissipation coefficient relation:
  
  $\frac{2C_D}{H^*} = \begin{cases}
  0.207 + 0.00205(4-H_k)^{5.5}, & H_k < 4 \\
  0.207 - 0.003(H_k-4)^2, & H_k \geq 4
  \end{cases}$
  
  where $H^*$ is the kinetic energy shape parameter.

- Kinematic shape parameter related to conventional shape parameter:
  
  $H_k = H - \frac{M_e^2 - 1}{2 + \sqrt{H-1}} \left(1 - \frac{1}{H}\right)$

- Shape parameter related to non-dimensional pressure gradient:
  
  $H = H(\lambda)$ where $\lambda = \frac{\theta^2}{\nu} \frac{du_e}{dx}$ 
  
  is the Pohlhausen parameter.

### `turbulent.py`

This module implements closure models for turbulent boundary layer flows.

#### Key Functions

- **`turbulent_skin_friction`**: Computes skin friction for turbulent flow.
  - Based on empirical correlations for turbulent profiles
  - Accounts for Reynolds number and shape parameter effects
  - Includes compressibility corrections

- **`turbulent_dissipation`**: Computes dissipation coefficient for turbulent flow.
  - Energy dissipation in turbulent shear layer
  - Lag model for non-equilibrium effects
  - Separated flow extensions

- **`turbulent_heat_transfer`**: Computes heat transfer for turbulent flow.
  - Reynolds analogy with modifications
  - Accounts for turbulent Prandtl number effects
  - Temperature recovery factor calculations

- **`turbulent_closure`**: Comprehensive closure model for turbulent flow.
  - Relates all required parameters for integral formulation
  - Consistent set of relationships based on advanced turbulence modeling
  - Validated against experimental data

#### Implementation Details

- Turbulent skin friction based on shape parameter and Reynolds number:
  
  $C_f = \frac{0.3}{(\log_{10} Re_\theta)^{1.74 + 0.31H}} \cdot F_c$
  
  where $F_c$ is a compressibility correction:
  
  $F_c = \frac{1}{(1 + 0.2M_e^2)^{0.5}}$

- Turbulent energy shape parameter correlation:
  
  $H^* = \frac{1.505 + \frac{400}{Re_\theta}}{H_k} \left(1 + \frac{(H_\tau - H)^2}{H_k} \left(0.04 + 0.0007\ln Re_\theta\right) \right)$
  
  where $H_\tau$ is the limiting shape factor for turbulent separation.

- Lag equation for Reynolds stress coefficient:
  
  $\delta \frac{dC_\tau}{dx} = K_C \left(\sqrt{C_{\tau,EQ}} - \sqrt{C_\tau}\right)$
  
  where $C_{\tau,EQ}$ is the equilibrium value and $K_C = 5.6$ is the lag constant.

- Dissipation coefficient formulation including lag effects:
  
  $C_D = C_f \cdot U_s^2 + C_\tau (1 - U_s)$
  
  where $U_s$ is a shear correlation parameter.

### `dissipation.py`

This module implements models for physical and artificial dissipation in the flow.

#### Key Functions and Classes

- **`PhysicalDissipation`**: Computes physical energy dissipation.
  - Viscous effects in boundary layer
  - Relates to skin friction and velocity profile
  - Used in energy equation

- **`ArtificialDissipation`**: Implements artificial dissipation models.
  - Adds numerical dissipation for shock capturing
  - Pressure-switched activation for transonic regions
  - Preserves physical accuracy in smooth regions

- **`ShockCapturing`**: Specialized shock capturing methods.
  - Detects shock locations
  - Adds appropriate dissipation near shocks
  - Maintains solution accuracy away from shocks

#### Implementation Details

- MISES-specific artificial dissipation model:
  
  $\mu = \begin{cases}
  0, & M < M_c \\
  C_d(1-(M_c/M)^2)/(2M^2), & M > M_c
  \end{cases}$
  
  where $M_c$ is a threshold Mach number (typically 0.9-0.95) and $C_d$ is a dissipation coefficient.

- Physical dissipation in the boundary layer:
  
  $C_D = \frac{1}{\theta} \int_0^\delta \frac{\mu}{\rho} \left(\frac{\partial u}{\partial y}\right)^2 dy$
  
  which is related to the skin friction through closure relations specific to laminar or turbulent flow.

- Matrix dissipation formulation for improved shock capturing:
  
  $D(U)_i = \alpha_i (U_{i+1} - U_i) - \alpha_{i-1} (U_i - U_{i-1})$
  
  where $\alpha_i$ is a blend of first and third-order differences:
  
  $\alpha_i = \alpha^{(1)}_i - \alpha^{(3)}_i$

- Shock detection through pressure gradient:
  
  $\sigma_i = \left|\frac{p_{i+1} - 2p_i + p_{i-1}}{p_{i+1} + 2p_i + p_{i-1}}\right|$
  
  with the first-order dissipation coefficient:
  
  $\alpha^{(1)}_i = \kappa^{(1)} \max(\sigma_{i-1}, \sigma_i, \sigma_{i+1})$

## Dependencies

- NumPy: For array operations
- SciPy: For interpolation and integration

## Example Usage

### Thermodynamic Calculations

```python
from pymises.physics.thermo import IdealGas

# Create ideal gas model with γ = 1.4 (air)
gas = IdealGas(gamma=1.4, gas_constant=287.0)

# Calculate properties at specified Mach number
mach = 0.8
p_ratio = gas.pressure_ratio(mach)  # p/p0
t_ratio = gas.temperature_ratio(mach)  # T/T0
d_ratio = gas.density_ratio(mach)  # ρ/ρ0

# Calculate properties for a specific state
p = 101325.0  # Pa
t = 288.15    # K
rho = gas.density(p, t)
a = gas.speed_of_sound(t)
h = gas.enthalpy(t)
s = gas.entropy(p, t, p_ref=101325.0, t_ref=288.15)

# Calculate isentropic stagnation properties
p0 = p / p_ratio
t0 = t / t_ratio

print(f"Mach: {mach}, p/p0: {p_ratio}, T/T0: {t_ratio}")
print(f"Density: {rho} kg/m³, Sound speed: {a} m/s")
print(f"Stagnation pressure: {p0} Pa, Stagnation temperature: {t0} K")
```

### Transition Prediction

```python
from pymises.physics.transition import ModifiedAGSTransition, TransitionPredictor

# Create a transition model
transition_model = ModifiedAGSTransition({
    'amplification_constant': 0.10,
    'ramp_width': 0.30,
    'turbulence_level': 0.01
})

# Create a transition predictor
predictor = TransitionPredictor(
    model_type='modified_ags',
    config={'turbulence_level': 0.01}
)

# Example boundary layer data along airfoil
x = np.linspace(0.001, 1.0, 101)  # Surface distance from leading edge
h = 2.3 * np.ones_like(x)         # Shape parameter (initially laminar)
re_theta = np.linspace(100, 1000, 101)  # Reynolds number based on momentum thickness
edge_velocity = np.ones_like(x)  # Edge velocity (normalized)

# Predict transition
result = predictor.predict_transition(x, h, re_theta, edge_velocity)

# Check results
if result['transition_occurred']:
    print(f"Transition at x/c: {result['transition_x']:.4f}")
    print(f"Transition Reynolds number: {re_theta[result['transition_index']]:.0f}")
else:
    print("No transition detected")

# Get n-factor distribution
n_factor = result['n_factor']
```

### Turbulent Boundary Layer Closure

```python
from pymises.physics.turbulent import turbulent_skin_friction, turbulent_dissipation

# Boundary layer parameters
re_theta = 1000.0  # Reynolds number based on momentum thickness
h = 1.5            # Shape parameter
h_k = 1.3          # Kinematic shape parameter
m_e = 0.2          # Edge Mach number

# Calculate skin friction coefficient
cf = turbulent_skin_friction(re_theta, h_k, m_e)

# Calculate dissipation coefficient
c_tau = 0.03  # Reynolds stress coefficient (from lag equation)
c_tau_eq = 0.035  # Equilibrium Reynolds stress coefficient
u_s = 0.7     # Shear parameter
cd = turbulent_dissipation(cf, c_tau, u_s)

print(f"Skin friction coefficient: {cf:.6f}")
print(f"Dissipation coefficient: {cd:.6f}")
```

### Artificial Dissipation

```python
from pymises.physics.dissipation import ArtificialDissipation

# Create artificial dissipation model
dissipation = ArtificialDissipation(
    threshold_mach=0.95,
    dissipation_coeff=0.5
)

# Flow field data
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

For more detailed descriptions of the physical models, refer to:
- Drela, M., "MISES Implementation of Modified Abu-Ghannam/Shaw Transition Criterion", MIT ACDL Report, 1998.
- Drela, M., "A Two-Dimensional Viscous Aerodynamic Design and Analysis Code", AIAA Paper 86-0424, 1986.
- Drela, M., "Implicit Implementation of the Full e^n Transition Criterion", AIAA Paper 2003-4066, 2003.
