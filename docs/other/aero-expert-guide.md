PyMISES Aerodynamic Expert Guide

# 1. Mathematical Foundation

# 1.1 Euler Equations

The core of PyMISES is based on the steady, two-dimensional Euler equations formulated in integral form:

Mass Conservation:

∮ ρq·n ds = 0

Momentum Conservation:

∮(ρ(q·n)q + pn) ds = 0

Energy Conservation:

∮(ρ(q·n)h) ds = 0

Where:

ρ is density

q is velocity vector

n is unit normal vector

p is pressure

h is total enthalpy

# 1.2 Streamline-Based Discretization

The unique aspect of PyMISES is the use of an intrinsic streamline grid where:

One set of coordinate lines corresponds to streamlines

There is no convection across the corresponding faces of conservation cells

Mass flux and stagnation enthalpy are constant along streamtubes

This leads to several advantages:

Reduction from 4 to 2 unknowns per grid node

No numerical diffusion of entropy or enthalpy across streamlines

Faithful reproduction of analytic information propagation mechanisms

Natural adaptation to the flow field

# 1.3 Artificial Dissipation

For transonic flows, artificial dissipation is implemented as a bulk viscosity-like term:

λ = { 0,                    M < Mc

{ Cd(1-(Mc/M)²)/(2M²), M > Mc

Where:

Mc is threshold Mach number (typically 0.9-0.95)

Cd is dissipation coefficient (typically 0.5-1.0)

This approach provides:

Stable shock capturing

Prevention of expansion shocks

Minimal numerical diffusion in subsonic regions

# 1.4 Integral Boundary Layer Equations

The viscous effects are modeled using integral boundary layer equations:

Momentum Integral Equation:

θ(dθ/dx) = (Cf/2) - (H+2-M²ₑ)(θ/Uₑ)(dUₑ/dx)

Shape Parameter Equation:

θ(dH*/dx) = 2CD/H* - Cf/2·(2H**/H* + 1 - H) - θ/Uₑ·(dUₑ/dx)

Where:

θ is momentum thickness

δ* is displacement thickness

H = δ*/θ is shape parameter

H* is kinetic energy shape parameter

H** is density thickness shape parameter

Cf is skin friction coefficient

CD is dissipation coefficient

Ue is edge velocity

Me is edge Mach number

# 1.5 Closure Relations

# 1.5.1 Laminar Closure

For laminar flow, closure relations are based on the Falkner-Skan profile family:

H* = 1.515 + 0.076(Hk-4)²/Hk,  Hk < 4

H* = 1.515 + 0.040(Hk-4)²/Hk,  Hk > 4

Cf·Reθ = -0.067 + 0.01977(7.4-Hk)²/(Hk-1)

2CD/H* = 0.207 + 0.00205(4-Hk)⁵·⁵, Hk < 4

2CD/H* = 0.207 - 0.003(Hk-4)², Hk > 4

Where Hk is the kinematic shape parameter.

# 1.5.2 Turbulent Closure

For turbulent flow:

Cf = 0.3e⁻¹·³³Hk/(log₁₀(Reθ))¹·⁷⁴⁺⁰·³¹Hk · Fc

Where Fc is a compressibility correction.

The energy shape parameter is given by:

H* = (1.505 + 400/Reθ) / Hk·(1 + (H₀-H)²(0.04 + 0.007ln(Reθ)/Hk))

The dissipation coefficient uses a lag-entrainment approach:

CD = Cf·Us² + CT(1-Us)

Where CT is a maximum shear stress coefficient that follows a lag equation to account for non-equilibrium effects:

δ dCT/dx = KC(CT_EQ^(1/2) - CT^(1/2))

With KC = 5.6 as the lag constant and δ as the boundary layer thickness.

# 1.6 Transition Prediction

Based on Drela's modified Abu-Ghannam/Shaw transition criterion, transition is predicted using a hybrid approach combining the eⁿ method with a bypass transition model:

# 1.6.1 Amplification Factor Method

The amplification factor n is governed by:

θ(dn/dx) = f(H,Reθ) + g(H,Reθ)

Where:

f(H,Reθ) is the Orr-Sommerfeld-based TS wave growth rate

g(H,Reθ) is a bypass transition term that activates as Reθ approaches a critical value

The g term is defined as:

g(H,Reθ) = { 0             r < 0

{ A(3r² - 2r³)  0 < r < 1

{ A             r > 1

r = (1/B)((Reθ/Reθs) - 1) + 1/2

Where:

A = 0.10 (amplification rate constant)

B = 0.30 (ramp width parameter)

Reθs is the critical Reynolds number

# 1.6.2 Critical Reynolds Number

The critical Reynolds number is modified from the original Abu-Ghannam/Shaw formulation:

Reθs(H,ñcrit) = 155 + 89.0[0.25tanh(10/(H-1) - 5.5) + 1](ñcrit)^1.25

Where ñcrit is related to the turbulence level τ by:

τ' = 2.7tanh(τ/2.7)

ñcrit(τ) = -8.43 - 2.4ln(τ'/100)

This formulation addresses the ill-posedness of the original AGS criterion by:

Using H instead of λ to parameterize boundary layer instability

Adding a smooth transition from stability to instability

Accounting for viscous-inviscid interaction effects

# 1.6.3 Initial Turbulence Level

For transition modeling in varying edge velocity fields (e.g., turbine cascades), the amplification equation is modified to:

(dn̄/dx) + (dlnue/dx) = (1/θ)[f(H,Reθ) + g(H,Reθ)]

This accounts for the acceleration or deceleration of the external flow.

# 1.7 Viscous-Inviscid Coupling

Coupling is achieved through the displacement body concept:

Inviscid flow sees a body displaced from the physical body by δ*

Boundary layer equations are integrated into the global Newton system

No separate iteration between solvers

This provides robust convergence even with separation bubbles present.

# 1.8 Inverse Design Framework

The inverse design capability includes:

Full Inverse:

Pressure specified on entire airfoil

Leading and trailing edge closure constraints

Free parameters for well-posed problem

Mixed Inverse:

Pressure specified on part of airfoil

Geometry fixed elsewhere

Continuity enforced at segment boundaries

Additional free parameters for segment continuity

# 2. Numerical Solution Procedure

# 2.1 Newton Method

The coupled system of equations is solved using Newton's method:

[dF/dU]·ΔU = -F(U)

U^(v+1) = U^(v) + ΔU

Where:

F(U) is the residual vector

U is the solution vector (ρ, δn at each node)

[dF/dU] is the Jacobian matrix

v is the iteration index

# 2.2 Jacobian Structure

The Jacobian matrix has a block structure with:

Four block diagonals (A, B, C, Z) for local variables

Additional rows/columns for global variables and constraints

Typical size ~4000-10000 equations

# 2.3 Linear System Solution

The linear system is solved using block Gaussian elimination:

Forward sweep to eliminate Z and B blocks

Backward sweep to eliminate modified C blocks

Separate solution for global variables

Back-substitution for local variables

# 2.4 Convergence Acceleration

Convergence is accelerated through:

Proper under-relaxation in early iterations (r < 1)

Full Newton steps (r = 1) once near solution

Grid redistribution when needed

Careful initialization strategies

## 3. Implementation Guidelines

## 3.1 Grid Generation

Start with elliptic grid generation (Thompson method)

Initialize streamlines for incompressible flow

Allow grid adaptation during Newton iterations

Redistribute grid periodically to maintain quality

Special treatment for leading edge stagnation point

## 3.2 Boundary Conditions

Wall Boundary Conditions:

No-flow-through (streamline coincides with wall)

Displacement thickness offset in viscous mode

Far-Field Conditions:

Pressure specification based on potential flow (vortex + source + doublet)

Consistent inflow/outflow angle specification

Treatment of wake streamline

Cascade Conditions:

Periodicity in tangential direction

Specified inlet flow angle

Kutta condition at trailing edge

## 3.3 Stability Considerations

Enhanced artificial dissipation near shocks

Pressure correction to eliminate grid sawtooth modes

Careful treatment of stagnation points

Edge velocity correction in boundary layer coupling

Special attention to transition region to avoid numerical instabilities

## 3.4 Performance Metrics

Lift and Moment:

Integrated from surface pressure distribution

Consistent with far-field circulation

Drag Components:

Profile drag from momentum deficit in wake

Wave drag from entropy increase across shock

Induced drag from far-field analysis

Boundary Layer Parameters:

Displacement thickness (δ*)

Momentum thickness (θ)

Shape factor (H)

Separation locations

## 4. Validation Test Cases

The following test cases are recommended for validation:

Joukowski airfoil: Compare with exact incompressible solution

RAE 2822 (Case 6): Attached transonic flow with shock

RAE 2822 (Case 10): Transonic flow with shock-induced separation

NACA 4412 at high angle: Trailing edge separation

LA203A at low Reynolds number: Transitional separation bubbles

Supercritical cascade: For periodicity and shock capturing

## 5. Common Issues and Remedies

Convergence Difficulties:

Increase artificial dissipation (Cd)

Lower threshold Mach number (Mc)

Use stronger under-relaxation

Improve initial grid quality

Shock Instabilities:

Refine grid near shock

Increase artificial dissipation locally

Smooth pressure distribution in early iterations

Separation Problems:

Carefully tune boundary layer closure parameters

Improve grid resolution in separation region

Modify transition criteria if needed

Transition Prediction Issues:

Use the modified transition model based on H instead of λ

Verify proper coupling between boundary layer and inviscid flow

Consider appropriate transition ramp parameters (A and B)

Inverse Design Issues:

Smooth target pressure distribution

Add more degrees of freedom

Enforce reasonable geometric constraints

Start from well-converged direct solution

## 6. Advanced Features

## 6.1 Multi-Point Design

For multi-point design:

Perform inverse design at primary design point

Check performance at off-design conditions

Modify target pressure to improve off-design performance

Repeat until satisfactory compromise is achieved

## 6.2 Performance Optimization

For direct optimization:

Parameterize geometry (e.g., control points)

Define objective function (e.g., L/D, Cd at fixed Cl)

Use gradient-based methods leveraging sensitivity information from Newton solver

Apply constraints on geometry and aerodynamic parameters

## 6.3 Sensitivity Analysis

Newton method automatically provides sensitivity information:

dCl/dα: Lift curve slope

dCd/dM: Drag-divergence characteristics

dCm/dCl: Stability derivatives

These can be used for design refinement without additional solutions.