# bubble_kla

**A Python library for gas–liquid mass transfer modeling in bubble column reactors.**

`bubble_kla` is a Python library designed to model gas–liquid mass transfer in bubble column reactors. It combines hydrodynamic correlations, physical mass transfer models, and reaction-enhanced mass transfer calculations to estimate volumetric mass transfer coefficients.

The library also includes numerical tools for solving reaction–diffusion equations and calibrating model parameters using experimental data.

## Features

- **Bubble hydrodynamics**
  - Calculation of terminal bubble rise velocity.
  - Estimation of gas holdup in bubble columns.

- **Physical mass transfer**
  - Calculation of the liquid-side mass transfer coefficient (`kL`).
  - Calculation of the volumetric mass transfer coefficient (`kLa`).

- **Reaction-enhanced mass transfer**
  - Numerical solution of reaction–diffusion boundary value problems.
  - Calculation of reaction enhancement effects.
  - Estimation of the liquid-side mass transfer coefficient in the presence of chemical reactions (`kL,r`).
  - Calculation of reaction-enhanced volumetric mass transfer coefficients (`kLa,r`).

- **Parameter calibration**
  - Estimation of bubble radius from experimental `kLa` measurements.
  - Calibration of reaction rate constants by minimizing the difference between experimental and predicted mass transfer coefficients.

- **Physical properties**
  - Temperature-dependent dynamic viscosity of water.

## Installation

Clone the repository:

```bash
git clone https://github.com/YOUR_USERNAME/bubble_kla.git
cd bubble_kla
```

Install the required dependencies:

```bash
pip install numpy scipy matplotlib
```

Alternatively, if a `requirements.txt` file is available:

```bash
pip install -r requirements.txt
```

## Dependencies

The library relies on the following Python packages:

| Package | Purpose |
|---|---|
| `numpy` | Numerical calculations and array operations |
| `scipy` | Numerical integration, boundary value problems, and optimization |
| `matplotlib` | Plotting and visualization |

The library uses standard Python modules such as `math` for mathematical operations.

## Theoretical Background

### 1. Physical Properties

The dynamic viscosity of water can be calculated as a function of temperature:

```python
mu(T)
```

where `T` is the temperature in Kelvin.

The correlation is based on Reid, Prausnitz, and Poling (1987), *The Properties of Gases and Liquids*.

### 2. Bubble Hydrodynamics

The terminal rise velocity of a bubble is estimated using:

```python
ut(rho, mu, g, Rb, sigma)
```

The current implementation uses the viscous-regime expression:

\[
u_t = \frac{g\rho R_b^2}{3\mu}
\]

where:

- \(\rho\): Liquid density (kg/m³).
- \(\mu\): Liquid dynamic viscosity (Pa·s).
- \(g\): Gravitational acceleration (m/s²).
- \(R_b\): Bubble radius (m).
- \(\sigma\): Surface tension (N/m).

The correlation is based on the hydrodynamic modeling discussed by Park et al. (2017).

### 3. Gas Holdup

The gas holdup is estimated using:

```python
holdup(Qg, Vl, H0, u)
```

\[
\varepsilon_g = \frac{Q_g}{V_l}\frac{H_0}{u_t}
\]

where:

- \(Q_g\): Gas flow rate (L/s).
- \(V_l\): Liquid volume (L).
- \(H_0\): Initial liquid height (m).
- \(u_t\): Bubble rise velocity (m/s).

The resulting gas holdup is dimensionless.

### 4. Physical Mass Transfer

The liquid-side mass transfer coefficient for an individual bubble is calculated using:

```python
kl(D, u, Rb)
```

\[
k_L = \sqrt{\frac{2Du}{3\pi R_b}}
\]

where \(D\) is the molecular diffusivity (m²/s).

The volumetric mass transfer coefficient is calculated using:

```python
kla(D, rho, mu, g, Qg, Vl, H0, Rb, sigma)
```

\[
k_La = k_L \frac{3}{R_b}
\frac{\varepsilon_g}{1-\varepsilon_g}
\]

### 5. Reaction–Diffusion Model

The library solves a second-order reaction–diffusion differential equation:

\[
\frac{d^2f}{dx^2}+2x\frac{df}{dx}-kf=0
\]

subject to the boundary conditions:

\[
f(0)=1
\]

\[
f(\infty)=0
\]

The numerical solution is obtained using SciPy's boundary value problem solver.

The main functions are:

- `ode_system()`: Defines the differential equation as a system of first-order equations.
- `bc()`: Defines the boundary conditions.
- `solve_ode()`: Solves the boundary value problem.
- `compute_integral()`: Evaluates the integral contribution to the reaction enhancement factor.

### 6. Reaction-Enhanced Mass Transfer

The library accounts for chemical reaction effects on mass transfer through the following functions:

| Function | Description |
|---|---|
| `fi(t)` | Auxiliary angular function |
| `de(t, D, u, Rb)` | Calculates the characteristic diffusion length |
| `Integrand(t, D, u, Rb, k)` | Evaluates the reaction-enhanced integration term |
| `kl_r(D, u, Rb, k)` | Calculates the liquid-side mass transfer coefficient with reaction |
| `kla_r(...)` | Calculates the volumetric mass transfer coefficient with reaction |

The resulting volumetric coefficient is:

\[
k_La_r = k_{L,r}\frac{3}{R_b}
\frac{\varepsilon_g}{1-\varepsilon_g}
\]

where \(k_{L,r}\) represents the reaction-enhanced liquid-side mass transfer coefficient.

## Usage

### Example 1: Calculate Physical Mass Transfer

```python
import bubble_kla as bk

# Physical properties
D = 2.2e-9       # Molecular diffusivity (m²/s)
rho = 997        # Liquid density (kg/m³)
T = 298.15       # Temperature (K)
g = 9.81         # Gravitational acceleration (m/s²)
sigma = 0.072    # Surface tension (N/m)

# Reactor parameters
Qg = 0.01        # Gas flow rate (L/s)
Vl = 1.0         # Liquid volume (L)
H0 = 0.15        # Liquid height (m)
Rb = 1e-3        # Bubble radius (m)

# Calculate water viscosity
viscosity = bk.mu(T)

# Bubble rise velocity
u = bk.ut(rho, viscosity, g, Rb, sigma)

# Gas holdup
epsilon = bk.holdup(Qg, Vl, H0, u)

# Liquid-side mass transfer coefficient
kL = bk.kl(D, u, Rb)

# Volumetric mass transfer coefficient
kLa = bk.kla(D, rho, viscosity, g, Qg, Vl, H0, Rb, sigma)

print(f"Bubble rise velocity: {u:.4e} m/s")
print(f"Gas holdup: {epsilon:.4f}")
print(f"kL: {kL:.4e} m/s")
print(f"kLa: {kLa:.4e} 1/s")
```

### Example 2: Calculate Reaction-Enhanced Mass Transfer

```python
import bubble_kla as bk

# Physical properties
D = 2.2e-9
rho = 997
T = 298.15
g = 9.81
sigma = 0.072

# Reactor parameters
Qg = 0.01
Vl = 1.0
H0 = 0.15
Rb = 1e-3

# Reaction rate constant
k = 1000

# Calculate viscosity
viscosity = bk.mu(T)

# Reaction-enhanced liquid-side coefficient
u = bk.ut(rho, viscosity, g, Rb, sigma)
kL_reaction = bk.kl_r(D, u, Rb, k)

# Reaction-enhanced volumetric coefficient
kLa_reaction = bk.kla_r(
    D, rho, viscosity, g,
    Qg, Vl, H0, Rb, k, sigma
)

print(f"Reaction-enhanced kL: {kL_reaction:.4e} m/s")
print(f"Reaction-enhanced kLa: {kLa_reaction:.4e} 1/s")
```

### Example 3: Estimate Bubble Radius from Experimental Data

The library provides two functions for estimating bubble radius from experimental mass transfer measurements.

**Without chemical reaction:**

```python
Rb = bk.solve_for_Rb(
    target_kla=0.05,
    D=D,
    rho=rho,
    mu=viscosity,
    g=g,
    Qg=Qg,
    Vl=Vl,
    H0=H0,
    sigma=sigma
)

print(f"Estimated bubble radius: {Rb:.4e} m")
```

**With chemical reaction:**

```python
Rb = bk.find_Rb_r(
    kla_exp=0.05,
    D=D,
    rho=rho,
    mu=viscosity,
    g=g,
    Qg=Qg,
    Vl=Vl,
    H0=H0,
    k=1000,
    sigma=sigma
)

print(f"Estimated bubble radius: {Rb:.4e} m")
```

### Example 4: Calibrate the Reaction Rate Constant

The library can estimate an optimal reaction rate constant by minimizing the difference between experimental and predicted reaction-enhanced mass transfer coefficients.

```python
k_opt, error_value = bk.find_k(
    D=D,
    rho=rho,
    mu=viscosity,
    g=g,
    Qg=Qg,
    Vl=Vl,
    H0=H0,
    sigma=sigma
)

print(f"Optimal reaction rate constant: {k_opt:.2f}")
print(f"Calibration error: {error_value:.4e}")
```

## API Reference

### Physical Properties

| Function | Description |
|---|---|
| `mu(T)` | Calculates water dynamic viscosity as a function of temperature. |

### Hydrodynamics

| Function | Description |
|---|---|
| `ut(rho, mu, g, Rb, sigma)` | Calculates terminal bubble rise velocity. |

### Physical Mass Transfer

| Function | Description |
|---|---|
| `holdup(Qg, Vl, H0, u)` | Calculates gas holdup. |
| `kl(D, u, Rb)` | Calculates the liquid-side mass transfer coefficient. |
| `kla(D, rho, mu, g, Qg, Vl, H0, Rb, sigma)` | Calculates the volumetric mass transfer coefficient. |

### Reaction–Diffusion

| Function | Description |
|---|---|
| `ode_system(x, Y, k)` | Defines the reaction–diffusion differential equation. |
| `bc(Y_a, Y_b)` | Defines the boundary conditions. |
| `solve_ode(k, x_max=5, num_points=100)` | Solves the reaction–diffusion boundary value problem. |
| `compute_integral(k, x_max=5, num_points=100)` | Calculates the integral contribution to the enhancement factor. |

### Reaction-Enhanced Mass Transfer

| Function | Description |
|---|---|
| `fi(t)` | Auxiliary angular function. |
| `de(t, D, u, Rb)` | Calculates the characteristic diffusion length. |
| `Integrand(t, D, u, Rb, k)` | Calculates the integration term. |
| `kl_r(D, u, Rb, k)` | Calculates the reaction-enhanced liquid-side coefficient. |
| `kla_r(D, rho, mu, g, Qg, Vl, H0, Rb, k, sigma)` | Calculates the reaction-enhanced volumetric coefficient. |

### Calibration

| Function | Description |
|---|---|
| `find_Rb_r(...)` | Estimates bubble radius using experimental reaction-enhanced `kLa` data. |
| `solve_for_Rb(...)` | Estimates bubble radius using physical `kLa` data. |
| `error(...)` | Calculates the calibration error for a reaction rate constant. |
| `find_k(...)` | Estimates the optimal reaction rate constant. |

## Units and Assumptions

Unless otherwise specified, the library uses SI units for physical properties and transport parameters.

| Parameter | Symbol | Unit |
|---|---|---|
| Temperature | \(T\) | K |
| Liquid density | \(\rho\) | kg/m³ |
| Dynamic viscosity | \(\mu\) | Pa·s |
| Molecular diffusivity | \(D\) | m²/s |
| Gravitational acceleration | \(g\) | m/s² |
| Bubble radius | \(R_b\) | m |
| Surface tension | \(\sigma\) | N/m |
| Gas flow rate | \(Q_g\) | L/s |
| Liquid volume | \(V_l\) | L |
| Liquid height | \(H_0\) | m |
| Reaction rate constant | \(k\) | Model-dependent |
| Liquid-side mass transfer coefficient | \(k_L\) | m/s |
| Volumetric mass transfer coefficient | \(k_La\) | s⁻¹ |

**Important:** Gas flow rate and liquid volume are both expressed in liters in the gas holdup calculation, so their ratio is dimensionally consistent with the remaining terms.

The models rely on simplified hydrodynamic and mass transfer correlations. Their applicability should be assessed for the specific reactor geometry, fluid properties, bubble size distribution, and operating conditions under investigation.

## References

- Reid, R. C., Prausnitz, J. M., & Poling, B. E. (1987). *The Properties of Gases and Liquids*.
- Park et al. (2017). *Nuclear Engineering and Technology*, 49.

## Limitations and Future Development

Potential extensions include:

- Incorporation of additional bubble rise velocity correlations.
- Modeling of bubble size distributions.
- Extension to non-Newtonian liquids.
- Validation against a broader range of experimental datasets.
- Support for different reaction kinetics and operating conditions.
- Improved numerical handling of limiting cases in the reaction–diffusion model.

## License

This project is distributed under the MIT License.

Copyright (c) 2026.

See the `LICENSE` file for the complete license text.
