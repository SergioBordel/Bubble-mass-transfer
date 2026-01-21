"""
bubble_kla
==========

A Python library for gas–liquid mass transfer modeling in bubble column
reactors, including hydrodynamics, physical mass transfer, and
reaction-enhanced mass transfer.

This module provides tools to:
- Compute bubble terminal velocities
- Estimate gas holdup
- Calculate kL and kLa
- Solve reaction–diffusion enhancement factors
- Calibrate bubble size and reaction rate constants from experimental data

License
-------
MIT License

Copyright (c) 2026

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions...
"""

# ============================================================================
# Imports
# ============================================================================

import math
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import root_scalar, minimize_scalar
from scipy.integrate import solve_bvp, simps

# ============================================================================
# Physical properties
# ============================================================================

def mu(T):
    """
    Dynamic viscosity of water as a function of temperature.

    Parameters
    ----------
    T : float
        Temperature in Kelvin.

    Returns
    -------
    float
        Dynamic viscosity of water (Pa·s).

    References
    ----------
    Reid, Prausnitz, Poling (1987),
    The Properties of Gases and Liquids.
    """
    mu = 1.856e-14 * math.exp((4209 / float(T)) + 0.04527 * T - 3.376e-5 * T**2)
    return mu


# ============================================================================
# Hydrodynamics
# ============================================================================

def ut(rho, mu, g, Rb, sigma):
    """
    Terminal rise velocity of a gas bubble in a liquid.

    Parameters
    ----------
    rho : float
        Liquid density (kg/m³)
    mu : float
        Liquid dynamic viscosity (Pa·s)
    g : float
        Gravitational acceleration (m/s²)
    Rb : float
        Bubble radius (m)
    sigma : float
        Surface tension (N/m)

    Returns
    -------
    float
        Terminal bubble velocity (m/s)

    Notes
    -----
    Based on correlations from:
    Park et al., Nuclear Engineering and Technology 49 (2017).
    Better predictions have been observed using the u_vis instead of the expansion proposed by Park et. al. This can be modified by the user.
    """
    u_vis = (g * rho * Rb**2) / (3 * mu)
    #u_in = 0.14425 * (g**(5/6)) * ((rho / mu)**(2/3)) * ((2 * Rb)**(3/2))
    #u_sp = np.sqrt(2.14 * sigma / (rho * 2 * Rb) + 0.505 * g * 2 * Rb)

    #u = min(u_vis, u_in, u_sp)
    return u_vis


# ============================================================================
# Physical mass transfer
# ============================================================================

def holdup(Qg, Vl, H0, u):
    """
    Gas holdup in a bubble column.

    Parameters
    ----------
    Qg : float
        Gas flow rate (L/s)
    Vl : float
        Liquid volume (L)
    H0 : float
        Liquid height (m)
    u : float
        Bubble rise velocity (m/s)

    Returns
    -------
    float
        Gas holdup (dimensionless)
    """
    epsilon = (Qg / Vl) * (H0 / u)
    return epsilon


def kl(D, u, Rb):
    """
    Liquid-side mass transfer coefficient for a single bubble.

    Parameters
    ----------
    D : float
        Molecular diffusivity (m²/s)
    u : float
        Bubble rise velocity (m/s)
    Rb : float
        Bubble radius (m)

    Returns
    -------
    float
        Liquid-side mass transfer coefficient (m/s)
    """
    kl_v = np.sqrt((2 * D * u) / (3 * Rb * np.pi))
    return kl_v


def kla(D, rho, mu, g, Qg, Vl, H0, Rb, sigma):
    """
    Volumetric mass transfer coefficient (kLa).

    Parameters
    ----------
    D : float
        Diffusivity (m²/s)
    rho : float
        Density (kg/m³)
    mu : float
        Viscosity (Pa·s)
    g : float
        Gravity (m/s²)
    Qg : float
        Gas flow rate (L/s)
    Vl : float
        Liquid volume (L)
    H0 : float
        Liquid height (m)
    Rb : float
        Bubble radius (m)
    sigma : float
        Surface tension (N/m)

    Returns
    -------
    float
        Volumetric mass transfer coefficient (1/s)
    """
    u = ut(rho, mu, g, Rb, sigma)
    epsilon = holdup(Qg, Vl, H0, u)
    kl_val = kl(D, u, Rb)
    return kl_val * (3 / Rb) * (epsilon / (1 - epsilon))


# ============================================================================
# Reaction–diffusion model
# ============================================================================

def ode_system(x, Y, k):
    """
    Converts the second-order ODE into two first-order ODEs.
    """
    f, df_dx = Y
    d2f_dx2 = -2 * x * df_dx + k * f
    return np.vstack([df_dx, d2f_dx2])


def bc(Y_a, Y_b):
    """
    Boundary conditions: f(0) = 1, f(∞) = 0 (approximated).
    """
    return np.array([Y_a[0] - 1, Y_b[0]])


def solve_ode(k, x_max=5, num_points=100):
    """
    Solves the reaction–diffusion boundary value problem.
    """
    x = np.linspace(0, x_max, num_points)
    Y_guess = np.zeros((2, x.shape[0]))
    Y_guess[0] = np.exp(-x**2)

    sol = solve_bvp(lambda x, Y: ode_system(x, Y, k), bc, x, Y_guess)

    if sol.success:
        return sol.x, sol.y[0]
    else:
        raise RuntimeError(f"ODE solver failed for k={k}")


def compute_integral(k, x_max=5, num_points=100):
    """
    Computes the integral term appearing in the reaction enhancement factor.
    """
    if k < 3000:
        x_vals, f_vals = solve_ode(k, x_max, num_points)
        integral_value = simps(f_vals, x_vals)
        result = (2 + k) * integral_value
    else:
        result = np.sqrt(k)
    return result


# ============================================================================
# Reaction-enhanced mass transfer
# ============================================================================

def fi(t):
    return (-3) * np.cos(t) + (np.cos(t))**3 + 2


def de(t, D, u, Rb):
    d = np.sqrt((8 / 3) * (D * Rb / u) * fi(t)) / (np.sin(t)**2)
    return d


def Integrand(t, D, u, Rb, k):
    d = de(t, D, u, Rb)
    In = np.sin(t) * compute_integral((k / D) * d**2) / d
    return In


def kl_r(D, u, Rb, k):
    """
    Liquid-side mass transfer coefficient with chemical reaction.
    """
    t_values = np.linspace(0, np.pi, 100)
    y = np.zeros_like(t_values)
    y[1:-1] = np.array([Integrand(ti, D, u, Rb, k)
                        for ti in t_values[1:-1]])

    kl = 0.5 * D * np.trapz(y, t_values)
    return kl


def kla_r(D, rho, mu, g, Qg, Vl, H0, Rb, k, sigma):
    """
    Volumetric mass transfer coefficient with chemical reaction.
    """
    u = ut(rho, mu, g, Rb, sigma)
    epsilon = holdup(Qg, Vl, H0, u)
    kl_val = kl_r(D, u, Rb, k)
    return kl_val * (3 / Rb) * (epsilon / (1 - epsilon))


# ============================================================================
# Calibration utilities
# ============================================================================

def solve_for_Rb(target_kla, D, rho, mu, g, Qg, Vl, H0,
                 sigma, Rb_min=1e-4, Rb_max=1e-1):
    """
    Solves kla(Rb) = target_kla for bubble radius.
    """
    def f(Rb):
        return kla(D, rho, mu, g, Qg, Vl, H0, Rb, sigma) - target_kla

    sol = root_scalar(f, bracket=[Rb_min, Rb_max], method='bisect')

    if not sol.converged:
        raise RuntimeError("Root finding did not converge.")

    return sol.root


def error(D, rho, mu, g, Qg, Vl, H0, k, sigma):
    """
    Least-squares error for calibration of reaction rate constant k.
    """
    target1 = 1.16 * 95 / 3600
    target2 = 1.16 * 61 / 3600
    target3 = 1.16 * 40 / 3600

    e3 = (kla_r(D, rho, mu, g, Qg, Vl, H0, 0.000507, k, sigma) / target3 - 1)**2
    e2 = (kla_r(D, rho, mu, g, Qg, Vl, H0, 0.000406, k, sigma) / target2 - 1)**2
    e1 = (kla_r(D, rho, mu, g, Qg, Vl, H0, 0.00031,  k, sigma) / target1 - 1)**2

    return e1 + e2 + e3


def find_k(D, rho, mu, g, Qg, Vl, H0, sigma):
    """
    Finds optimal reaction rate constant k by minimizing the calibration error.
    """
    res = minimize_scalar(
        lambda k: error(D, rho, mu, g, Qg, Vl, H0, k, sigma),
        bounds=(100, 2000),
        method='bounded'
    )
    return res.x, res.fun

