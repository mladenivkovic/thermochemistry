#!/usr/bin/env python3

# ----------------------------------------------------
# Given a photon number and a minimal frequency,
# compute the corresponding luminosity in units of
# stellar luminosities assuming a blackbody spectrum
# ----------------------------------------------------

import numpy as np
from matplotlib import pyplot as plt
import unyt
import scipy.integrate as integrate

# temperature for blackbody spectrum
T = 1e5  # K
# photon number emission rate you want
Ndot = 5e48  # s^-1
# lowest frequency
nu_min = 3.288e15  # Hz

kB = 1.3806488e-23  # J/K
h_planck = 6.62606957e-34  # Js
c = 299792458.0  # m/s
L_Sol = 3.826e26  # J/s


def B_nu(nu, T, kB, h_planck, c):
    """
    Return the blackbody energy density for
    a temperature T and frequency nu
    """
    res = 2.0 * h_planck * nu ** 3 / c ** 2 / (np.exp(h_planck * nu / kB / T) - 1.0)
    return res


def B_nu_over_h_nu(nu, T, kB, h_planck, c):
    return B_nu(nu, T, kB, h_planck, c) / (h_planck * nu)


def nu_peak(T, kB, h_planck):
    """
    Return the (approximate) frequency where the peak of the
    blackbody energy density spectrum should be
    """
    return 2.82144 * kB * T / h_planck


nu_max = 10 * nu_peak(T, kB, h_planck)
print("nu_min: {0:12.3e} nu_peak: {1:12.3e} [Hz]".format(nu_min, nu_max))

E_scipy_numax = integrate.quad(B_nu, nu_min, nu_max, args=(T, kB, h_planck, c))
N_scipy_numax = integrate.quad(
    B_nu_over_h_nu, nu_min, nu_max, args=(T, kB, h_planck, c)
)
print(
    "Luminosity up to  10 nu_peak: {0:12.6e} [L_Sol]".format(
        E_scipy_numax[0] / N_scipy_numax[0] * Ndot / L_Sol
    )
)

nu_max = 100 * nu_peak(T, kB, h_planck)
E_scipy_numax = integrate.quad(B_nu, nu_min, nu_max, args=(T, kB, h_planck, c))
N_scipy_numax = integrate.quad(
    B_nu_over_h_nu, nu_min, nu_max, args=(T, kB, h_planck, c)
)
print(
    "Luminosity up to 100 nu_peak: {0:12.6e} [L_Sol]".format(
        E_scipy_numax[0] / N_scipy_numax[0] * Ndot / L_Sol
    )
)
