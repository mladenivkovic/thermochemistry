#!/usr/bin/env python3


# ------------------------------------------
# Plot and test integration accuracy for
# a blackbody spectrum.
# ------------------------------------------

import numpy as np

from matplotlib import pyplot as plt
import unyt

skip_plots = False

T = 1e5 * unyt.K

kB = unyt.boltzmann_constant
h_planck = unyt.h
c = unyt.c


def BT(nu, T):
    """
    Return the blackbody energy density for
    a temperature T and frequency nu
    """
    res = (
        8.0
        * np.pi
        * h_planck
        * nu ** 3
        / c ** 3
        / (np.exp(h_planck * nu / kB / T) - 1.0)
    )
    return res.to("erg/cm**3*s")


def nu_max(T):
    """
    Return the (approximate) frequency where the peak of the
    blackbody energy density spectrum should be
    """
    return 2.82144 * kB * T / h_planck


# ---------------
# Plotting
# ---------------

if not skip_plots:

    nu = np.linspace(0, 10 * nu_max(T), 100000)  # already in units Hz
    spectrum = BT(nu, T)
    peakval = spectrum[np.isfinite(spectrum)].max()

    ionization_frequencies = [
        3.288e15 * unyt.Hz,
        5.946e15 * unyt.Hz,
        1.316e16 * unyt.Hz,
    ]
    ionization_names = ["HI ionization", "HeI ionization", "HeII ionization"]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot([1.0, 1.0], [0.0, peakval * 1.01], label=r"$\nu_{max}$", c="k")
    ax.plot(nu / nu_max(T), spectrum, label="T = {0:.1e}".format(T))
    for i in range(len(ionization_frequencies)):
        f = ionization_frequencies[i] / nu_max(T)
        if f < nu.max():
            ax.plot([f, f], [0, peakval], label=ionization_names[i])
    ax.set_ylabel(r"[erg s / cm$^3$]")
    ax.set_xlabel(r"$\nu$ / $\nu_{max}$")
    ax.legend()
    plt.show()


# --------------------------------
# Integrating the Spectrum
# --------------------------------

# NOTE: deprecated. This is not a very good way of integrating
# functions like these. Use scipy.integrate.quad instead.
# See ./convert_photon_number_to_luminosity.py


def integrate_spectrum(nu_start, nu_end, npoints, T):

    dnu = (nu_end - nu_start) / npoints

    integral = 0
    nu_current = nu_start
    left = BT(nu_start, T)
    for i in range(npoints):
        nu_current += dnu
        right = BT(nu_current, T)
        add = 0.5 * (left + right) * dnu
        if np.isfinite(add):
            integral += add
        else:
            print("Got", add, "for nu=", nu_current, "i=", i)

        left = right

    return integral


npoints = [100, 1000, 10000, 100000]  # , 1000000] #, 1000000]
nu_ends = [5, 10, 20]


results = [[0 for e in nu_ends] for n in npoints]
nu_peak = nu_max(T)

for i, n in enumerate(npoints):
    for j, end in enumerate(nu_ends):
        print("Computing", n, end)
        # Don't start at zero Hz, it'll lead to infinities
        res = integrate_spectrum(1 * unyt.Hz, end * nu_peak, n, T)
        results[i][j] = res

print()
topline = "{0:20} || {1:22.1f} | {2:22.1f} | {3:22.1f}".format(
    "nu_end/nu_peak = ", nu_ends[0], nu_ends[1], nu_ends[2]
)
print(topline)
for char in topline:
    print("-", end="")
print()

for i, n in enumerate(npoints):
    nextline = "n = {0:16d} || {1:12.6e} | {2:12.6e} | {3:12.6e}".format(
        n, results[i][0], results[i][1], results[i][2]
    )
    print(nextline)
