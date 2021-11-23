#!/usr/bin/env python3

# --------------------------------------
# Plot the ionization and recombination
# rates for various temperatures
# --------------------------------------

import numpy as np
from matplotlib import pyplot as plt

T = np.logspace(0, 11, 10000)

# Recombination rate for H+ in units of cm^3 s^-1
A_Hp = 8.40e-11 / np.sqrt(T) * (T * 1e-3) ** (-0.2) / (1.0 + (T * 1e-6) ** 0.7)
# TODO: DOC
A_d = 1.9e-3 / T ** 1.5 * np.exp(-470000.0 / T) * (1.0 + 0.3 * np.exp(-94000.0 / T))
# Recombination rate for He+ in units of cm^3 s^-1
A_Hep = 1.5e-10 / T ** 0.6353
# Recombination rate for He++ in units of cm^3 s^-1
A_Hepp = 3.36e-10 / np.sqrt(T) * (T * 1e-3) ** (-0.2) / (1.0 + (T * 1e-6) ** 0.7)
# ionization rate for H0 in units of cm^3 s^-1
G_H0 = 5.85e-11 * np.sqrt(T) * np.exp(-157809.1 / T) / (1.0 + np.sqrt(T * 1e-5))
# TODO: DOC
G_He0 = 2.38e-11 * np.sqrt(T) * np.exp(-285335.4 / T) / (1.0 + np.sqrt(T * 1e-5))
# TODO: DOC
G_Hep = 5.68e-12 * np.sqrt(T) * np.exp(-631515.0 / T) / (1.0 + np.sqrt(T * 1e-5))


fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1)
ax1.loglog(T, A_Hp, label="$\\alpha_{H^+}$")
ax1.loglog(T, A_d, label="$\\alpha_{d}$")
ax1.loglog(T, A_Hep, label="$\\alpha_{He^{+}}$")
ax1.loglog(T, A_Hepp, label="$\\alpha_{He^{++}}$")
ax1.loglog(T, G_H0, "--", label="$\\Gamma_{H^{0}}$")
ax1.loglog(T, G_He0, "--", label="$\\Gamma_{He^{0}}$")
ax1.loglog(T, G_Hep, "--", label="$\\Gamma_{He^{+}}$")
ax1.legend()
ax1.set_title("Rates")
ax1.set_xlabel("T [K]")
ax1.set_ylabel("rates")
ax1.set_ylim(1e-14, 1e-7)

ax2 = fig.add_subplot(1, 2, 2)
ax2.loglog(
    T,
    A_Hp / (A_Hp + G_H0),
    label="$\\alpha_{H^{+}} / (\\alpha_{H^{+}} + \\Gamma_{H^{0}})$",
)
ax2.loglog(
    T[T > 1e3],
    ((A_Hep + A_d) / G_He0)[T > 1e3],
    label="$(\\alpha_{He^{+}} + \\alpha_{d}) / \\Gamma_{He^{0}})$ for T > 1000 K",
)
ax2.loglog(T, G_Hep / A_Hepp, label="$(\\Gamma_{He^{+}} / \\alpha_{He^{++}})$")
ax2.legend()
G_Hep / A_Hepp

#  temp = (A_Hep + A_d) / G_He0
#  for i, t in enumerate(T):
#      print(t, temp[i])


plt.show()
