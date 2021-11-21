#!/usr/bin/env python3

# ------------------------------------
# test (plot) ideal gas functions
# ------------------------------------

import numpy as np
import unyt
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from gas_functions import mean_molecular_weight
from gas_functions import gas_temperature
from gas_functions import internal_energy


def test_mean_molecular_weight():
    """
    Create a figure for every H+ ionization state
    Then plot for every He+ and He++ ionization state in
    each figure an image for H and He along the axes
    """

    nbins = 100

    # fully neutral gas
    # ---------------------
    Hp = 0.0
    Hep = 0.0
    Hepp = 0.0

    fully_neutral = np.empty((nbins, nbins))
    fully_neutral[:, :] = np.nan

    for i, H0 in enumerate(np.linspace(0, 1, nbins)):
        if Hp + Hep + Hepp + H0 > 1.0:
            break

        He0 = 1.0 - Hp - Hep - Hepp - H0
        j = int(He0 * nbins + 0.5) - 1

        mu = mean_molecular_weight(H0, Hp, He0, Hep, Hepp)
        fully_neutral[i, j] = mu

    # slightly ionized
    # ------------------------

    slightly_ionized = np.empty((nbins, nbins))
    slightly_ionized[:, :] = np.nan

    for ion in np.arange(0, 0.35, 0.01):
        Hp = ion
        Hep = ion
        Hepp = ion

        for i, H0 in enumerate(np.linspace(0, 1, nbins)):
            if Hp + Hep + Hepp + H0 > 1.0:
                break

            He0 = 1.0 - Hp - Hep - Hepp - H0
            j = int(He0 * nbins + 0.5) - 1

            mu = mean_molecular_weight(H0, Hp, He0, Hep, Hepp)
            slightly_ionized[i, j] = mu

    # increasingly ionized hydrogen
    # -------------------------------
    Hep = 0.0
    Hepp = 0.0

    increasingly_ionized_hydrogen = np.empty((nbins, nbins))
    increasingly_ionized_hydrogen[:, :] = np.nan

    for i, H0 in enumerate(np.linspace(0, 1, nbins)):
        for j, He0 in enumerate(np.linspace(0, 1, nbins)):
            if Hep + Hepp + H0 + He0 > 1.0:
                break

            Hp = 1.0 - He0 - Hep - Hepp - H0

            mu = mean_molecular_weight(H0, Hp, He0, Hep, Hepp)
            increasingly_ionized_hydrogen[i, j] = mu

    # increasingly ionized Helium
    # --------------------------------
    H0 = 0.0
    Hepp = 0.0

    increasingly_ionized_helium = np.empty((nbins, nbins))
    increasingly_ionized_helium[:, :] = np.nan

    for i, He0 in enumerate(np.linspace(0, 1, nbins)):
        for j, Hp in enumerate(np.linspace(0, 1, nbins)):
            if He0 + Hp + Hepp + H0 > 1.0:
                break

            Hep = 1.0 - He0 - Hepp - H0 - Hp

            mu = mean_molecular_weight(H0, Hp, He0, Hep, Hepp)
            increasingly_ionized_helium[i, j] = mu

    # Make plots!
    # -----------------------

    imshow_kwargs = {"extent": (0, 1, 0, 1), "origin": "lower"}

    fig = plt.figure(dpi=200, figsize=(10, 10))
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4)

    im1 = ax1.imshow(fully_neutral.T, **imshow_kwargs)
    ax1.set_title("Fully neutral gas")
    ax1.set_xlabel("mass fraction H0")
    ax1.set_ylabel("mass fraction He0")
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im1, cax=cax, orientation="vertical")

    im2 = ax2.imshow(slightly_ionized.T, **imshow_kwargs)
    ax2.set_title(
        "slightly ionized gas: $X_{H+} = X_{He+} = X_{He++} = n \\times 0.01$"
    )
    ax2.set_xlabel("mass fraction H0")
    ax2.set_ylabel("mass fraction He0")
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im2, cax=cax, orientation="vertical")

    im3 = ax3.imshow(increasingly_ionized_hydrogen.T, **imshow_kwargs)
    ax3.set_title("$X_{He} = X_{He+} = X_{He++} = 0$, $X_{H+} = 1 - X_{H0} - X_{He0}$")
    ax3.set_xlabel("mass fraction H0")
    ax3.set_ylabel("mass fraction He0")
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im3, cax=cax, orientation="vertical")

    im4 = ax4.imshow(increasingly_ionized_helium.T, **imshow_kwargs)
    ax4.set_title("$X_{H0} = X_{He++} = 0$, $X_{He+} = 1 - X_{He0} - X_{H+}$")
    ax4.set_xlabel("mass fraction He")
    ax4.set_ylabel("mass fraction H+")
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im4, cax=cax, orientation="vertical")

    #  plt.show()
    plt.tight_layout()
    plt.savefig("test_gas_functions.png")
    return


#  def test_gas_temperature():
#      """
#
#      """
#
#      nsamples = 100
#      Temps = np.random.uniform(0, 1e7, nsamples) * unyt.K
#      mus = np.random.uniform(0, 4, nsamples)
#
#      for i in range(nsamples):
#          T = Temps[i]
#          mu = mus[i]
#
#          u_calc = internal_energy(T, mu)
#          print(u_calc, T, mu)
#


if __name__ == "__main__":
    test_mean_molecular_weight()
    #  test_gas_temperature()
