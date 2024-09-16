# Smith-Purcell Radiation Wavelength calculator and grapher
# For the AHS Physics Club

import numpy as np
import matplotlib.pyplot as plt


# Get velocity as a fraction of c from electron energy in GeV
def beta(E):
    # c = 1  # Speed of light
    m = 5.11e-4  # Electron mass in GeV
    # return np.sqrt(1 - (m**2 * c**4) / E**2)
    return np.sqrt(1 - (m**2) / E**2)


# Calculate wavelength of Smith-Purcell Radiation
# D in m
# theta, phi, in radians
# beam_energy in GeV
def lambda_this(theta, phi, D, beam_energy, n=-1):
    beta = beta(beam_energy)
    return (D / np.abs(n)) * (1 / beta - np.cos(theta) * np.sin(phi))


if __name__ == '__main__':
    D = float(input('Enter grating period in nm: ')) * 1e-9  # in m
    beam_energy = float(input('Enter beam energy in GeV: '))  # in GeV
    n = int(input('Enter diffraction index: '))  # Diffraction index

    phis, thetas = np.meshgrid(np.linspace(0, np.pi, 101), np.linspace(np.pi / -2, np.pi / 2, 101))
    lambdas = np.array([[1e9 * lambda_this(theta[0], phi, D, beam_energy, n) for theta in thetas] for phi in phis[0]])

    print(thetas)
    print(phis)

    thetas = thetas * 180 / np.pi
    phis = phis * 180 / np.pi

    cmap = 'Spectral'

    figure, axes = plt.subplots()
    chart = axes.pcolormesh(thetas, phis, lambdas, cmap=cmap)

    cbar = figure.colorbar(chart)
    cbar.set_label('Wavelength (nm)', rotation=270, labelpad=12)

    axes.set_title(f'Expected wavelength of Smith-Purcell radiation \nof diffraction order {n} \
as a function of emission angles θ and Φ \
\nwith grating period {int(D * 1e9)} nm and beam energy {beam_energy} GeV')
    axes.set_xlabel('$\\phi$ (°)')
    axes.set_ylabel('$\\theta$ (°)')

    plt.show()
