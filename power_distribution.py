# Calculates, graphs angular power distribution of emitted SPR


import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors


# Get velocity as a fraction of c from electron energy in GeV
def calc_beta(E):
    # c = 1  # Speed of light
    m = 5.11e-4  # Electron mass in GeV
    # return np.sqrt(1 - (m**2 * c**4) / E**2)
    return np.sqrt(1 - (m ** 2) / E ** 2)


def calc_lambda(theta, phi, D, beta, n):
    return (D / n) * (1 / beta - np.cos(theta) * np.sin(np.pi / 2 - phi))


def calc_lambda_e(theta, phi, beta, gamma, lambda_this):
    term1 = lambda_this / (2 * np.pi)
    term2 = (beta * gamma) / np.sqrt(1 + beta**2 * gamma**2 * np.sin(theta)**2 * np.sin(phi)**2)
    return term1 * term2


# Calculate R_n when phi = 0
# theta - Angle 'up' wrt. beam
# phi - Angle 'left/right' wrt. beam (pi / 2 is parallel)
# Equations from https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.8.091301
def calc_R2(theta, phi, beta, gamma, N, L, n, alpha):

    D = L / N # grating period
    h = D * np.tan(alpha) # grating height

    lambda_this = calc_lambda(theta, phi, D, beta, n)
    lambda_e = calc_lambda_e(theta, phi, beta, gamma, lambda_this)

    k = 2 * np.pi / lambda_this # wavenumber
    k_x = k * np.cos(phi) * np.sin(theta)
    k_y = k * np.sin(phi)
    k_z = k * np.cos(phi) * np.cos(theta)

    D_j = k / beta - k_z - k_x * np.tan(alpha) - 1j * np.tan(alpha) / lambda_e

    iDD_j = 1j *D * D_j

    # G_bar_term_1 = np.array([np.tan(alpha), 2j * k_y * lambda_e * np.tan(alpha), 1])
    # G_bar_term_2 = np.exp((1 / lambda_e - 1j * k_x) * h + 1j * (k / beta - k_z) * D)
    # G_bar_term_3 = (np.exp(-1j * D_j * D) - 1) / (1j * D_j * D)

    G_bar_term_1 = np.array([np.tan(alpha), 2j * k_y * lambda_e * np.tan(alpha), 1])
    G_bar_term_2 = np.exp(iDD_j)
    G_bar_term_3 = (np.exp(-1 * iDD_j) - 1) / iDD_j

    G_bar = G_bar_term_1 * G_bar_term_2 * G_bar_term_3

    eps_bar_par = np.array([np.cos(theta) * np.cos(phi), np.cos(theta) * np.sin(phi), -1 * np.sin(theta)])
    eps_bar_perp = np.array([-1 * np.sin(phi), np.cos(phi), 0])

    R2_par = np.abs(np.dot(eps_bar_par, G_bar))**2
    R2_perp = np.abs(np.dot(eps_bar_perp, G_bar))**2

    if np.isnan(R2_par):
        return 0
    else:
        return R2_par + R2_perp

    # Old
    # l = h / np.tan(alpha)
    # k = 2 * np.pi * n / (l * (1 / beta - np.cos(theta)))
    #
    # # Equation (7)
    # term1 = (1 / (np.sin(theta) - np.tan(alpha) * (1 + np.cos(theta)))) ** 2
    # term2 = np.tan(alpha) ** 2 * np.sin(theta) ** 2
    # term3 = np.sinc(k * h * np.sin(theta) / (2 * np.pi)) ** 2
    # return term1 * term2 * term3


# Calculate angular power distribution of SPR radiation along normal plane (phi = 0)
# theta - Angle 'up' wrt. beam direction (rad)
# phi - Angle 'left/right' wrt. beam
# I - Beam current (A)
# n - Diffraction order
# L - Total length of grating (m)
# N - Number of grating periods
# alpha - Blaze angle of echelle grating (rad)
# h - Height of grating (peak to trough)
# E - Beam energy (GeV)
# d - Height of beam above grating (m)
def calc_distribution(theta, phi, n, L, N, alpha, E, d):
    D = L / N  # grating period

    beta = calc_beta(E)
    gamma = 1 / np.sqrt(1 - beta ** 2)

    lambda_this = calc_lambda(theta, np.pi / 2, D, beta, n)

    h_int = lambda_this * beta * gamma / (4 * np.pi)

    term1 = 1 / 137 * np.abs(n) * N
    term2 = np.sin(theta)**2 * np.cos(phi)**2 / ( # cos instead of sin because different phi
                1 / beta - np.cos(theta) * np.cos(phi))**2
    term3 = calc_R2(theta, phi, beta, gamma, N, L, n, alpha)
    term4 = np.exp(-1 * d / h_int * np.sqrt(1 + (beta * gamma * np.sin(phi))**2))

    return term1 * term2 * term3 * term4 if term4 != 0 else 1e-100


if __name__ == '__main__':

    # filepath_img = input('Enter file path for graph: ')
    filepath_data = input('Enter file path for JSON data: ')

    n = int(input('Enter diffraction order: '))
    L = float(input('Enter total length of grating in m: '))
    N = int(input('Enter total number of grating periods: '))
    E = float(input('Enter beam energy in GeV: '))
    d = float(input('Enter height of beam above grating in m: '))
    alpha = float(input('Enter blaze angle of echelle grating in °: '))

    thetas, phis = np.meshgrid(np.linspace(0, np.pi, 501), np.linspace(-np.pi / 60, np.pi / 60, 501))
    power_dist = np.array([[calc_distribution(theta, phi[0], n, L, N, alpha * np.pi / 180, E, d) for theta in thetas[0]] for phi in phis])

    print(thetas)
    print(phis)
    print(power_dist)

    cmap = 'magma'

    thetas = thetas * 180 / np.pi
    phis = phis * 180 / np.pi

    fig, ax = plt.subplots()
    chart = ax.pcolormesh(phis, thetas, power_dist, norm=colors.LogNorm(vmin=1e-10, vmax=1e5), cmap=cmap)

    cbar = fig.colorbar(chart)
    cbar.set_label('$\\frac{dN_\\gamma}{d\\Omega}$', rotation=0, labelpad=12)

    ax.set_title("Expected angular distribution of \nSPR per electron along normal plane")

    ax.set_xlabel('$\\phi$ (°)')
    ax.set_ylabel('$\\theta$ (°)')

    plt.suptitle(f'$n = {n}$, $L = {L * 1e3}$mm, $N_ = {N}$,\n$E = {E}$ GeV, $d= {d * 1e3}$mm, $\\alpha = {alpha}°$', x=0.432, y=0.85, color='white', fontsize=10)

    plt.show()

    # fig.savefig(filepath_img)

    with open(filepath_data, 'w') as f:
        f.write(json.dumps(power_dist.tolist()))
