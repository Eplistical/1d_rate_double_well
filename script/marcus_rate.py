#!/usr/bin/env python3

import sys
import numpy as np
from scipy import integrate

omega = 2e-4
mass = 2000.0
g = 20.6097
kT = 9.500432590749929e-4
dG = -0.0038;


def fermi(x):
    """fermi function"""
    return 1.0 / (1.0 + np.exp(x))


def Marcus_coef():
    """calculate Marcus rate constant"""
    Er = 0.5 * mass * omega**2 * g**2
    _4ErkT = 4 * Er * kT
    def integrand_1_to_0(e):
        return (1 - fermi(e / kT)) * np.exp(-(Er - dG + e)**2 / _4ErkT) / np.sqrt(_4ErkT * np.pi)

    def integrand_0_to_1(e):
        return fermi(e / kT) * np.exp(-(Er + dG - e)**2 / _4ErkT) / np.sqrt(_4ErkT * np.pi)

    k10, err10 = integrate.quad(integrand_1_to_0, -np.inf, np.inf, epsabs=1e-12)
    k01, err01 = integrate.quad(integrand_0_to_1, -np.inf, np.inf, epsabs=1e-12)
    return k01 + k10



if __name__ == '__main__':
    Glist = [6.4e-3, 1.6e-3, 4e-4, 1e-4]
    coef = Marcus_coef()
    print('coef = ', coef)
    for G in Glist:
        print("%.4e %.4f" % (G, np.log10(coef * G)))
