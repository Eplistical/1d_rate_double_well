#!/usr/bin/env python3
import re
from datetime import datetime as dt
import numpy as np
import scipy.integrate as integrate
import h5py

np.seterr(over='ignore')

with open('config.hpp', 'r') as f:
    txt = f.read()
    kT = float(re.findall('const double kT = (.*);', txt)[0])
    gamma0 = float(re.findall('const double gamma0 = (.*);', txt)[0])

W = 1.0
hmin, hmax = -0.5, 0.5
Nh = int(10000)
harr, dh = np.linspace(hmin, hmax, Nh, retstep=True)

outdir = '.'
outfile = outdir + '/' + 'inttable.hdf5'

def A(e, h):
    return gamma0 / ((e - h)**2 + gamma0**2 / 4)

def cal_f(E):
    return 1.0 / (np.exp(E / kT) + 1.0)

def cal_dfde(E):
    f = cal_f(E)
    return -f * (1.0 - f) / kT

def Af_integrand(e, h):
    return A(e, h) * cal_f(e)

def A2dfde_integrand(e, h):
    return A(e, h)**2 * cal_dfde(e)

def quadint():
    """calc integrals by quad"""
    Af_arr = np.zeros(harr.size, dtype=np.float64)
    A2dfde_arr = np.zeros(harr.size, dtype=np.float64)
    for i, h in enumerate(harr):
        Af_arr[i] = integrate.quad(Af_integrand, -W, W, args=(h, ))[0]
        A2dfde_arr[i] = integrate.quad(A2dfde_integrand, -np.inf, np.inf, args=(h, ))[0]
    return Af_arr, A2dfde_arr

def naiveint():
    """calc integrals by summation"""
    Af_arr = np.zeros(harr.size, dtype=np.float64)
    A2dfde_arr = np.zeros(harr.size, dtype=np.float64)
    Ne = 100000
    e_arr, de = np.linspace(-2*W, 2*W, Ne, Ne, retstep=True)
    for i, h in enumerate(harr):
        Af_arr[i] = np.sum(A(e_arr, h) * cal_f(e_arr)) * de
        A2dfde_arr[i] = np.sum(A(e_arr, h)**2 * cal_dfde(e_arr)) * de
    return Af_arr, A2dfde_arr



if __name__ == '__main__':
    logf = open('.Af_int.log', 'w+')
    print('#  kT = ', kT, file=logf)
    print('#  gamma0 = ', gamma0, file=logf)
    start_t = dt.now()

    # calculate integral
    Af_arr, A2dfde_arr = naiveint()

    f = h5py.File(outfile, 'w')
    # parameters
    para = f.create_dataset('para', data=np.float64(0.0))
    for key in ('kT', 'gamma0', 'hmax', 'hmin', 'dh', 'Nh', 'W'):
        para.attrs[key] = locals()[key]
    # data
    f.create_dataset('Af', data=Af_arr)
    f.create_dataset('A2dfde', data=A2dfde_arr)
    f.close()

    for h, x, y in zip(harr, Af_arr, A2dfde_arr):
        print(('%18.6f' * 3) % (h, x, y), file=logf)

    end_t = dt.now()
    print("# time elasped: " + str(end_t - start_t), file=logf)
