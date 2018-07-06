#!/usr/bin/env python3
import re
from datetime import datetime as dt
import numpy as np
import scipy.integrate as integrate
import argparse
import h5py

np.seterr(over='ignore')

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gamma0', required=True, type=float, help='coupling')
    parser.add_argument('--kT', required=True, type=float, help='temperature')
    parser.add_argument('--W', type=float, default=1.0, help='bandwidth, the band is defined as [-W, W]')
    parser.add_argument('--hmax', type=float, default=0.5, help='max h(x) to integrate')
    parser.add_argument('--hmin', type=float, default=-0.5, help='hmin h(x) to integrate')
    parser.add_argument('--Nh', type=int, default=10000, help='number of h to process')
    parser.add_argument('--integrator', type=str, default="quad", help='integrator to use, quad or naive')
    parser.add_argument('--outfile', type=str, default="inttable.h5", help='output hdf5 file')
    return parser.parse_args()


args = parse_args()
W = args.W
hmin, hmax = args.hmin, args.hmax
Nh = args.Nh
harr, dh = np.linspace(hmin, hmax, Nh, retstep=True)
outfile = args.outfile
gamma0 = args.gamma0
kT = args.kT

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
    print(args)
    start_t = dt.now()

    # calculate integral
    if args.integrator == 'quad':
        Af_arr, A2dfde_arr = quadint()
    elif args.integrator == 'naive':
        Af_arr, A2dfde_arr = naiveint()
    else:
        raise ValueError("invaild integrator type")

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
        print(('%18.6f' * 3) % (h, x, y))

    end_t = dt.now()
    print("# time elasped: " + str(end_t - start_t))
