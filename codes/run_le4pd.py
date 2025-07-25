#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import le4pd as le4pd
#import matplotlib.pyplot as plt
#import physt
import os
import sys

path = sys.argv[1]
TOP = sys.argv[2]
XTC = sys.argv[3]
SYS = XTC.split(".")[0]
T = 300
print(path, SYS, XTC, TOP)

ftraj = le4pd.traj_from_xtc(XTC, TOP)
print(ftraj.shape[1])
NFRS = ftraj.shape[1]
N = ftraj.shape[0] // 3
NATOMS = N
print(N, NATOMS, NFRS)
covar, R, mu = le4pd.calc_Cmatrix(ftraj, path = path)
Rinv = le4pd.Rij(ftraj, path = path)

#Calculate the friction coefficients; viscosity of solvent is in units of Pa s .
fratio, sigma, fric, avfr = le4pd.fric_calc(TOP, N, NFRS, NATOMS, mu, T, intv = 2.71828, viscosity = 1e-3, fd20 = 0.0, path_to_resarea = path)

#Calculate the H, M, a, L, and the eigendecomposition of the LU matrix
H, AIHI, Q, QINV, lambda_eig = le4pd.eigendecomp(covar, Rinv, avfr, fratio, fric, path = path)
#barlist, xi, dummy = le4pd.mode_mad(ftraj, N, NFRS, Q, QINV, T, nmodes = 10, path = path, HA = True)
#tau, tau_scaled = le4pd.tau_convert(lambda_eig, sigma, barlist, T, path = './', HA = True)
lml = le4pd.LML(Q, mu, path = path, HA = True)

#np.savetxt(path + 'tau.dat', tau)
#np.savetxt(path + 'tau_scaled.dat', tau_scaled)
#np.savetxt(path + 'barriers_kcal.dat', barlist)
np.savetxt(path + 'lml.dat', lml)

with open(path + 'protname.txt', "w+") as f:
	f.write(SYS + '\n')
	f.write(str(N) + '\n')
	f.write(str(NFRS) + '\n')
	f.write(str(NATOMS) + '\n')
f.close()


np.save(path + 'Rinv.npy', Rinv)
np.save(path + 'Cmatrix.npy', covar)
np.save(path + 'Qmatrix.npy', Q)
np.save(path + 'QINVmatrix.npy', QINV)
np.savetxt(path + 'lambda_eig', lambda_eig)
np.savetxt(path + 'mu_eig', mu)
np.savetxt(path + 'sigma', np.array([sigma]))
np.savetxt(path + 'fric', fric)
np.savetxt(path + 'avfr', np.array([avfr]))
np.save(path + 'AIHImatrix.npy', AIHI)
