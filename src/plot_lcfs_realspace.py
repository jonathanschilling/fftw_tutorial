#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 17:53:03 2021

@author: jonathan
"""

import os
import numpy as np
from mayavi import mlab

folder = ".."

# load output from FFTW example 'app_flux_surface.c'
lcfs_R = np.loadtxt(os.path.join(folder, "lcfs_R.dat"))
lcfs_Z = np.loadtxt(os.path.join(folder, "lcfs_Z.dat"))
n_zeta, n_theta = np.shape(lcfs_R)

# assume W7-X
nfp = 5

# regular grid in phi = zeta/nfp
phiGrid = 2.0*np.pi*np.arange(n_zeta)/(n_zeta*nfp)

# compute X,Y from R,phi
lcfs_X = np.zeros([n_zeta, n_theta])
lcfs_Y = np.zeros([n_zeta, n_theta])
for i in range(n_zeta):
    for j in range(n_theta):
        lcfs_X[i,j] = lcfs_R[i,j]*np.cos(phiGrid[i])
        lcfs_Y[i,j] = lcfs_R[i,j]*np.sin(phiGrid[i])

# plot real-space geometry of flux surface using Mayavi
mlab.mesh(lcfs_X, lcfs_Y, lcfs_Z)
mlab.show()
