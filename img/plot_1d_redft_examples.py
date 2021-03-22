#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 18:59:42 2021

@author: jonathan
"""

import numpy as np
from numpy.random import default_rng
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

#%% REDFT00 

# compute the REDFT00 as defined in the FFTW reference manual
def redft00(arr):
    n = len(arr)
    out = np.zeros([n])
    for k in range(n):
        out[k] += arr[0]
        out[k] += (-1)**k * arr[-1]
        for j in range(1, n-1):
            out[k] += 2.0 * arr[j] * np.cos(np.pi*j*k/(n-1))
    return out

# evaluate the REDFT00 at all points given in x
def eval_redft00(arr, x):
    y = np.zeros_like(x)
    y += arr[0]
    for j in range(1,n-1):
        y += 2.0*arr[j]*np.cos(np.pi*j*x)
    y += arr[-1]*np.cos(np.pi*(n-1)*x)
    return y

# input size for REDFT00
n = 5

# logical size of equivalent DFT (for REDFT00)
N = 2*(n-1)

# random Fourier coefficients
rng = default_rng(seed=42)
r2r_in = rng.uniform(size=n)

# compute output of REDFT00
r2r_out = redft00(r2r_in)/n

# extend to equivalent DFT size
full_out_x = np.arange(n, 2*n-2)
full_out_y = np.zeros([n-2])
for i in range(n-2):
    # even symmetry around n-1
    full_out_y[i] = r2r_out[-i-2]

# sample at finer intervals
nFine = 1024
x = np.linspace(-n+1, N, nFine)
y = eval_redft00(r2r_in, x/(n-1))/n

plt.figure(figsize=(5,3))

plt.plot(x, y, '-', color='gray', linewidth=0.5)
plt.axhline(0, ls='-', color='k')
plt.axvline(0, ls='-', color='k')

# symmetry lines
plt.axvline(0, ls='--', color='r')
plt.axvline(n-1, ls='--', color='r')

# highlight array contents
plt.axvspan(0-0.25, (n-1)+0.25, alpha=0.3, color='gray')

# plot actual REDFT00 output
plt.plot(r2r_out, 'bo')

# label individual points
for i in range(n):
    plt.text(i+0.2, r2r_out[i]-0.03, chr(ord("a")+i))

# plot data extending REDFT00 to full DFT
plt.plot(full_out_x, full_out_y, 'bo')

for i in range(n-2):
    plt.text(full_out_x[i]+0.2, full_out_y[i]-0.03, chr(ord("a")+(n-2-i)))

# only integer xaxis ticks at intervals of 1
plt.gca().xaxis.set_major_locator(MaxNLocator(steps=(1,10), integer=True))
plt.xlim((-4.5, 8.5))

plt.grid(True)
plt.xlabel("j")
plt.title("REDFT00 N=%d n=%d"%(N, n))
plt.tight_layout()
plt.savefig("redft00.png")

#%% REDFT10
n = int(N/2)

plt.figure(figsize=(5,3))
plt.plot(t_ref[int(n_ref/N/2):-int(n_ref/N/2)]*N-0.5, y_ref[int(n_ref/N/2):-int(n_ref/N/2)], '-', color='gray', linewidth=0.5)
plt.axhline(0, ls='-', color='k')
plt.axvline(0, ls='-', color='k')

# symmetry lines
plt.axvline(-0.5, ls='--', color='r')
plt.axvline(n-0.5, ls='--', color='r')

# highlight array contents
plt.axvspan(0-0.25, (n-1)+0.25, alpha=0.3, color='gray')

plt.plot(t, y1, 'bo')

for i in range(N):
    if i==n-1:
        plt.text(t[int(i+N/2)]-0.7, y1[int(i+N/2)], chr(ord("a")+i))
    elif i==N-1:
        plt.text(t[int(i+N/2)]-0.7, y1[int(i+N/2)], chr(ord("a")+(N-i-1)))
    elif i>=n:
        plt.text(t[int(i+N/2)]+0.2, y1[int(i+N/2)], chr(ord("a")+(N-i-1)))
    else:
        plt.text(t[int(i+N/2)]+0.2, y1[int(i+N/2)], chr(ord("a")+i))

# only integer xaxis ticks
plt.gca().xaxis.set_major_locator(MaxNLocator(steps=(1,10), integer=True))
plt.xlim((-4.5, 8.5))

plt.grid(True)
plt.xlabel("j")
plt.title("REDFT10 N=%d n=%d"%(N, n))
plt.tight_layout()
plt.savefig("redft10.png")

#%% REDFT01

plt.figure(figsize=(5,3))
plt.plot(t_ref*N, y_ref_2, '-', color='gray', linewidth=0.5)
plt.axhline(0, ls='-', color='k')
plt.axvline(0, ls='-', color='k')

# symmetry lines
plt.axvline(0, ls='--', color='r')
plt.axvline(n, ls='--', color='b')

# highlight array contents
plt.axvspan(0-0.25, (n-1)+0.25, alpha=0.3, color='gray')

plt.plot(t, y2, 'bo')

for i in range(N):
    if i>n:
        plt.text(t[int(i+N/2)]+0.2, y2[int(i+N/2)], "-"+chr(ord("a")+(N-i)))
    elif i<n:
        plt.text(t[int(i+N/2)]+0.2, y2[int(i+N/2)], chr(ord("a")+i))

# only integer xaxis ticks
plt.gca().xaxis.set_major_locator(MaxNLocator(steps=(1,10), integer=True))
plt.xlim((-4.5, 8.5))

plt.grid(True)
plt.xlabel("j")
plt.title("REDFT01 N=%d n=%d"%(N, n))
plt.tight_layout()
plt.savefig("redft01.png")


#%% REDFT11

plt.figure(figsize=(5,3))
plt.plot(t_ref[int(n_ref/N/2):-int(n_ref/N/2)]*N-0.5, y_ref_2[int(n_ref/N/2):-int(n_ref/N/2)], '-', color='gray', linewidth=0.5)
plt.axhline(0, ls='-', color='k')
plt.axvline(0, ls='-', color='k')

# symmetry lines
plt.axvline(-0.5, ls='--', color='r')
plt.axvline(n-0.5, ls='--', color='b')

# highlight array contents
plt.axvspan(0-0.25, (n-1)+0.25, alpha=0.3, color='gray')

plt.plot(t, y3, 'bo')

for i in range(N):
    if i==n-1:
        plt.text(t[int(i+N/2)]-0.6, y3[int(i+N/2)], chr(ord("a")+i))
    elif i>=n:
        plt.text(t[int(i+N/2)]+0.2, y3[int(i+N/2)], "-"+chr(ord("a")+(N-i-1)))
    else:
        plt.text(t[int(i+N/2)]+0.2, y3[int(i+N/2)], chr(ord("a")+i))

# only integer xaxis ticks
plt.gca().xaxis.set_major_locator(MaxNLocator(steps=(1,10), integer=True))
plt.xlim((-4.5, 8.5))

plt.grid(True)
plt.xlabel("j")
plt.title("REDFT11 N=%d n=%d"%(N, n))
plt.tight_layout()
plt.savefig("redft11.png")