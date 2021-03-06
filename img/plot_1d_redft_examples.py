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

# compute the REDFT10 as defined in the FFTW reference manual
def redft10(arr):
    n = len(arr)
    out = np.zeros([n])
    for k in range(n):
        for j in range(n):
            out[k] += 2.0 * arr[j] * np.cos(np.pi*(j+0.5)*k/n)
    return out

# compute the REDFT01 as defined in the FFTW reference manual
def redft01(arr):
    n = len(arr)
    out = np.zeros([n])
    for k in range(n):
        out[k] += arr[0]
        for j in range(1, n):
            out[k] += 2.0 * arr[j] * np.cos(np.pi*j*(k+0.5)/n)
    return out

# compute the REDFT11 as defined in the FFTW reference manual
def redft11(arr):
    n = len(arr)
    out = np.zeros([n])
    for k in range(n):
        for j in range(n):
            out[k] += 2.0 * arr[j] * np.cos(np.pi*(j+0.5)*(k+0.5)/n)
    return out

# evaluate the REDFT00 at all points given in x
def eval_redft00(arr, x):
    n = len(arr)
    y = np.zeros_like(x)
    y += arr[0]
    for j in range(1,n-1):
        y += 2.0*arr[j]*np.cos(np.pi*j*x)
    y += arr[-1]*np.cos(np.pi*(n-1)*x)
    return y

# evaluate the REDFT10 at all points given in x
def eval_redft10(arr, x):
    n = len(arr)
    y = np.zeros_like(x)
    for j in range(n):
        y += 2.0*arr[j]*np.cos(np.pi*(j+0.5)*x)
    return y

# evaluate the REDFT01 at all points given in x
def eval_redft01(arr, x):
    n = len(arr)
    y = np.zeros_like(x)
    y += arr[0]
    for j in range(1, n):
        y += 2.0*arr[j]*np.cos(np.pi*j*(x+0.5/n))
    return y

# evaluate the REDFT11 at all points given in x
def eval_redft11(arr, x):
    n = len(arr)
    y = np.zeros_like(x)
    for j in range(n):
        y += 2.0*arr[j]*np.cos(np.pi*(j+0.5)*(x+0.5/n))
    return y

#%% REDFT00 

# input size for REDFT00
n = 5

# logical size of equivalent DFT (for REDFT00)
N = 2*(n-1)

# random Fourier coefficients
rng = default_rng(seed=42)
r2r_out = rng.uniform(size=n)

# compute input of REDFT00 by REDFT00 (the inverse of REDFT00)
r2r_in = redft00(r2r_out)/n

# extend to left for highlighting symmetry at start of array
left_in_x = np.arange(-n, 0)
left_in_y = np.zeros([n])
for i in range(n):
    # even symmetry around 0
    left_in_y[i] = r2r_in[-i]

# extend to equivalent DFT size to highlight symmetry at end of array
full_in_x = np.arange(n, N)
full_in_y = np.zeros([N-n])
for i in range(N-n):
    # even symmetry around n-1
    full_in_y[i] = r2r_in[-i-2]

# sample at finer intervals
nFine = 1024
x = np.linspace(-n+1, N, nFine)
y = eval_redft00(r2r_out, x/(n-1))/n

plt.figure(figsize=(5,3))

plt.plot(x, y, '-', color='gray', linewidth=0.5)
plt.axhline(0, ls='-', color='k')
plt.axvline(0, ls='-', color='k')

# symmetry lines
plt.axvline(0, ls='--', color='r')
plt.axvline(n-1, ls='--', color='r')

# highlight array contents
plt.axvspan(0-0.25, (n-1)+0.25, alpha=0.3, color='gray')

# plot actual REDFT00 input
plt.plot(r2r_in, 'bo')

# label individual points
for i in range(n):
    plt.text(i+0.2, r2r_in[i]-0.03, chr(ord("a")+i))

plt.plot(left_in_x, left_in_y, 'bo')

# plot data extending REDFT00 to full DFT
plt.plot(full_in_x, full_in_y, 'bo')

for i in range(n-2):
    plt.text(full_in_x[i]+0.2, full_in_y[i]-0.03, chr(ord("a")+(n-2-i)))

# only integer xaxis ticks at intervals of 1
plt.gca().xaxis.set_major_locator(MaxNLocator(steps=(1,10), integer=True))
plt.xlim((-4.5, 8.5))

plt.grid(True)
plt.xlabel("j")
plt.title("REDFT00 N=%d n=%d"%(N, n))
plt.tight_layout()
plt.savefig("redft00.png")

#%% REDFT10
n = 4

# logical size of equivalent DFT (for REDFT10)
N = 2*n

# random Fourier coefficients
rng = default_rng(seed=40)
r2r_out = rng.uniform(size=n)

# compute input of REDFT10 by using REDFT01 (the inverse of REDFT10)
r2r_in = redft01(r2r_out)/n

# extend to left for highlighting symmetry at start of array
left_in_x = np.arange(-n, 0)
left_in_y = np.zeros([n])
for i in range(n):
    # even symmetry around -0.5
    left_in_y[i] = r2r_in[-i-1]

# extend to equivalent DFT size for highlighting symmetry at end of array
full_in_x = np.arange(n, N)
full_in_y = np.zeros([N-n])
for i in range(N-n):
    # even symmetry around n-0.5
    full_in_y[i] = r2r_in[-i-1]

# sample "inputs" at finer intervals
nFine = 1024
x = np.linspace(-n, N, nFine)
y = eval_redft01(r2r_out, x/n)/n

plt.figure(figsize=(5,3))

plt.plot(x, y, '-', color='gray', linewidth=0.5)

plt.axhline(0, ls='-', color='k')
plt.axvline(0, ls='-', color='k')

# symmetry lines
plt.axvline(-0.5, ls='--', color='r')
plt.axvline(n-0.5, ls='--', color='r')

# highlight array contents
plt.axvspan(0-0.25, (n-1)+0.25, alpha=0.3, color='gray')

plt.plot(r2r_in, 'bo')

# label individual points
for i in range(n-1):
    plt.text(i+0.2, r2r_in[i]-0.03, chr(ord("a")+i))
plt.text(n-1-0.6, r2r_in[n-1]-0.03, chr(ord("a")+n-1))

plt.plot(left_in_x, left_in_y, 'bo')

plt.plot(full_in_x, full_in_y, 'bo')

for i in range(n):
    plt.text(full_in_x[i]+0.2, full_in_y[i]-0.03, chr(ord("a")+(n-1-i)))

# only integer xaxis ticks
plt.gca().xaxis.set_major_locator(MaxNLocator(steps=(1,10), integer=True))
plt.xlim((-4.5, 8.5))

plt.grid(True)
plt.xlabel("j")
plt.title("REDFT10 N=%d n=%d"%(N, n))
plt.tight_layout()
plt.savefig("redft10.png")

#%% REDFT01
n = 4

# logical size of equivalent DFT (for REDFT01)
N = 2*n

# random Fourier coefficients
rng = default_rng(seed=41)
r2r_out = rng.uniform(size=n)

# compute input of REDFT01 by using REDFT10 (the inverse of REDFT01)
r2r_in = redft10(r2r_out)/n

# extend to left for highlighting symmetry at start of array
left_in_x = np.arange(-n, 0)
left_in_y = np.zeros([n])
for i in range(1,n):
    # even symmetry around 0
    left_in_y[-i] = r2r_in[i]

# extend to equivalent DFT size for highlighting symmetry at end of array
full_in_x = np.arange(n, N)
full_in_y = np.zeros([N-n])
for i in range(1,N-n):
    # odd symmetry around n
    full_in_y[i] = -r2r_in[-i]

# sample "inputs" at finer intervals
nFine = 1024
x = np.linspace(-n, N, nFine)
y = eval_redft10(r2r_out, x/n)/n

plt.figure(figsize=(5,3))

plt.plot(x, y, '-', color='gray', linewidth=0.5)

plt.axhline(0, ls='-', color='k')
plt.axvline(0, ls='-', color='k')

# symmetry lines
plt.axvline(0, ls='--', color='r')
plt.axvline(n, ls='--', color='b')

# highlight array contents
plt.axvspan(0-0.25, (n-1)+0.25, alpha=0.3, color='gray')

plt.plot(r2r_in, 'bo')

# label individual points
i=0
plt.text(i+0.2, r2r_in[i]-0.03, chr(ord("a")+i))
for i in range(1,n):
    plt.text(i-0.5, r2r_in[i]-0.03, chr(ord("a")+i))

plt.plot(left_in_x, left_in_y, 'bo')

plt.plot(full_in_x, full_in_y, 'bo')

plt.text(full_in_x[0]+0.2, full_in_y[0]+0.04, "0")
for i in range(1,n):
    plt.text(full_in_x[i]+0.2, full_in_y[i]-0.03, "-"+chr(ord("a")+(n-i)))

# only integer xaxis ticks
plt.gca().xaxis.set_major_locator(MaxNLocator(steps=(1,10), integer=True))
plt.xlim((-4.5, 8.5))

plt.grid(True)
plt.xlabel("j")
plt.title("REDFT01 N=%d n=%d"%(N, n))
plt.tight_layout()
plt.savefig("redft01.png")


#%% REDFT11

# input size for REDFT11
n = 4

# logical size of equivalent DFT (for REDFT11)
N = 2*n

# random Fourier coefficients
rng = default_rng(seed=42)
r2r_out = rng.uniform(size=n)

# compute input of REDFT11 by REDFT11 (the inverse of REDFT11)
r2r_in = redft11(r2r_out)/n

# extend to left for highlighting symmetry at start of array
left_in_x = np.arange(-n, 0)
left_in_y = np.zeros([n])
for i in range(n):
    # even symmetry around -0.5
    left_in_y[i] = r2r_in[-i-1]

# extend to equivalent DFT size to highlight symmetry at end of array
full_in_x = np.arange(n, N)
full_in_y = np.zeros([N-n])
for i in range(N-n):
    # odd symmetry around n-0.5
    full_in_y[i] = -r2r_in[-i-1]

# sample at finer intervals
nFine = 1024
x = np.linspace(-n, N, nFine)
y = eval_redft11(r2r_out, x/n)/n

plt.figure(figsize=(5,3))

plt.plot(x, y, '-', color='gray', linewidth=0.5)
plt.axhline(0, ls='-', color='k')
plt.axvline(0, ls='-', color='k')

# symmetry lines
plt.axvline(-0.5, ls='--', color='r')
plt.axvline(n-0.5, ls='--', color='b')

# highlight array contents
plt.axvspan(0-0.25, (n-1)+0.25, alpha=0.3, color='gray')

# plot actual REDFT11 input
plt.plot(r2r_in, 'bo')

# label individual points
for i in range(n-1):
    plt.text(i+0.2, r2r_in[i]-0.03, chr(ord("a")+i))
i=n-1
plt.text(i-0.5, r2r_in[i]-0.18, chr(ord("a")+i))

plt.plot(left_in_x, left_in_y, 'bo')

# plot data extending REDFT11 to full DFT
plt.plot(full_in_x, full_in_y, 'bo')

for i in range(N-n):
    plt.text(full_in_x[i]+0.2, full_in_y[i]-0.03, "-"+chr(ord("a")+(n-1-i)))

# only integer xaxis ticks at intervals of 1
plt.gca().xaxis.set_major_locator(MaxNLocator(steps=(1,10), integer=True))
plt.xlim((-4.5, 8.5))

plt.grid(True)
plt.xlabel("j")
plt.title("REDFT11 N=%d n=%d"%(N, n))
plt.tight_layout()
plt.savefig("redft11.png")