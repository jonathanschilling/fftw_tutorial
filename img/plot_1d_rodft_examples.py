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

# compute the RODFT00 as defined in the FFTW reference manual
def rodft00(arr):
    n = len(arr)
    out = np.zeros([n])
    for k in range(n):
        for j in range(n):
            out[k] += 2.0 * arr[j] * np.sin(np.pi*(j+1)*(k+1)/(n+1))
    return out

# compute the RODFT10 as defined in the FFTW reference manual
def rodft10(arr):
    n = len(arr)
    out = np.zeros([n])
    for k in range(n):
        for j in range(n):
            out[k] += 2.0 * arr[j] * np.sin(np.pi*(j+0.5)*(k+1)/n)
    return out

# compute the RODFT01 as defined in the FFTW reference manual
def rodft01(arr):
    n = len(arr)
    out = np.zeros([n])
    for k in range(n):
        out[k] += (-1)**k * arr[-1]
        for j in range(n-1):
            out[k] += 2.0 * arr[j] * np.sin(np.pi*(j+1)*(k+0.5)/n)
    return out

# compute the RODFT11 as defined in the FFTW reference manual
def rodft11(arr):
    n = len(arr)
    out = np.zeros([n])
    for k in range(n):
        for j in range(n):
            out[k] += 2.0 * arr[j] * np.sin(np.pi*(j+0.5)*(k+0.5)/n)
    return out

# evaluate the RODFT00 at all points given in x
def eval_rodft00(arr, x):
    n = len(arr)
    y = np.zeros_like(x)
    for j in range(n):
        y += 2.0*arr[j]*np.sin(np.pi*(j+1)*(x+1/(n+1)))
    return y

# evaluate the RODFT10 at all points given in x
def eval_rodft10(arr, x):
    n = len(arr)
    y = np.zeros_like(x)
    for j in range(n):
        y += 2.0*arr[j]*np.sin(np.pi*(j+0.5)*(x+1/n))
    return y

# evaluate the RODFT01 at all points given in x
def eval_rodft01(arr, x):
    n = len(arr)
    y = np.zeros_like(x)
    for j in range(n-1):
        y += 2.0*arr[j]* np.sin(np.pi*(j+1)*(x+0.5/n))
    y += arr[-1]* np.sin(np.pi*n*(x+0.5/n))
    return y

# evaluate the RODFT11 at all points given in x
def eval_rodft11(arr, x):
    n = len(arr)
    y = np.zeros_like(x)
    for j in range(n):
        y += 2.0*arr[j]*np.sin(np.pi*(j+0.5)*(x+0.5/n))
    return y

#%% RODFT00 

# input size for RODFT00
n = 3

# logical size of equivalent DFT (for RODFT00)
N = 2*(n+1)

# random Fourier coefficients
rng = default_rng(seed=42)
r2r_out = rng.uniform(size=n)

# compute input of RODFT00 by RODFT00 (the inverse of RODFT00)
r2r_in = rodft00(r2r_out)/n

# extend to left for highlighting symmetry at start of array
left_in_x = np.arange(-N+n, 0)
left_in_y = np.zeros([N-n])
for i in range(n):
    # odd symmetry around -1
    left_in_y[-i-2] = -r2r_in[i]

# extend to equivalent DFT size to highlight symmetry at end of array
full_in_x = np.arange(n, N)
full_in_y = np.zeros([N-n])
for i in range(n):
    # even symmetry around n-1
    full_in_y[i+1] = -r2r_in[-i-1]

# sample at finer intervals
nFine = 1024
x = np.linspace(-N/2, N, nFine)
y = eval_rodft00(r2r_out, x/(n+1))/n

plt.figure(figsize=(5,3))

plt.plot(x, y, '-', color='gray', linewidth=0.5)
plt.axhline(0, ls='-', color='k')
plt.axvline(0, ls='-', color='k')

# symmetry lines
plt.axvline(-1, ls='--', color='b')
plt.axvline(n, ls='--', color='b')

# highlight array contents
plt.axvspan(0-0.25, (n-1)+0.25, alpha=0.3, color='gray')

# plot actual RODFT00 input
plt.plot(r2r_in, 'bo')

# label individual points
for i in range(n):
    plt.text(i+0.2, r2r_in[i]-0.13, chr(ord("a")+i))
plt.text(-1+0.2, 0.0+0.05, "0")

plt.plot(left_in_x, left_in_y, 'bo')

# plot data extending REDFT00 to full DFT
plt.plot(full_in_x, full_in_y, 'bo')

plt.text(full_in_x[0]+0.2, full_in_y[0]+0.05, "0")
for i in range(1,n+1):
    plt.text(full_in_x[i]+0.2, full_in_y[i]+0.07, "-"+chr(ord("a")+(n-i)))
plt.text(full_in_x[-1]+0.2, full_in_y[-1]+0.05, "0")

# only integer xaxis ticks at intervals of 1
plt.gca().xaxis.set_major_locator(MaxNLocator(steps=(1,10), integer=True))
plt.xlim((-4.5, 8.5))

plt.grid(True)
plt.xlabel("j")
plt.title("RODFT00 N=%d n=%d"%(N, n))
plt.tight_layout()
plt.savefig("rodft00.png")

#%% RODFT10
n = 4

# logical size of equivalent DFT (for RODFT10)
N = 2*n

# random Fourier coefficients
rng = default_rng(seed=41)
r2r_out = rng.uniform(size=n)

# compute input of RODFT10 by using RODFT01 (the inverse of RODFT10)
r2r_in = rodft01(r2r_out)/n

# extend to left for highlighting symmetry at start of array
left_in_x = np.arange(-n, 0)
left_in_y = np.zeros([n])
for i in range(n):
    # odd symmetry around -0.5
    left_in_y[i] = -r2r_in[-i-1]

# extend to equivalent DFT size for highlighting symmetry at end of array
full_in_x = np.arange(n, N)
full_in_y = np.zeros([N-n])
for i in range(N-n):
    # odd symmetry around n-0.5
    full_in_y[i] = -r2r_in[-i-1]

# sample "inputs" at finer intervals
nFine = 1024
x = np.linspace(-n, N, nFine)
y = eval_rodft01(r2r_out, x/n)/n

plt.figure(figsize=(5,3))

plt.plot(x, y, '-', color='gray', linewidth=0.5)

plt.axhline(0, ls='-', color='k')
plt.axvline(0, ls='-', color='k')

# symmetry lines
plt.axvline(-0.5, ls='--', color='b')
plt.axvline(n-0.5, ls='--', color='b')

# highlight array contents
plt.axvspan(0-0.25, (n-1)+0.25, alpha=0.3, color='gray')

plt.plot(r2r_in, 'bo')

# # label individual points
# for i in range(n-1):
#     plt.text(i+0.2, r2r_in[i]-0.03, chr(ord("a")+i))
# plt.text(n-1-0.6, r2r_in[n-1]-0.03, chr(ord("a")+n-1))

plt.plot(left_in_x, left_in_y, 'bo')

plt.plot(full_in_x, full_in_y, 'bo')

# for i in range(n):
#     plt.text(full_in_x[i]+0.2, full_in_y[i]-0.03, chr(ord("a")+(n-1-i)))

# only integer xaxis ticks
plt.gca().xaxis.set_major_locator(MaxNLocator(steps=(1,10), integer=True))
plt.xlim((-4.5, 8.5))

plt.grid(True)
plt.xlabel("j")
plt.title("RODFT10 N=%d n=%d"%(N, n))
plt.tight_layout()
plt.savefig("rodft10.png")

#%% RODFT01
n = 4

# logical size of equivalent DFT (for RODFT01)
N = 2*n

# random Fourier coefficients
rng = default_rng(seed=42)
r2r_out = rng.uniform(size=n)

# compute input of RODFT01 by using RODFT10 (the inverse of RODFT01)
r2r_in = rodft10(r2r_out)/n

# extend to left for highlighting symmetry at start of array
left_in_x = np.arange(-n, 0)
left_in_y = np.zeros([n])
for i in range(n-1):
    # odd symmetry around -1
    left_in_y[-i-2] = -r2r_in[i]

# extend to equivalent DFT size for highlighting symmetry at end of array
full_in_x = np.arange(n, N)
full_in_y = np.zeros([N-n])
for i in range(N-n-1):
    # even symmetry around n-1
    full_in_y[i] = r2r_in[-i+2]

# sample "inputs" at finer intervals
nFine = 1024
x = np.linspace(-n, N, nFine)
y = eval_rodft10(r2r_out, x/n)/n

plt.figure(figsize=(5,3))

plt.plot(x, y, '-', color='gray', linewidth=0.5)

plt.axhline(0, ls='-', color='k')
plt.axvline(0, ls='-', color='k')

# symmetry lines
plt.axvline(-1, ls='--', color='b')
plt.axvline(n-1, ls='--', color='r')

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

for i in range(n-1):
    plt.text(full_in_x[i]+0.2, full_in_y[i]-0.03, chr(ord("a")+(n-2-i)))
plt.text(full_in_x[-1]+0.2, full_in_y[-1]+0.03, "0")

# only integer xaxis ticks
plt.gca().xaxis.set_major_locator(MaxNLocator(steps=(1,10), integer=True))
plt.xlim((-4.5, 8.5))

plt.grid(True)
plt.xlabel("j")
plt.title("RODFT01 N=%d n=%d"%(N, n))
plt.tight_layout()
plt.savefig("rodft01.png")


#%% RODFT11

# input size for RODFT11
n = 4

# logical size of equivalent DFT (for RODFT11)
N = 2*n

# random Fourier coefficients
rng = default_rng(seed=40)
r2r_out = rng.uniform(size=n)

# compute input of RODFT11 by RODFT11 (the inverse of RODFT11)
r2r_in = rodft11(r2r_out)/n

# extend to left for highlighting symmetry at start of array
left_in_x = np.arange(-n, 0)
left_in_y = np.zeros([n])
for i in range(n):
    # odd symmetry around -0.5
    left_in_y[i] = -r2r_in[-i-1]

# extend to equivalent DFT size to highlight symmetry at end of array
full_in_x = np.arange(n, N)
full_in_y = np.zeros([N-n])
for i in range(N-n):
    # even symmetry around n-0.5
    full_in_y[i] = r2r_in[-i-1]

# sample at finer intervals
nFine = 1024
x = np.linspace(-n, N, nFine)
y = eval_rodft11(r2r_out, x/n)/n

plt.figure(figsize=(5,3))

plt.plot(x, y, '-', color='gray', linewidth=0.5)
plt.axhline(0, ls='-', color='k')
plt.axvline(0, ls='-', color='k')

# symmetry lines
plt.axvline(-0.5, ls='--', color='b')
plt.axvline(n-0.5, ls='--', color='r')

# highlight array contents
plt.axvspan(0-0.25, (n-1)+0.25, alpha=0.3, color='gray')

# plot actual RODFT11 input
plt.plot(r2r_in, 'bo')

# label individual points
for i in range(n-1):
    plt.text(i+0.2, r2r_in[i]-0.05, chr(ord("a")+i))
i=n-1
plt.text(i-0.5, r2r_in[i]-0.03, chr(ord("a")+i))

plt.plot(left_in_x, left_in_y, 'bo')

# plot data extending REDFT11 to full DFT
plt.plot(full_in_x, full_in_y, 'bo')

for i in range(N-n):
    plt.text(full_in_x[i]+0.2, full_in_y[i]-0.05, chr(ord("a")+(n-1-i)))

# only integer xaxis ticks at intervals of 1
plt.gca().xaxis.set_major_locator(MaxNLocator(steps=(1,10), integer=True))
plt.xlim((-4.5, 8.5))

plt.grid(True)
plt.xlabel("j")
plt.title("RODFT11 N=%d n=%d"%(N, n))
plt.tight_layout()
plt.savefig("rodft11.png")