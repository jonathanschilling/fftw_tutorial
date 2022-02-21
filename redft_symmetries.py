#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 13:44:56 2022

@author: jonathan

Explore DFTs with only even-m or only odd-m harmonics
based on Shao & Johnson 2008

"""

import numpy as np
import matplotlib.pyplot as plt

def dft(x):
    n = len(x)
    y = np.zeros([n], dtype=np.complex128)
    for k in range(n):
        for j in range(n):
            y[k] += x[j] * np.exp(-2j*np.pi* j * k / n)
    return y

def redft00(x):
    n = len(x)
    y = np.zeros([n])
    for k in range(n):
        y[k] = x[0] + (-1)**k * x[n-1]
        for j in range(1, n-1):
            y[k] += 2.0 * x[j] * np.cos(np.pi*j*k/(n-1))
    return y

def redft10(x):
    N = len(x)
    c = np.zeros([N])
    for k in range(N):
        for n in range(N):
            c[k] += 2.0 * x[n] * np.cos(np.pi/N * (n+0.5) * k)
    return c
    
N = 8

np.random.seed(0)
x = np.random.rand(N)

x_odd = []
for i in range(N):
    x_odd.append(0.0)
    x_odd.append(x[i])

x_odd_full = []
for xo in x_odd:
    x_odd_full.append(xo)
x_odd_full.append(0.0)
for xo in x_odd[:0:-1]:
    x_odd_full.append(xo)

x_odd.append(0.0)

plt.figure()

plt.subplot(3,1,1)
plt.plot(x, 'o')
plt.grid(True)
plt.title("x")

plt.subplot(3,1,2)
plt.plot(x_odd, 'o')
plt.grid(True)
plt.title("x_odd")

plt.subplot(3,1,3)
plt.plot(x_odd_full, 'o')
plt.grid(True)
plt.title("x_odd_full")

plt.tight_layout()

#%%

c = redft10(x)

# c_odd = redft00(x_odd)

c_odd_full = dft(x_odd_full)


#%%
# plt.figure()
# plt.plot(np.arange(2*N), c_odd_full[:2*N].real, 'o', label='first half')
# plt.plot(c_odd, '*', label='c_odd')
# plt.plot(1+np.arange(2*N-1), c_odd_full[4*N:2*N:-1].real, 'x', label='second half')
# plt.grid(True)
# plt.legend(loc='upper right')
# plt.title("symmetry of dft")

# plt.tight_layout()




#%%




plt.figure()

plt.plot(c, 'o', label='c')
plt.plot(c_odd_full[:N].real, 'x', label='c_odd_full')
#plt.plot(c_odd_full.real, 'x', label='c_odd_full')
plt.grid(True)
plt.legend(loc='upper right')

plt.tight_layout()



#%%

np.random.seed(0)
x2 = np.random.rand(N+1)

x2_evn = []
x2_evn.append(x2[0])
for i in range(1,N+1):
    x2_evn.append(0.0)
    x2_evn.append(x2[i])

x2_evn_full = []
for xe in x2_evn:
    x2_evn_full.append(xe)
for xe in x2_evn[-2:0:-1]:
    x2_evn_full.append(xe)


plt.figure()

plt.subplot(3,1,1)
plt.plot(x2, 'o')
plt.grid(True)
plt.title("x2")

plt.subplot(3,1,2)
plt.plot(x2_evn, 'o')
plt.grid(True)
plt.title("x2_evn")

plt.subplot(3,1,3)
plt.plot(x2_evn_full, 'o')
plt.grid(True)
plt.title("x2_evn_full")

plt.tight_layout()

#%%

c2 = redft00(x2)

#c2_evn = redft00(x2_evn)

c2_evn_full = dft(x2_evn_full)

#%%

plt.figure()

plt.plot(c2, 'o', label='c2')
plt.plot(c2_evn_full[:N+1].real, 'x', label='c2_evn_full')
plt.grid(True)
plt.legend(loc='upper right')

plt.tight_layout()



