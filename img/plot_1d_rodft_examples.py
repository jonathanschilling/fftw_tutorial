#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 18:59:42 2021

@author: jonathan
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

omega = 2.0*np.pi

# eval cosine at n_ref points for background
n_ref = 1024

# logical DFT size
N = 8

t_ref = np.arange(-n_ref/2*(1+1/N), n_ref)/n_ref+1/N
y_ref = np.sin(omega*t_ref)
y_ref_1 = np.sin(omega*(t_ref+0.5))
y_ref_2 = np.sin(omega*(t_ref/2))

t = np.arange(-N/2, N)
y0 = np.sin(omega*(t+1)/N)
y1 = np.sin(omega*(t+0.5)/N)
y2 = np.sin(omega*(t+1)/N/2)
y3 = np.sin(omega*(t+0.5)/N/2)

#%% RODFT00 
n = int(N/2)-1

plt.figure(figsize=(5,3))

plt.plot(t_ref[int(n_ref/N/2):]*N-1, y_ref[int(n_ref/N/2):], '-', color='gray', linewidth=0.5)
plt.axhline(0, ls='-', color='k')
plt.axvline(0, ls='-', color='k')

# symmetry lines
plt.axvline(-1, ls='--', color='b')
plt.axvline(n, ls='--', color='b')

# highlight array contents
plt.axvspan(0-0.25, (n-1)+0.25, alpha=0.3, color='gray')

plt.plot(t, y0, 'bo')

for i in range(n):
    plt.text(t[int(i+N/2)]+0.2, y0[int(i+N/2)], chr(ord("a")+i))
    plt.text(t[int(N-i+2)]+0.2, y0[int(N-i+2)]+0.05, "-"+chr(ord("a")+i))
    

# only integer xaxis ticks
plt.gca().xaxis.set_major_locator(MaxNLocator(steps=(1,10), integer=True))
plt.xlim((-4.5, 8.5))

plt.grid(True)
plt.xlabel("j")
plt.title("RODFT00 N=%d n=%d"%(N, n))
plt.tight_layout()
plt.savefig("rodft00.png")

#%% RODFT10
n = int(N/2)

plt.figure(figsize=(5,3))
plt.plot(t_ref[:-int(n_ref/N/2)]*N-0.5, y_ref[:-int(n_ref/N/2)], '-', color='gray', linewidth=0.5)
plt.axhline(0, ls='-', color='k')
plt.axvline(0, ls='-', color='k')

# symmetry lines
plt.axvline(-0.5, ls='--', color='b')
plt.axvline(n-0.5, ls='--', color='b')

# highlight array contents
plt.axvspan(0-0.25, (n-1)+0.25, alpha=0.3, color='gray')

plt.plot(t, y1, 'bo')

for i in range(N):
    if i==n-1:
        plt.text(t[int(i+N/2)]-0.7, y1[int(i+N/2)], chr(ord("a")+i))
    elif i==N-1:
        plt.text(t[int(i+N/2)]-0.7, y1[int(i+N/2)], "-"+chr(ord("a")+(N-i-1)))
    elif i>=n:
        plt.text(t[int(i+N/2)]+0.2, y1[int(i+N/2)], "-"+chr(ord("a")+(N-i-1)))
    else:
        plt.text(t[int(i+N/2)]+0.2, y1[int(i+N/2)], chr(ord("a")+i))

# only integer xaxis ticks
plt.gca().xaxis.set_major_locator(MaxNLocator(steps=(1,10), integer=True))
plt.xlim((-4.5, 8.5))

plt.grid(True)
plt.xlabel("j")
plt.title("RODFT10 N=%d n=%d"%(N, n))
plt.tight_layout()
plt.savefig("rodft10.png")

#%% RODFT01

plt.figure(figsize=(5,3))
plt.plot(t_ref[int(n_ref/N/2):]*N-1, y_ref_2[int(n_ref/N/2):], '-', color='gray', linewidth=0.5)
plt.axhline(0, ls='-', color='k')
plt.axvline(0, ls='-', color='k')

# symmetry lines
plt.axvline(-1, ls='--', color='b')
plt.axvline(n-1, ls='--', color='r')

# highlight array contents
plt.axvspan(0-0.25, (n-1)+0.25, alpha=0.3, color='gray')

plt.plot(t, y2, 'bo')

for i in range(N):
    if i>=n and i<N-1:
        plt.text(t[int(i+N/2)]+0.2, y2[int(i+N/2)], chr(ord("a")+(N-i-2)))
    elif i<n:
        plt.text(t[int(i+N/2)]+0.2, y2[int(i+N/2)], chr(ord("a")+i))

# only integer xaxis ticks
plt.gca().xaxis.set_major_locator(MaxNLocator(steps=(1,10), integer=True))
plt.xlim((-4.5, 8.5))

plt.grid(True)
plt.xlabel("j")
plt.title("RODFT01 N=%d n=%d"%(N, n))
plt.tight_layout()
plt.savefig("rodft01.png")


#%% RODFT11

plt.figure(figsize=(5,3))
plt.plot(t_ref[:-int(n_ref/N/2)]*N-0.5, y_ref_2[:-int(n_ref/N/2)], '-', color='gray', linewidth=0.5)
plt.axhline(0, ls='-', color='k')
plt.axvline(0, ls='-', color='k')

# symmetry lines
plt.axvline(-0.5, ls='--', color='b')
plt.axvline(n-0.5, ls='--', color='r')

# highlight array contents
plt.axvspan(0-0.25, (n-1)+0.25, alpha=0.3, color='gray')

plt.plot(t, y3, 'bo')

for i in range(N):
    if i==n-1:
        plt.text(t[int(i+N/2)]-0.6, y3[int(i+N/2)], chr(ord("a")+i))
    elif i>=n:
        plt.text(t[int(i+N/2)]+0.2, y3[int(i+N/2)], chr(ord("a")+(N-i-1)))
    else:
        plt.text(t[int(i+N/2)]+0.2, y3[int(i+N/2)], chr(ord("a")+i))

# only integer xaxis ticks
plt.gca().xaxis.set_major_locator(MaxNLocator(steps=(1,10), integer=True))
plt.xlim((-4.5, 8.5))

plt.grid(True)
plt.xlabel("j")
plt.title("RODFT11 N=%d n=%d"%(N, n))
plt.tight_layout()
plt.savefig("rodft11.png")