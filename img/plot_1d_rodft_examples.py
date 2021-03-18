#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 18:59:42 2021

@author: jonathan
"""

import numpy as np
import matplotlib.pyplot as plt

omega = 2.0*np.pi

# eval cosine at 1000 points for background
t_ref = np.arange(1000)/1000
y_ref = np.sin(omega*t_ref)

# discrete sample at n points
n = 8

t0 = np.arange(n)/n
y0 = np.sin(omega*t0)

t1 = (np.arange(n)+0.5)/n
y1 = np.sin(omega*t1)

plt.figure(figsize=(5,3))
plt.plot(t_ref, y_ref, '-', color='gray', linewidth=0.5)
plt.axhline(0.0, ls='-', color='k')
plt.axvline(0.5, ls='--', color='k')
plt.plot(t0, y0, 'bo', label='input')
plt.plot(t0, y0, 'rx', label='output')
plt.grid(True)
plt.title("RODFT00")
plt.legend(loc=(0.525, 0.738))
plt.tight_layout()
plt.savefig("rodft00.png")

plt.figure(figsize=(5,3))
plt.plot(t_ref, y_ref, '-', color='gray', linewidth=0.5)
plt.axhline(0.0, ls='-', color='k')
plt.axvline(0.5, ls='--', color='k')
plt.plot(t1, y1, 'bo', label='input')
plt.plot(t0, y0, 'rx', label='output')
plt.grid(True)
plt.title("RODFT10")
plt.legend(loc=(0.525, 0.738))
plt.tight_layout()
plt.savefig("rodft10.png")

plt.figure(figsize=(5,3))
plt.plot(t_ref, y_ref, '-', color='gray', linewidth=0.5)
plt.axhline(0.0, ls='-', color='k')
plt.axvline(0.5, ls='--', color='k')
plt.plot(t0, y0, 'bo', label='input')
plt.plot(t1, y1, 'rx', label='output')
plt.grid(True)
plt.title("RODFT01")
plt.legend(loc=(0.525, 0.738))
plt.tight_layout()
plt.savefig("rodft01.png")

plt.figure(figsize=(5,3))
plt.plot(t_ref, y_ref, '-', color='gray', linewidth=0.5)
plt.axhline(0.0, ls='-', color='k')
plt.axvline(0.5, ls='--', color='k')
plt.plot(t1, y1, 'bo', label='input')
plt.plot(t1, y1, 'rx', label='output')
plt.grid(True)
plt.title("RODFT11")
plt.legend(loc=(0.525, 0.738))
plt.tight_layout()
plt.savefig("rodft11.png")