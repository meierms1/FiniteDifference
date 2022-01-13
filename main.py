#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 20:11:43 2022

@author: maycon

Use this code to compute the trajectory of the vector field.
"""


import numpy as np
from Methods import FiniteDifference as FD, SquaredMesh as SM
from matplotlib import pyplot as plt

## -------------------------------------------------------------------------
## Numeric parameters
## -------------------------------------------------------------------------
dt = 0.01

#c_vector, dx = SM().onedimension(0,1,100)
c_vector = np.random.rand(100)
dx = 0.1
c = FD().onedimension(c_vector, dt, dx, L = 0.4, alpha = 0.01, beta = 0.01, Nt = 100)
c = np.array(c)


## -------------------------------------------------------------------------
## Plots
## -------------------------------------------------------------------------

plt.figure()
plt.plot(c[0])
for i in range(10):
    plt.plot(c[i+10])
plt.xlabel('position in x')
plt.ylabel('concentration c')
plt.title('Concentration for all points over 20 time steps')

plt.figure()
for i in range(10):
    plt.plot(c[:,2*i])
plt.xlabel('time t')
plt.ylabel('concentration c')
plt.title('Concentration evolution in time for 10 points')