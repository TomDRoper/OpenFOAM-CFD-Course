# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 10:12:24 2022

@author: colin
"""

# =============================================================================
# 
# Explicit Finite Difference Method Code to Solve the 1D Diffusion/Heat Equation
# Adapted by: Cameron Armstrong (2019)
# Source: Lorena Barba, 12 Steps to NS in Python
# Institution: Virginia Commonwealth University
# 
# =============================================================================

# Required Modules
import numpy as np
from matplotlib import pyplot as plt

# Grid and Numerical Parameters
xl = 2                                      # x length
nx = 100                                    # number of grid points
x = np.linspace(0,2,nx)                     # x grid 
dx = xl/(nx-1)                              # x stepsize
nt = 500                                    # number of time steps
nu = 0.01                                  # diffusivity term
sigma = 0.4                                 # parameter for numerical stability
dt = sigma*dx**2/nu                         # time step defined using nu and sigma

# Array Initialization and Plot Conditions
u = np.zeros(nx)                            # initializing solution array
un = np.zeros(nx)                           # initializing temporary solution array
u[int(.8/dx):int(1.2/dx+1)]=1               # initial condition 
plt.plot(x,u)                               # plotting initial condition (IC)
plotcond = np.zeros(nt)                     # plot condition, determines when to plot solution
p = 1                                       # plot condition counter
f = 50                                      # plot frequency

for n in range(nt): # main time-marching loop
    un = u.copy() # update temporary solution array as IC
    # u[0]=u[1] #neumann BC
    # u[-1]=u[-2] #neumann BC
    u[0]=0 #dirichlet BC
    u[nx-1]=0 #dirichlet BC
    u[1:-1] = un[1:-1] + nu*dt/dx**2*(un[2:]-2*un[1:-1]+un[:-2]) # vectorized FDM solution using forward Euler and CDS
    plotcond[n] = np.round((n/f/p),1) # logging plot condition
    plotcond[n] = plotcond[n].astype(int) # converting to an integer
    if plotcond[n] == 1: #checking plot condition
        plt.figure(1)
        plt.plot(x,u,label='t = %s s' %(np.round(n*dt,1)))
        plt.legend()
        #plt.ylim(0,1.2)
        p += 1 # updating plot condition counter