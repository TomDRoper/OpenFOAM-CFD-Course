# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 14:38:59 2022

@author: colin
"""

# =============================================================================
# 
# Explicit Finite Difference Method Code
# Solves the 2D Diffusion Equation
# Adapted by: Cameron Armstrong (2020)
# Source: Lorena Barba, 12 Steps to NS in Python
# Institution: Virginia Commonwealth University
#
# =============================================================================

# Required Modules
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D 

nx = 100                                        # number of x grid-points
ny = 100                                        # number of y grid-points
nt = 17                                         # number of time steps
nu = .05                                        # diffusivity term (m2/s)
xl = 2                                          # x-length
yl = 2                                          # y-length
dx = xl/(nx - 1)                                # x stepsize
dy = yl/(ny - 1)                                # y stepsize
sigma = .25                                     # numerical stability parameter
dt = sigma*dx*dy/nu                             # time stepsize 
T0 = 2                                          # high value/hot temp
Tw = 1                                          # low value/cool temp

# grid formation
x = np.linspace(0, xl, nx)
y = np.linspace(0, yl, ny)
X, Y = np.meshgrid(x, y)

# array initialization
u = np.ones((ny, nx))  
un = np.ones((ny, nx))
u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = T0  

# plotting IC as 3D bump
fig1 = plt.figure()
ax1 = fig1.gca(projection='3d')
surf = ax1.plot_surface(X, Y, u, rstride=1, cstride=1, cmap=cm.viridis,
        linewidth=0, antialiased=False)


ax1.set_xlim(0, xl)
ax1.set_ylim(0, yl)
ax1.set_zlim(1, 2.5)

ax1.set_xlabel('$x$')
ax1.set_ylabel('$y$');

# plotting IC as 2D
fig2  = plt.figure(figsize=(6,6), dpi=80)
ax2 = fig2.gca()
disp = plt.imshow(u[:])
ax2.set_xlabel('$x$')
ax2.set_ylabel('$y$');

# main function loop
# vectorized CDS 
def diffuse(nt):
    u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = T0    
    for n in range(nt + 1): 
        un = u.copy()
        u[1:-1, 1:-1] = (un[1:-1,1:-1] + 
                        nu * dt / dx**2 * 
                        (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                        nu * dt / dy**2 * 
                        (un[2:,1: -1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1]))
        # Dirichlet BCs
        u[0, :] = Tw
        u[-1, :] = Tw
        u[:, 0] = Tw
        u[:, -1] = Tw
    
    # plotting solution in 3D/2D
    fig1 = plt.figure()
    ax1 = fig1.gca(projection='3d')
    surf = ax1.plot_surface(X, Y, u[:], rstride=1, cstride=1, cmap=cm.viridis,
        linewidth=0, antialiased=True)
    ax1.set_zlim(1, 2.5)
    ax1.set_xlabel('$x$')
    ax1.set_ylabel('$y$');
    fig2  = plt.figure(figsize=(6,6), dpi=80)
    ax2 = fig2.gca()
    disp = plt.imshow(u[:])
    ax2.set_xlabel('$x$')
    ax2.set_ylabel('$y$');