# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 21:28:17 2022

@author: colin
"""

# =============================================================================
# 
# Explicit Finite Difference Method Code to Solve the 1D Linear Transport Equation
# Adapted by: Cameron Armstrong (2019)
# Source: Lorena Barba, 12 Steps to NS in Python
# Institution: Virginia Commonwealth University
# 
# =============================================================================

# Required Modules
import numpy as np
from matplotlib import pyplot as plt
import time

xl = 2                                      # x length
nx = 600                                    # number of grid points
x = np.linspace(0,xl,nx)                    # x grid evenly spaced between 0 and 2
dx = xl/(nx-1)                              # x stepsize
nt = 200                                    # number of timesteps
dt = 0.0025                                 # time stepsize
c = 1                                       # wave speed
g = .01                                     # gaussian variance parameter (peak width)
theta = x/(0.5*xl)                          # gaussian mean parameter (peak position)
cfl = round(c*dt/dx,2)                      # cfl condition 2 decimal places

# Fun little CFL condition check and print report
if cfl >= 1:
    print('Hold your horses! The CFL is %s, which is over 1' %(cfl))
else:
    print('CFL = %s' %(cfl))

# Array Initialization
u = np.ones(nx)                             # initializing solution array
un = np.ones(nx)                            # initializing temporary solution array
u = (1/(2*np.sqrt(np.pi*(g))))*np.exp(-(1-theta)**2/(4*g)) # initial condition (IC) as a gaussian
ui = u.copy()
plt.plot(x,u); # plots IC

# BDS/Upwind with inner for-loop with example on process timing
start = time.process_time()
for n in range(nt):
    un = u.copy()
    for i in range(1,nx-1):
        u[i] = un[i] - c*dt/(dx)*(un[i]-un[i-1])
        # periodic BC's
        u[0] = u[nx-2] 
        u[nx-1] = u[1]

end = time.process_time()
print(end-start)

# BDS/Upwind with vectorization
# start = time.process_time()
# for n in range(nt):
#     un = u.copy()
#     u[1:-1] = un[1:-1] - c*dt/(dx)*(un[1:-1]-un[:-2])
#     # periodic BC's
#     u[0] = u[nx-2]
#     u[nx-1] = u[1]

# end = time.process_time()
# print(end-start)

# # CDS with inner for-loop
#for n in range(nt):
#    un = u.copy()
#    for i in range(1,nx-1):
#        u[i] = un[i] - c*dt/(2*dx)*(un[i+1]-un[i-1])
#        # periodic BC's
#        u[0] = u[nx-2]
#        u[nx-1] = u[1]

# # CDS with vectorization
#for n in range(nt):
#    un = u.copy()
#    u[1:-1] = un[1:-1] - c*dt/(2*dx)*(un[2:]-un[:-2])
#    # periodic BC's
#    u[0] = u[nx-2]
#    u[nx-1] = u[1]

plt.plot(x,u);
plt.show()
plt.plot(u,x);
plt.show()