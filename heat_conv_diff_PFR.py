# =============================================================================
#
# Explicit FDM Code to Solve 1D Convection-Diffusion Temperature Eqn.
# Assumes Tubular Plug-Flow-Reactor in Laminar Regime 
# Heat Source-Sink Included Uses Laminar Nusselt Correlation for "h"
# Written by: Cameron Armstrong (2019)
# Institution: Virginia Commonwealth University
# 
# =============================================================================

import numpy as np
from matplotlib import pyplot as plt
import math

xL = 30/100                         # tubing length in m
D = 0.0015875                       # tubing diameter in m
Vr = math.pi*(D/2)**2*xL            # tubing volume in m
Qmlm = 0.5                          # volumetric flowrate in mL/min
Q = (Qmlm*10**-6)/60                # volumetric flowrate in m3/s
Ac = math.pi*(D/2)**2               # cross-sectional area in m2
u = Q/Ac                            # average velocity m/s
k= .12                              # thermal conductvity W/(m*K)
p = 1750                            # density (kg/m^3)
Cp = 1172                           # constant pressure specifc heat (J/kg/K)
a = k/(p*Cp)                        # alpha - thermal diffusivity m2/s
x0 = 0                              # tubing entrance
nx=200                              # number of patial grid points
dx = xL/(nx-1)                      # spatial step-size
dt = .01                            # time step
lam = a*dt/dx**2                    # lumped coefficient
Nu = 3.66                           # nusselt number laminar flow in tube
h = Nu*k/D                          # convective heat transfer coefficient (W/m2/K)
T0 = 130+273.15                     # stream inlet temperature in degK
Tw = 25+273.15                      # wall temperature in degK
Ttol = np.zeros(nx)                 # Tolerance check array initialization
reltol = 1e-9                       # tolerance for convergence
Qnum = 4                            # number of volumetric flowrates to model
Qfact = 2                           # factor to change each flowrate each loop
x = np.linspace(x0,xL,nx)           # spatial grid

for m in range(Qnum): # main loop for running various flowrates
    u=Q/Ac # recalculates average velocity per loop
    cfl = u*dt/dx # CFL condition
    if cfl >= 1:
        print('CFL condition not met, CFL = %s, Q = %s mL/min' %(cfl,Qmlm))
    Ttol = np.zeros(nx) # re-initializes tolerance check loop
    T = np.ones(nx)*Tw # array initialization / Initial condition
    Tn = np.ones(nx)*Tw # temporary array initialization / Initial condition
    termcond = np.abs((np.linalg.norm(Ttol)-np.linalg.norm(Tn))) # termination condition
    while termcond >= reltol: # model loop runs until termination condition met
        termcond = np.abs((np.linalg.norm(Ttol)-np.linalg.norm(Tn))) # termination condition
        T[0]=T0 # impose dirichlet BC
        T[nx-1]=T[nx-2] # impose neumann BC
        Tn = T.copy() # update temporary array
        # next line is vectorized FDM spatial solution
        T[1:-1] = Tn[1:-1]-(u*(dt/dx)*(Tn[1:-1]                     \
                    -Tn[:-2]))+lam*(Tn[2:]-2*Tn[1:-1]+Tn[:-2])      \
                    -h*D*math.pi*(Tn[1:-1]-Tw)*dt/p/Cp*xL/Vr
        Ttol = T.copy() # update tolerance check array
    # plotting each loop result on same figure as part of main loop
    plt.figure(1)  
    plt.plot(x*100,T-273.15,label='Q = %s mL/min' %(Qmlm))
    plt.xlabel('tubing length (cm)')
    plt.ylabel('Fluid temperature degC')
    plt.legend()
    Q = Q*Qfact # update flowrate by constant factor for next loop
    Qmlm = Q*60/(10**-6) # update florwate in mL/min for legend