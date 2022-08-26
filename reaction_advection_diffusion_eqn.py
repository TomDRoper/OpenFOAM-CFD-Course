# =============================================================================
# 
# Explicit Finite Difference Method Code
# Solves the 1D Isothermal Reaction-Advection-Diffusion Equation
# Assumes A + B -> C with k-forward = k1
# Assumes Tubular Plug-Flow-Reactor in Laminar Regime 
# Written by: Cameron Armstrong (2020)
# Institution: Virginia Commonwealth University
# 
# =============================================================================

# Required Modules
import numpy as np
from matplotlib import pyplot as plt
import time

# Allows for division by zero for initial termination condition
np.seterr(divide='ignore')

F_a = 600                                   # initial conc. species A (mol/m3)               
F_b = 500                                   # initial conc. species B (mol/m3)  
D1 = .005                                   # axial dispersion coefficient (m2/s)
Q1 = 8.3*10**-8                             # volumetric flowrate (m3/s)
Ac1 = 2.0*10**-6                            # cross-sectional area (m2)
u1 = Q1/Ac1                                 # average velocity (m/s)
V1 = 1.0*10**-5                             #reactor volume (m3)
xl1 = V1/Ac1                                #tubing length (m)
nx1 = 200                                   # number spatial gridpoints
dx1 = xl1/(nx1-1)                           # spatial step size (m)
dt1 = .2*dx1                                # time step size (s) defined by dx
x1 = np.linspace(0,xl1,nx1)                 # x grid
Ea1 = 53106.5                               # activation energy (J/mol)
k01 = 1175.26                               # arrhenius factor (m3/mol/s)
R = 8.314                                   # universal gas constant (J/mol/K)
T1 = 150 + 273.15                           # set-point temperature (K)
k1 = k01*np.exp(-Ea1/R/T1)                  # rate constant calc (L/mol/s)
p = 0                                       # iteration counter

# array initialization for each species (also the ICs)
# assumes ICs for reagents to be inlet concentration
A = np.ones(nx1)*F_a
An = np.ones(nx1)*F_a
B = np.ones(nx1)*F_b
Bn = np.ones(nx1)*F_b
C = np.zeros(nx1)
Cn = np.zeros(nx1)

# tolerance array initlization and assignment
tolcheck = np.ones(nx1)
reltol = 1e-7                          

# function for checking termination condition by relative tolerance
# compares norms of current and previous solution arrays
def check_yourself(tolcheck,Cn):
    out = (np.abs((np.linalg.norm(tolcheck)-np.linalg.norm(Cn))))/np.linalg.norm(Cn) 
    return out

# Vectorized FDM schemes function calls (Central and Backward/Upwind)
# main = species to be modelled
# s1,s2 = reagent species "sources"
# stoic = -1 for reagents, 1 for products
def CDS(u1,dt1,dx1,D1,k1,main,s1,s2,stoic):
    out = main[1:-1] -(u1)*(dt1/(2*dx1))*(main[2:]-main[:-2])   \
    +D1*dt1/dx1**2*(main[2:]-2*main[1:-1]+main[:-2])            \
    +stoic*k1*s1[1:-1]*s2[1:-1]*dt1   
    return out

def UDS(u1,dt1,dx1,D1,k1,main,s1,s2,stoic):
    out = main[1:-1] -(u1)*(dt1/(dx1))*(main[1:-1]-main[:-2])   \
    +D1*dt1/dx1**2*(main[2:]-2*main[1:-1]+main[:-2])            \
    +stoic*k1*s1[1:-1]*s2[1:-1]*dt1
    return out

# main loop for time iteration
# initializes the termination conditions variable first
# runs until product (species C) has converged using relative tolerance
start = time.process_time()
termcond = check_yourself(tolcheck,Cn)
while termcond >= reltol: 
        # updates termination condition each loop
        termcond = check_yourself(tolcheck,Cn)
        # updating temporary arrays to current solution
        An = A.copy()
        Bn = B.copy()
        Cn = C.copy()
        # Dirichlet BCs for Reactor Inlet
        A[0] = F_a
        B[0] = F_b
        # Neumann BCs for Reactor Outlet
        A[nx1-1] = A[nx1-2] 
        B[nx1-1] = B[nx1-2] 
        C[nx1-1] = C[nx1-2] 
        # species solvers
        A[1:-1] = UDS(u1,dt1,dx1,D1,k1,An,An,Bn,-1)
        B[1:-1] = UDS(u1,dt1,dx1,D1,k1,Bn,An,Bn,-1)
        C[1:-1] = UDS(u1,dt1,dx1,D1,k1,Cn,An,Bn,1)
        # updates tolerance check array
        tolcheck = C.copy()
        p += 1

end = time.process_time()
cputime = round(end-start,2)
# plots species concentrations against normalized reactor length (x/xl)           
plt.figure(1)
plt.plot(x1/xl1,C)
plt.plot(x1/xl1,A)
plt.plot(x1/xl1,B)
plt.xlabel('Normalized Reactor Length (-)')
plt.ylabel('Molar Concentration R1 (mol/m3)') 

# determines limiting reagent and calculates predicted yield of product
# and then prints the yield
lim_reagent = np.min([F_a,F_b])
prod_yield = round(C[nx1-1]/lim_reagent*100,2)
print('Predicted Product Yield: %s %%' %(prod_yield))
print('The model solved in %ss and required %s iterations' %(cputime,p))

