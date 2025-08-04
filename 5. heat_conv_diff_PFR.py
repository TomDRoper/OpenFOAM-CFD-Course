# ===============================
# SECTION 1: Required Modules
# ===============================
import numpy as np
from matplotlib import pyplot as plt
import math

# ===============================
# SECTION 2: Geometry and Flow Parameters
# ===============================
xL = 30/100                         # tubing length in m
D = 0.0015875                       # tubing diameter in m
Vr = math.pi*(D/2)**2*xL            # tubing volume in m
Qmlm = 0.5                          # volumetric flowrate in mL/min
Q = (Qmlm*10**-6)/60                # volumetric flowrate in m3/s
Ac = math.pi*(D/2)**2               # cross-sectional area in m2
u = Q/Ac                            # average velocity m/s

# ===============================
# SECTION 3: Thermal and Material Properties
# ===============================
k = .12                             # thermal conductvity W/(m*K)
p = 1750                            # density (kg/m^3)
Cp = 1172                           # constant pressure specifc heat (J/kg/K)
a = k/(p*Cp)                        # alpha - thermal diffusivity m2/s

# ===============================
# SECTION 4: Discretization Parameters
# ===============================
x0 = 0                              # tubing entrance
nx = 200                            # number of spatial grid points
dx = xL/(nx-1)                      # spatial step-size
dt = .01                            # time step
lam = a*dt/dx**2                    # lumped coefficient

# ===============================
# SECTION 5: Heat Transfer and Boundary Conditions
# ===============================
Nu = 3.66                           # nusselt number laminar flow in tube
h = Nu*k/D                          # convective heat transfer coefficient (W/m2/K)
T0 = 130+273.15                     # stream inlet temperature in degK
Tw = 25+273.15                      # wall temperature in degK
Ttol = np.zeros(nx)                 # tolerance check array initialization
reltol = 1e-9                       # tolerance for convergence

# ===============================
# SECTION 6: Flowrate Loop Parameters
# ===============================
Qnum = 4                            # number of volumetric flowrates to model
Qfact = 2                           # factor to change each flowrate each loop
x = np.linspace(x0,xL,nx)           # spatial grid

# ===============================
# SECTION 7: Main Loop for Multiple Flowrates
# ===============================
for m in range(Qnum): # main loop for running various flowrates
    u = Q/Ac                         # recalculate average velocity per loop
    cfl = u*dt/dx                    # CFL condition
    if cfl >= 1:
        print('CFL condition not met, CFL = %s, Q = %s mL/min' %(cfl,Qmlm))
    
    # ===============================
    # SECTION 8: Array Initialization
    # ===============================
    Ttol = np.zeros(nx)              # re-initializes tolerance check loop
    T = np.ones(nx)*Tw               # array initialization / Initial condition
    Tn = np.ones(nx)*Tw              # temporary array initialization / Initial condition
    
    termcond = np.abs((np.linalg.norm(Ttol)-np.linalg.norm(Tn))) # termination condition
    
    # ===============================
    # SECTION 9: Iterative Solver Loop
    # ===============================
    while termcond >= reltol: # model loop runs until termination condition met
        termcond = np.abs((np.linalg.norm(Ttol)-np.linalg.norm(Tn))) # termination condition
        T[0] = T0                  # impose dirichlet BC (inlet)
        T[nx-1] = T[nx-2]          # impose neumann BC (outlet)
        Tn = T.copy()              # update temporary array
        
        # Vectorized FDM spatial solution (advection + diffusion + wall heat loss)
        T[1:-1] = Tn[1:-1]-(u*(dt/dx)*(Tn[1:-1]                     \
                    -Tn[:-2]))+lam*(Tn[2:]-2*Tn[1:-1]+Tn[:-2])      \
                    -h*D*math.pi*(Tn[1:-1]-Tw)*dt/p/Cp*xL/Vr
        Ttol = T.copy()            # update tolerance check array
    
    # ===============================
    # SECTION 10: Plot Results
    # ===============================
    plt.figure(1)  
    plt.plot(x*100, T-273.15, label='Q = %s mL/min' %(Qmlm))
    plt.xlabel('tubing length (cm)')
    plt.ylabel('Fluid temperature (degC)')
    plt.legend()
    
    # Update flowrate for next iteration
    Q = Q * Qfact
    Qmlm = Q * 60 / (10**-6)        # update flowrate in mL/min for legend
