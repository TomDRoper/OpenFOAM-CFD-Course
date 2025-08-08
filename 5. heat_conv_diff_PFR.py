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
Vr = math.pi*(D/2)**2 * xL          # tubing volume in m^3
Ac = math.pi*(D/2)**2               # cross-sectional area in m^2

# ===============================
# SECTION 3: Thermal and Material Properties
# ===============================
k  = 0.12                           # thermal conductivity W/(m*K)
p  = 1750                           # density (kg/m^3)
Cp = 1172                           # specific heat (J/kg/K)
a  = k/(p*Cp)                       # thermal diffusivity m^2/s

# ===============================
# SECTION 4: Discretization Parameters
# ===============================
x0  = 0                             # tubing entrance
nx  = 200                           # number of spatial grid points
dx  = xL/(nx-1)                     # spatial step-size
dt  = 0.01                          # time step
lam = a*dt/dx**2                    # diffusion coefficient for FDM

# ===============================
# SECTION 5: Heat Transfer and Boundary Conditions
# ===============================
Nu  = 3.66                          # Nusselt number (laminar, constant wall T)
h   = Nu*k/D                        # convective heat transfer coefficient (W/m^2/K)
T0  = 130 + 273.15                  # inlet temperature (K)
Tw  = 25  + 273.15                  # wall temperature (K)
reltol = 1e-9                       # tolerance for convergence

# ===============================
# SECTION 6: Flowrate Loop Parameters
# ===============================
Qnum  = 4                           # number of flowrates to model (0.5,1,2,4 mL/min)
Qfact = 2                           # factor between successive flowrates
Qmlm_start = 0.5                    # starting flowrate in mL/min
x = np.linspace(x0, xL, nx)         # spatial grid

plt.figure()                        # one figure for all curves

# ===============================
# SECTION 7: Main Loop for Multiple Flowrates
# ===============================
Qmlm = Qmlm_start
for m in range(Qnum):
    # --- convert current flowrate to SI for simulation ---
    Q = (Qmlm * 1e-6) / 60          # volumetric flowrate in m^3/s
    u = Q / Ac                       # average velocity m/s

    # CFL check (advection)
    cfl = u*dt/dx
    if cfl >= 1:
        print(f'CFL condition not met (CFL={cfl:.3f}) at Q={Qmlm:.1f} mL/min')

    # ===============================
    # SECTION 8: Array Initialization
    # ===============================
    Ttol = np.zeros(nx)             # tolerance check array
    T    = np.ones(nx)*Tw           # solution array / initial condition
    Tn   = np.ones(nx)*Tw           # temp array / initial condition

    termcond = abs(np.linalg.norm(Ttol) - np.linalg.norm(Tn))

    # ===============================
    # SECTION 9: Iterative Solver Loop
    # ===============================
    while termcond >= reltol:
        termcond = abs(np.linalg.norm(Ttol) - np.linalg.norm(Tn))
        # boundary conditions
        T[0]     = T0               # Dirichlet at inlet
        T[-1]    = T[-2]            # Neumann at outlet
        Tn = T.copy()

        # advection + diffusion + wall heat loss
        T[1:-1] = (Tn[1:-1]
                   - u*(dt/dx)*(Tn[1:-1] - Tn[:-2])
                   + lam*(Tn[2:] - 2*Tn[1:-1] + Tn[:-2])
                   - h*D*math.pi*(Tn[1:-1]-Tw)*dt/(p*Cp)*xL/Vr)

        Ttol = T.copy()

    # ===============================
    # SECTION 10: Plot Results
    # ===============================
    plt.plot(x*100, T-273.15, label=f'Q = {Qmlm:.1f} mL/min')

    # update flowrate for next iteration
    Qmlm *= Qfact

# ===============================
# SECTION 11: Display All Plots
# ===============================
plt.xlabel('Tubing length (cm)')
plt.ylabel('Fluid temperature (°C)')
plt.legend()
plt.tight_layout()
plt.show()
