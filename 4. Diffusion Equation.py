# ============================================
# SECTION 1: Import Required Python Modules
# ============================================

import numpy as np
from matplotlib import pyplot as plt

# ============================================
# SECTION 2: Grid and Numerical Parameters
# ============================================

xl = 2                        # Total length of spatial domain (x)
nx = 100                      # Number of spatial grid points
x = np.linspace(0, xl, nx)    # 1D spatial grid from 0 to xl
dx = xl / (nx - 1)            # Spatial step size

nt = 500                      # Total number of time steps
nu = 0.01                     # Diffusivity coefficient (thermal or molecular)
sigma = 0.4                   # Stability factor for time step calculation
dt = sigma * dx**2 / nu       # Time step size (based on stability condition)

# ============================================
# SECTION 3: Initial Condition and Array Setup
# ============================================

u = np.zeros(nx)              # Initialize solution array
un = np.zeros(nx)             # Temporary array for updates

# Set initial condition: pulse between x = 0.8 and x = 1.2
u[int(0.8/dx):int(1.2/dx+1)] = 1

# --- Ensure a figure is created ---
plt.figure("1D Diffusion Example")

# Plot the initial condition
plt.plot(x, u, label='Initial Condition')

# Setup for controlling when intermediate plots appear
plotcond = np.zeros(nt)       # Array to track plotting intervals
p = 1                         # Counter for plot conditions
f = 50                        # Frequency of plotting (every f steps)

# ============================================
# SECTION 4: Time-Stepping Loop (Explicit FDM)
# ============================================

for n in range(nt):  
    un = u.copy()  # Store previous time step
    
    # === Boundary Conditions ===
    # Dirichlet Boundary Conditions (fixed to 0 at both ends)
    u[0] = 0
    u[-1] = 0
    # Neumann BC (Flux fixed at both ends)
    # u[0] = u[1] 
    # u[-1] = u[-2] 

    # === Update interior points using Explicit Finite Difference ===
    # Forward Euler in time + Central Difference in space
    u[1:-1] = un[1:-1] + nu * dt / dx**2 * (un[2:] - 2*un[1:-1] + un[:-2])

    # === Plot intermediate results every f steps ===
    plotcond[n] = np.round((n / (f * p)), 1).astype(int)
    if plotcond[n] == 1:
        plt.plot(x, u, label=f't = {np.round(n * dt, 1)} s')
        p += 1  # Update plot counter

# ============================================
# Final output will be shown with multiple time curves
# ============================================

plt.title("1D Diffusion Equation (Explicit FDM)")
plt.xlabel("x")
plt.ylabel("u")
plt.legend()
plt.grid(True)

# --- Force reliable display across Spyder, Jupyter, Colab, etc. ---
try:
    plt.tight_layout()
except Exception:
    pass
plt.show(block=True)
