# ===========================================
# 1D Linear Advection (Upwind Scheme Example)
# ===========================================

# ============
# Imports
# ============
import numpy as np
from matplotlib import pyplot as plt
import time

# =======================================
# SECTION 1: Simulation Setup Parameters
# =======================================

xl = 2                      # Length of spatial domain (x)
nx = 600                    # Number of spatial grid points
x = np.linspace(0, xl, nx)  # Evenly spaced grid from 0 to xl
dx = xl / (nx - 1)          # Spatial step size

nt = 200                    # Number of time steps
dt = 0.0025                 # Time step size
c = 1                       # Constant wave speed

g = 0.01                    # Gaussian variance (peak width)
theta = x / (0.5 * xl)      # Gaussian mean (peak position scaling)
cfl = round(c * dt / dx, 2) # CFL number (rounded to 2 decimal places)

# =======================================
# SECTION 2: CFL Condition Check
# =======================================
# Stability condition for explicit schemes

if cfl >= 1:
    print(f'⚠️ Hold your horses! The CFL is {cfl}, which is over 1 — unstable!')
else:
    print(f'CFL = {cfl} — good to go.')

# =======================================
# SECTION 3: Initial Conditions
# =======================================

u = np.ones(nx)  # Solution array
un = np.ones(nx) # Temporary array for time stepping

# Gaussian initial condition
u = (1 / (2 * np.sqrt(np.pi * g))) * np.exp(-(1 - theta)**2 / (4 * g))
ui = u.copy()  # Save initial condition for comparison

# --- Ensure a figure exists in every environment ---
plt.figure("1D Linear Advection")  # (NEW) make sure a figure is created

# Plot the initial condition
plt.plot(x, u, label='Initial Condition')

# =======================================
# SECTION 4: Upwind (Backward Difference) Scheme with For-Loop
# =======================================

start = time.process_time()

for n in range(nt):
    un = u.copy()
    for i in range(1, nx - 1):
        u[i] = un[i] - c * dt / dx * (un[i] - un[i - 1])
    
    # Apply periodic boundary conditions
    u[0] = u[nx - 2]
    u[nx - 1] = u[1]

end = time.process_time()
print(f"Upwind (loop) execution time: {end - start:.4f} s")

# =======================================
# SECTION 5: Optional Vectorized Version (Faster)
# =======================================

# Uncomment below for faster vectorized version

# start = time.process_time()
# for n in range(nt):
#     un = u.copy()
#     u[1:-1] = un[1:-1] - c * dt / dx * (un[1:-1] - un[:-2])
#     u[0] = u[nx - 2]
#     u[nx - 1] = u[1]
# end = time.process_time()
# print(f"Upwind (vectorized) execution time: {end - start:.4f} s")

# =======================================
# SECTION 6: Central Difference Scheme
# =======================================

# # Re-initialize from Gaussian
# def initial_conditions():
#     return (1 / (2 * np.sqrt(np.pi * g))) * np.exp(-(1 - theta)**2 / (4 * g))

# u = initial_conditions()
# uold = u.copy()   # u^{t-1}
# un = u.copy()     # storage for updates

# # Bootstrap first step with Forward Euler (since CDS needs t-1 and t)
# uold = u.copy()
# u[1:-1] = uold[1:-1] - (c * dt / (2 * dx)) * (uold[2:] - uold[:-2])
# u[0] = u[-2]
# u[-1] = u[1]

# # Time-stepping with CDS
# for n in range(2, nt):
#     un[1:-1] = uold[1:-1] - (c * dt / dx) * (u[2:] - u[:-2])
#     # Periodic BCs
#     un[0] = un[-2]
#     un[-1] = un[1]
#     # Shift arrays
#     uold, u = u, un.copy()

# # Plot result
# plt.plot(x, u, '-.', label="Central Difference")



# =======================================
# SECTION 7: Final Plot
# =======================================

plt.plot(x, u, label='Final Condition')
plt.title("1D Linear Advection (Upwind Scheme)")
plt.xlabel("x")
plt.ylabel("u")
plt.grid(True)
plt.legend()

# --- Make layout tidy and force a blocking show across IDEs ---
try:
    plt.tight_layout()
except Exception:
    pass
plt.show(block=True)  # (NEW) block=True helps some IDEs reliably display the figure
