# =============================================================================
# 
# Explicit Finite Difference Method Code
# Solves the 2D Temperature Convection-Diffusion Equation
# Assumes Tubular Plug-Flow-Reactor in Laminar Regime 
# Assumes hagen poiseuille velocity profile
# Heat Source-Sink Included Uses Laminar Nusselt Correlation for "h"
# Written by: Cameron Armstrong (2020)
# Institution: Virginia Commonwealth University
# 
# =============================================================================

# Required Modules
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D 
import math
from array import array

D = 0.0015875                                   # tubing diameter in m
xl = 30/100                                     # tubing length in m & x range
yl = D                                          # tubing diameter & y range
nx = 300                                        # x grid points
ny = 50                                         # y grid points
dx = xl/(nx-1)                                  # x stepsize
dy = yl/(ny-1)                                  # y stepsize
k= .12                                          # thermal conductvity W/(m*K)
p = 1750                                        # density (kg/m3)
Cp = 1172                                       # specifc heat (J/kg/K)
a = k/(p*Cp)                                    # thermal diffusivity (m2/s)
sigma = .001                                    # time step factor
dt = sigma * dx * dy / a                        # time stepsize 
Vr = math.pi*(D/2)**2*xl                        # tubing volume (m3)
Qmlm = 1                                        # volumetric flowrate (mL/min)
Q = (Qmlm*10**-6)/60                            # volumetric flowrate (m3/s)
Ac = math.pi*(D/2)**2                           # cross-sectional area (m2)
lamx = a*dt/dx**2                               # lumped coefficient
lamy = a*dt/dy**2                               # lumped coefficient
Nu = 3.66                                       # nusselt laminar flow in tube
h = Nu*k/D                                      # convective heat transfer coefficient (W/m2/K)
T0 = 130+273.15                                 # stream inlet temperature (degK)
Tw = 25+273.15                                  # wall temperature (degK)
reltol = 1e-8                                   # tolerance for convergence  

# grid formation
x = np.linspace(0, xl, nx) 
y = np.linspace(0, yl, ny) 
X, Y = np.meshgrid(x, y) 

# hagen poiseuille velocity field generation
uAvg = Q/Ac                                     # average velocity (m/s)
uMax = 2*uAvg                                   # max velocity (m/s)
u = np.zeros(ny)                                # array initilization
u[:] = np.linspace(-(D/2),(D/2),ny)             # array intialization
u[:] = uMax*(1-(u[:]/(D/2))**2)                 # hagan-poiselle profile
u[0]=u[-1]=0                                    # no slip BC
u = np.array([u,]*nx)                           # velocity field
u = u.T                                         # transpose/align field
maxCFL = np.max(u*dt/dx)                        # CFL condition calc.
print('The max CFL is %s'%(maxCFL))

# main function loop
def lets_get_tubular(): 
    # array initialization
    Ttol = np.zeros((ny,nx)) 
    T = np.ones((ny, nx))*Tw  
    Tn = np.ones((ny, nx))*Tw 
    # initialize termination condition
    # compares norms of current and previous solution arrays
    termcond = (np.abs((np.linalg.norm(Ttol)-np.linalg.norm(Tn))))/np.linalg.norm(Tn)
    stepcount = 1 # step counter
    while termcond >= reltol: 
        termcond = np.abs((np.linalg.norm(Ttol)-np.linalg.norm(Tn)))/np.linalg.norm(Tn)
        Tn = T.copy()
        # FDM vectorized solution using explicit euler and CDS
        T[1:-1, 1:-1] = (Tn[1:-1,1:-1] - (u[1:-1,1:-1]*(dt/(2*dx))*(Tn[1:-1,2:]  \
                     -Tn[1:-1,:-2]))                                             \
                     + lamx *(Tn[1:-1, 2:] - 2 * Tn[1:-1, 1:-1] + Tn[1:-1, :-2]) \
                     + lamy* (Tn[2:,1:-1] - 2 * Tn[1:-1, 1:-1] + Tn[:-2, 1:-1])) \
                     - h*D*math.pi*(Tn[1:-1,1:-1]-Tw)*dt/p/Cp*xl/Vr
        # BCs
        T[0, :] = Tw # tubing wall temp dirichlet BC
        T[-1, :] = Tw # tubing wall temp dirichlet BC
        T[:, 0] = T0 # inlet flow temp dirichlet BC
        T[:, -1] = T[:,-2] # outlet flow temp neumann BC
        Ttol=T.copy() # update solution
        stepcount += 1 # update counter
    
#    fig = plt.figure()
#    ax = fig.gca(projection='3d')
#    surf = ax.plot_surface(X, Y, T[:], rstride=1, cstride=1, cmap=cm.viridis,
#        linewidth=0, antialiased=True)
#    ax.set_xlabel('$x$')
#    ax.set_ylabel('$y$');
    T[:]=T[:]-273.15                            # converts back to degC
    
    # generates plots
    # top plot is 2D filled contour plot 
    # bottom plot is centerline and near-wall line data points
    fig1 = plt.subplot(211)
#    ax = fig1.gca()
#    plt.imshow(T[:])
    cont = plt.contourf(X,Y,T[:],50)
    ax = plt.gca()
    ax.axis('scaled')
    ax.axes.get_yaxis().set_visible(False)
    plt.xlim(0,.05)
    plt.xlabel('Tubing Length (m)')
    cbar = plt.colorbar(cont)
    cbar.ax.set_ylabel('Temperature (degC)')
    
    centerline = ny/2
    wallline = ny-5
    centerline = int(centerline)
    wallline = int(wallline)
    centerT = T[centerline,:]
    wallT = T[wallline,:]
    fig2 = plt.subplot(212)
    plt.plot(x, centerT,label='center')
    plt.plot(x,wallT,label='wall')
    plt.legend()
    plt.ylabel('Temperature (degC)')
    plt.xlabel('Tubing Length (m)')
    
    plt.show()
    print('Stepcount = %s' %(stepcount))
    
if __name__ == "__main__":
    lets_get_tubular()
