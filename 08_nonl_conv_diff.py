# 1: Import libraries
##############################
import numpy # numerics
from matplotlib import pyplot # plotting
from matplotlib import cm # colormal
from mpl_toolkits.mplot3d import Axes3D # 3d plot

######################################################################
# USER INPUT                                                         #
######################################################################

# 2: Physical parameters
##############################
#c = 1 # not needed for nonlinear conv
nu = .1

# 3: Geometry and spatial discretization
##############################
nx = 51 # Remark: Must not be too large --> small dx --> instability?
ny = 51 
xmin = 0
xmax = 2
ymin = 0
ymax = 2
dx = (xmax-xmin)/(nx-1)
dy = (ymax-ymin)/(ny-1)
x = numpy.linspace(xmin, xmax, nx)
y = numpy.linspace(ymin, ymax, ny)

# 4: Temporal discretization
##############################
nt = 201 # Remark: nt had to be increased compared to 07 again
tmax = 0.5
dt = tmax/(nt-1)
#sigma = 0.2
#dt = sigma*dx

# 5: Initial condition
##############################
# Initialize dependent variables
u = numpy.zeros((nx, ny))
v = numpy.zeros((nx, ny))
u[:,:] = 1
up = u.copy()    # from previous time-step
v[:,:] = 1
vp = v.copy()    # from previous time-step

# Hat function
u[int(.5/dy):int(1/dy)+1, int(0.5/dx):int(1/dx)+1] = 2
v[int(.5/dy):int(1/dy)+1, int(0.5/dx):int(1/dx)+1] = 2


# check initial condition
######################################################################
# VISUALIZATION                                                      #
######################################################################
def plot_surf(x, y, u, v, title):
    fig = pyplot.figure(figsize=(10, 3), dpi=100)
    X, Y = numpy.meshgrid(x, y, indexing='ij')
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    surf = ax.plot_surface(X, Y, u, rstride=2, cstride=2, cmap=cm.coolwarm, \
            linewidth=0, antialiased=False)
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    surf = ax.plot_surface(X, Y, v, rstride=2, cstride=2, cmap=cm.coolwarm, \
            linewidth=0, antialiased=False)
    pyplot.suptitle(title)
#    pyplot.show()

plot_surf(x, y, u, v, 'Initial')

# 5: Numerics
######################################################################
# NUMBER CRUNCHING                                                   #
######################################################################
# Advance in time: Euler forward (Forward Differencing)
# discretize in space: Backward Differencing
def advance(ts):
    for n in range(0,ts):   # loop over all time-steps
        up = u.copy()
        vp = v.copy()
        
        u[1:-1,1:-1] = up[1:-1,1:-1] + dt*(-up[1:-1,1:-1]/dx*(up[1:-1,1:-1] - up[0:-2,1:-1]) \
                - vp[1:-1,1:-1]/dy*(up[1:-1,1:-1] - up[1:-1,0:-2]) \
                + nu*(1/dx**2*(up[2:,1:-1] - 2*up[1:-1,1:-1] + up[0:-2,1:-1]) \
                + 1/dy**2*(up[1:-1,2:] - 2*up[1:-1,1:-1] + up[1:-1,0:-2])))
        v[1:-1,1:-1] = vp[1:-1,1:-1] + dt*(-up[1:-1,1:-1]/dx*(vp[1:-1,1:-1] - vp[0:-2,1:-1]) - \
                vp[1:-1,1:-1]/dy*(vp[1:-1,1:-1] - vp[1:-1,0:-2]) + \
                nu*(1/dx**2*(vp[2:,1:-1] - 2*vp[1:-1,1:-1] + vp[1:-1,0:-2]) + \
                1/dy**2*(vp[1:-1,2:] - 2*vp[1:-1,1:-1] + vp[0:-2,1:-1])))
	
        # Boundary condition: Dirichlet
        u[0,:] = 1
        v[0,:] = 1
        u[-1,:] = 1
        v[-1,:] = 1
        u[:,0] = 1
        v[:,0] = 1
        u[:,-1] = 1
        v[:,-1] = 1


# 6: Results
##############################
advance(nt)
plot_surf(x, y, u, v, 'After'+str(nt)+' timesteps')


pyplot.show()
