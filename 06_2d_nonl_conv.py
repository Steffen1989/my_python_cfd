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

# 3: Geometry and spatial discretization
##############################
nx = 81
ny = 81 
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
nt = 101
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
    surf = ax.plot_surface(X, Y, u, rstride=2, cstride=2, cmap=cm.viridis, \
            linewidth=0, antialiased=False)
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    surf = ax.plot_surface(X, Y, v, rstride=2, cstride=2, cmap=cm.viridis, \
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
        u[1:-1,1:-1] = up[1:-1,1:-1] - up[1:-1,1:-1]*dt/dx*(up[1:-1,1:-1] - up[0:-2,1:-1]) -\
                vp[1:-1,1:-1]*dt/dy*(up[1:-1,1:-1] - up[1:-1,0:-2])
        v[1:-1,1:-1] = vp[1:-1,1:-1] - up[1:-1,1:-1]*dt/dx*(vp[1:-1,1:-1] - vp[0:-2,1:-1]) -\
                vp[1:-1,1:-1]*dt/dy*(vp[1:-1,1:-1] - vp[1:-1,0:-2])
#        u[1:-1,1:-1] = up[1:-1,1:-1] - up[1:-1,1:-1]*dt/dx*(up[1:-1,1:-1] - up[1:-1,0:-2]) -\
#                vp[1:-1,1:-1]*dt/dy*(up[1:-1,1:-1] - up[0:-2,1:-1]) 
#        v[1:-1,1:-1] = vp[1:-1,1:-1] - up[1:-1,1:-1]*dt/dx*(vp[1:-1,1:-1] - vp[1:-1,0:-2]) -\
#                vp[1:-1,1:-1]*dt/dy*(vp[1:-1,1:-1] - vp[0:-2,1:-1]) 

#        # Boundary condition: periodic
#        # along y=0 but not at corner
#        u[0,1:-1] = up[0,1:-1] - up[0,1:-1]*dt/dx*(up[0,1:-1] - up[0,0:-2]) - \
#                vp[0,1:-1]*dt/dy*(up[0,1:-1] - up[-1,1:-1]) 
#        # along y=2 but not at corner
#        u[-1,1:-1] = up[-1,1:-1] - up[-1,1:-1]*dt/dx*(up[-1,1:-1] - up[-1,0:-2]) - \
#                vp[-1,1:-1]*dt/dy*(up[-1,1:-1] - up[-2,1:-1]) 
#        # along y=0 but not at corner
#        v[0,1:-1] = vp[0,1:-1] - up[0,1:-1]*dt/dx*(vp[0,1:-1] - vp[0,0:-2]) - \
#                vp[0,1:-1]*dt/dy*(vp[0,1:-1] - vp[-1,1:-1]) 
#        # along y=2 but not at corner
#        v[-1,1:-1] = vp[-1,1:-1] - up[-1,1:-1]*dt/dx*(vp[-1,1:-1] - vp[-1,0:-2]) - \
#                vp[-1,1:-1]*dt/dy*(vp[-1,1:-1] - vp[-2,1:-1]) 
#        # along x=0 but not at corner
#        u[1:-1,0] = up[1:-1,0] - up[1:-1,0]*dt/dx*(up[1:-1,0] - up[1:-1,-1]) - \
#                vp[1:-1,0]*dt/dy*(up[1:-1,0] - up[0:-2,0]) 
#        # along x=2 but not at corner
#        u[1:-1,-1] = up[1:-1,-1] - up[1:-1,1]*dt/dx*(up[1:-1,-1] - up[1:-1,-2]) - \
#                vp[1:-1,-1]*dt/dy*(up[-1,1:-1] - up[0:-2,-1]) 
#        # along x=0 but not at corner
#        v[1:-1,0] = vp[1:-1,0] - up[1:-1,0]*dt/dx*(vp[1:-1,0] - vp[1:-1,-1]) - \
#                vp[1:-1,0]*dt/dy*(vp[1:-1,0] - vp[0:-2,0]) 
#        # along x=2 but not at corner
#        v[1:-1,-1] = vp[1:-1,-1] - up[1:-1,-1]*dt/dx*(vp[1:-1,-1] - vp[1:-1,-2]) - \
#                vp[1:-1,-1]*dt/dy*(vp[-1,1:-1] - vp[0:-2,-1]) 
#        # Corners
#        # at x=y=0
#        u[0,0]= up[0,0] - up[0,0]*dt/dx*(up[0,0] - up[0,-1]) - \
#                vp[0,0]*dt/dy*(up[0,0] - up[-1,0])
#        v[0,0]= vp[0,0] - up[0,0]*dt/dx*(vp[0,0] - vp[0,-1]) - \
#                vp[0,0]*dt/dy*(vp[0,0] - vp[-1,0])
#        # at x=0, y=2
#        u[-1,0]= up[-1,0] - up[-1,0]*dt/dx*(up[-1,0] - up[-1,-1]) - \
#                vp[-1,0]*dt/dy*(up[-1,0] - up[-2,0])
#        v[-1,0]= vp[-1,0] - up[-1,0]*dt/dx*(vp[-1,0] - vp[-1,-1]) - \
#                vp[-1,0]*dt/dy*(vp[-1,0] - vp[-2,0])
#        # at x=2, y=2
#        u[-1,-1]= up[-1,-1] - up[-1,-1]*dt/dx*(up[-1,-1] - up[-1,-2]) - \
#                vp[-1,-1]*dt/dy*(up[-1,-1] - up[-2,-1])
#        v[-1,-1]= vp[-1,-1] - up[-1,-1]*dt/dx*(vp[-1,-1] - vp[-1,-2]) - \
#                vp[-1,-1]*dt/dy*(vp[-1,-1] - vp[-2,-1])
#        # at x=2, y=0
#        u[0,-1]= up[0,-1] - up[0,-1]*dt/dx*(up[0,-1] - up[0,-2]) - \
#                vp[0,-1]*dt/dy*(up[0,-1] - up[-1,-1])
#        v[0,-1]= vp[0,-1] - up[0,-1]*dt/dx*(vp[0,-1] - vp[0,-2]) - \
#                vp[0,-1]*dt/dy*(vp[0,-1] - vp[-1,-1])


# 6: Results
##############################
advance(nt)
plot_surf(x, y, u, v, 'After'+str(nt)+' timesteps')
#advance(nt)
#plot_surf(x, y, u, v, 'After'+str(2*nt)+' timesteps')
#advance(nt)
#plot_surf(x, y, u, v, 'After'+str(nt)+' timesteps')
#advance(nt)
#plot_surf(x, y, u, v, 'After'+str(2*nt)+' timesteps')


pyplot.show()
