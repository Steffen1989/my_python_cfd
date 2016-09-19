# Import libraries
###############################
import numpy # numerics
from matplotlib import pyplot # plotting
from matplotlib import cm # colormap
from mpl_toolkits.mplot3d import Axes3D # 3d plot

######################################################################
# FUNCTION DEFINITIONS                                               #
######################################################################
def poisson2d(p, pp, b, dx, dy, steps):
    for i in range(steps):
        p[1:-1,1:-1] = ((pp[2:,1:-1]+pp[0:-2,1:-1])*dy**2+\
                        (pp[1:-1,2:]+pp[1:-1,0:-2])*dx**2-\
                        b[1:-1,1:-1]*dx**2*dy**2)/\
                        (2*(dx**2+dy**2))

        # BCs
        p[0,:] = 0   # @ x=0
        p[-1,:] = 0  # @ x=2
        p[:,0] = 0   # @ y=0
        p[:,-1] = 0  # @ y=1

        # Update pp
        pp = p

    return p


def plot2d(x, y, p, title):
    fig = pyplot.figure(figsize=(8, 5), dpi=100)
    ax = fig.gca(projection='3d')
    X, Y = numpy.meshgrid(x, y, indexing='ij')
    surf = ax.plot_surface(X, Y, p, rstride=1, cstride=1, cmap=cm.coolwarm, \
            linewidth=0, antialiased=False)
    ax.set_xlim(0,2)
    ax.set_ylim(0,1)
    ax.view_init(30,225)
    pyplot.title(title)

######################################################################
# USER INPUT                                                         #
######################################################################

# Physical parameters
##############################

# Geometry and spatial discretization
##############################
nx = 50 # Remark: Must not be too large --> small dx --> instability?
ny = 50 
xmin = 0
xmax = 2
ymin = 0
ymax = 1
dx = (xmax-xmin)/(nx-1)
dy = (ymax-ymin)/(ny-1)
x = numpy.linspace(xmin, xmax, nx)
y = numpy.linspace(ymin, ymax, ny)

# Temporal discretization
##############################
#nt = 201 # Remark: nt had to be increased compared to 07 again

# Initial condition
##############################
# Initialize dependent variables
p = numpy.zeros((nx, ny))
pp = p.copy() # p at previous step
b = numpy.zeros((nx, ny))
b[nx/4,ny/4] = 100
b[3/4*nx,3/4*ny] = -100

######################################################################
# SOME FUN                                                           #
######################################################################
plot2d(x, y, p, 'p initial')
plot2d(x, y, b, 'b initial')

p = poisson2d(p, pp, b, dx, dy, 100)
plot2d(x, y, p, 'after 100')
p = poisson2d(p, pp, b, dx, dy, 900)
plot2d(x, y, p, 'after 1000')

pyplot.show()
