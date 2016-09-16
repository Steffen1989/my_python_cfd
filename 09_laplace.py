# Import libraries
###############################
import numpy # numerics
from matplotlib import pyplot # plotting
from matplotlib import cm # colormap
from mpl_toolkits.mplot3d import Axes3D # 3d plot

######################################################################
# FUNCTION DEFINITIONS                                               #
######################################################################
def laplace2d(p, y, dx, dy, l1norm_target):
    l1norm = 1
    pp = numpy.empty_like(p) # p at previous step

    while l1norm > l1norm_target:
        pp = p.copy()
        p[1:-1,1:-1] = (dy**2*(pp[2:,1:-1]+pp[0:-2,1:-1])+\
                        dx**2*(pp[1:-1,2:]+pp[1:-1,0:-2]))/(2*(dx**2+dy**2))

        p[0,:] = 0  # p = 0 @ x = 0
        p[-1,:] = y # p = y @ x = 2
        p[:,0] = p[:,1] # dp/dy = 0 @ y = 0
        p[:,-1] = p[:,-2]   # dp/dy = 0 @ y = 1
        l1norm = (numpy.sum(numpy.abs(p[:])-numpy.abs(pp[:])))\
                 /numpy.sum(numpy.abs(pp[:]))

    return p
        

def plot2d(x, y, p, title):
    fig = pyplot.figure(figsize=(11, 7), dpi=100)
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
nx = 31 # Remark: Must not be too large --> small dx --> instability?
ny = 31 
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
p[0,:] = 0  # p = 0 @ x = 0
p[-1,:] = y # p = y @ x = 2
p[:,0] = p[:,1] # dp/dy = 0 @ y = 0
p[:,-1] = p[:,-2]   # dp/dy = 0 @ y = 1

######################################################################
# SOME FUN                                                           #
######################################################################
plot2d(x, y, p, 'initial')
p = laplace2d(p, y, dx, dy, 1e-4)
plot2d(x, y, p, 'L1=1e-4')
pyplot.show()
