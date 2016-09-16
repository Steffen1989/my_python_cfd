import numpy
from matplotlib import pyplot
from matplotlib import cm

## Variable declaration
# parameters
nu = .1

# time 
t_max = 0.5
nt = 151
dt = t_max/nt

# mesh
nx = 51
xmin = 0
xmax = 2
dx = (xmax-xmin)/(nx-1)
x = numpy.linspace(xmin, xmax, nx)

# Initialize dependent variable
u = numpy.ones(nx)
up = numpy.ones(nx)    # from previous time-step

# Initial condition: hat function
u[int(.5/dx):int(1/dx)+1]=2


def advance(t_steps):
    for n in range(0,t_steps): # loop over all time-steps
        up = u.copy()
        # with periodic BC all are internal points
        u[1:-1] = up[1:-1] + nu*dt/dx**2*(up[2:] - 2*up[1:-1] + up[0:-2])
        u[0] = up[0] + nu*dt/dx**2*(up[1] - 2*up[0] + up[-1]) 
        u[-1] = up[-1] + nu*dt/dx**2*(up[0] - 2*up[-1] + up[-2])


def plot_multi(x, u, incr):
    fig = pyplot.figure(figsize=(7, 5), dpi=100)
    colour = iter(cm.rainbow(numpy.linspace(0,incr,nt)))
    for i in range(0, nt, incr):
        c = next(colour)
        pyplot.plot(x,u, c=c)
        advance(incr)
        pyplot.xlabel('x (m)')
        pyplot.ylabel('u (m/s)')
        pyplot.xlim(0, 2)
        pyplot.ylim(0, 2.5)
        pyplot.title('Non-linear convection')
    pyplot.show()

# Let us plot some curves to check if everything is alright
plot_multi(x, u, 10)
