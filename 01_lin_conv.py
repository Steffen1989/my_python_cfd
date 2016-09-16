import numpy
from matplotlib import pyplot

## Variable declaration
# parameters
c = 1

# time 
dt = .025

# mesh
nx = 81
xmin = 0
xmax = 2
dx = (xmax-xmin)/(nx-1)
x = numpy.linspace(xmin, xmax, nx)

# Initialize dependent variable
u = numpy.zeros(nx)
up = numpy.zeros(nx)    # from previous time-step

# Initial condition: hat function
u[int(.5/dx):int(1/dx)+1]=2

# check initial condition
def plot2d(x_plot, u_plot):
    fig = pyplot.figure(figsize=(5, 3), dpi=100)
    pyplot.xlim(0,2)
    pyplot.ylim(0,2.5)
    pyplot.plot(x_plot, u_plot)
    fig.show()

plot2d(x, u)

def advance(t_steps):
    for n in range(0,t_steps): # loop over all time-steps
        up = u.copy()
        # with periodic BC all are internal points
        u[1:] = up[1:] - c*dt/dx*(up[1:] - up[0:-1])
        u[0] = up[0] - c*dt/dx*(up[0] - up[-1])


# Let us plot some curves to check if everything is alright
advance(25)
plot2d(x, u)

advance(25)
plot2d(x, u)

advance(25)
plot2d(x, u)
