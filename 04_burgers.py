import numpy
from matplotlib import pyplot
from matplotlib import cm

## Variable declaration
# parameters
nu = .1

# time 
tmax = 0.5
nt = 151
dt = tmax/(nt-1)

# mesh
nx = 151
xmin = 0
xmax = 2*numpy.pi
dx = (xmax-xmin)/(nx-1)
x = numpy.linspace(xmin, xmax, nx)

# Initialize dependent variable
u = numpy.zeros(nx)
up = numpy.zeros(nx)    # from previous time-step

# Initial condition: saw tooth
phi = numpy.exp(-x**2/(4*nu)) + numpy.exp(-(x-2*numpy.pi)**2/(4*nu))
dphi = - 0.5*x/nu * numpy.exp(-x**2/(4*nu)) - 0.5*(x-2*numpy.pi)/nu * numpy.exp(-(x-2*numpy.pi)**2/(4*nu))

u = -2*nu*dphi/phi + 4

# Numerical solution
def anal_sol(t):
    phi2 = numpy.exp(-(x-4*t)**2/(4*nu*(t+1))) + numpy.exp(-(x-4*t-2*numpy.pi)**2/(4*nu*(t+1)))
    dphi2 = - 0.5*(x-4*t)/(nu*(t+1)) * numpy.exp(-(x-4*t)**2/(4*nu*(t+1))) - 0.5*(x-4*t-2*numpy.pi)/(nu*(t+1)) *\
            numpy.exp(-(x-4*t-2*numpy.pi)**2/(4*nu*(t+1)))
    u_anal = -2*nu*dphi2/phi2 + 4

    return u_anal


def advance(t_steps):
    for n in range(0,t_steps): # loop over all time-steps
        up = u.copy()
        # with periodic BC all are internal points
        u[1:-1] = up[1:-1] - up[1:-1]*dt/dx*(up[1:-1] - up[0:-2]) + \
                nu*dt/dx**2*(up[2:]-2*up[1:-1]+up[0:-2])

        u[0] = up[0] - up[0]*dt/dx*(up[0] - up[-1]) + \
                nu*dt/dx**2*(up[1]-2*up[0]+up[-1])
        u[-1] = up[-1] - up[-1]*dt/dx*(up[-1] - up[-2]) + \
                nu*dt/dx**2*(up[0]-2*up[-1]+up[-2])


def plot_multi(x, u, incr):
    fig = pyplot.figure(figsize=(9, 5), dpi=100)
    ax = pyplot.subplot(111)
    colour = iter(cm.rainbow(numpy.linspace(0,incr,nt)))
    for i in range(0, nt, incr):
        c = next(colour)
        ax.plot(x, u, label='i='+str(i)+' numerical')
        ax.plot(x, anal_sol(i*dt), 'ko', c=c, markerfacecolor='none', label='i='+str(i)+' analytical')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.97, box.height])
        legend = ax.legend( bbox_to_anchor=(1.02,1), loc=2, fontsize=12)
        advance(incr)
        pyplot.xlabel('x (m)')
        pyplot.ylabel('u (m/s)')
        pyplot.xlim(0, 2*numpy.pi)
        pyplot.title('Burgers equation')
    pyplot.show()



# Let us plot some curves to check if everything is alright
plot_multi(x, u, 20)
