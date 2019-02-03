import numpy as np
from matplotlib.pylab import *

default_kwars = dict(
    n_nodes = 100,
    x0 = 0,
    x1 = 2,
    n_ts = 100,
    t0 = 0,
    t1 = 3,
    vel = 1,
    p0 = 2  
)
def get_data(**kwargs):
    n_nodes = kwargs['n_nodes']
    x0 = kwargs['x0']
    x1 = kwargs['x1']
    n_ts = kwargs['n_ts']
    t0 = kwargs['t0']
    t1 = kwargs['t1']
    vel = kwargs['vel']
    p0 = kwargs['p0']
    xs = []
    ps = []
    ts = []
    s = None
    for t in np.linspace(t0, t1, n_ts):
        ts.append(t)
        x2 = x1 + vel*t
        xss = np.linspace(x0, x2, n_nodes)
        xs.append(xss)
        pss = np.sin(xss+t) + p0
        s1 = np.sum(pss[:-1]*(xss[1:] - xss[:-1]))
        if s is None:
            s = s1
        pss *= s/s1
        ps.append(pss)
    return ts, xs, ps 

ts, xs, ps = get_data(**default_kwars)

# f0 = figure(num = 0, figsize = (6, 4))#, dpi = 100)
f0 = figure(num = 0, figsize = (7, 3))
ax01 = subplot2grid((3, 7), (0, 0), colspan=2, rowspan=2)
ax02 = subplot2grid((3, 7), (0, 5), colspan=2, rowspan=2)
ax03 = subplot2grid((3, 7), (0, 2), colspan=3, rowspan=2)
ax04 = subplot2grid((3, 7), (2, 2), colspan=3)

ax01.set_title('Position vs Time')
ax02.set_title('Velocity vs Time')
ax03.set_title('Position and Velocity vs Time')

plt.setp(ax03.get_xticklabels(), visible=False)
plt.setp(ax03.get_yticklabels(), visible=False)
plt.setp(ax04.get_yticklabels(), visible=False)

ax03.tick_params(axis='both', which='both', length=0)
ax04.tick_params(axis='y', length=0)

ax02.yaxis.tick_right()

ax01.grid(True)
ax02.grid(True)
ax03.grid(True)
ax04.grid(True)

plt.show()