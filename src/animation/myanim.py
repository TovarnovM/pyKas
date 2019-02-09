import numpy as np
from matplotlib.pylab import *
import os
import imageio
import glob
from tqdm import tqdm

default_kwars = dict(
    n_nodes = 100,
    x0 = 0.5,
    x1 = 2,
    n_ts = 100,
    t0 = 0,
    t1 = 20,
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

def savefig(ts, xs, ps, i, fname=None, **kwargs):
    # f0 = figure(num = 0, figsize = (6, 4))#, dpi = 100)
    f0 = figure(num = 0, figsize = (14, 6))
    ax01 = subplot2grid((3, 7), (0, 0), colspan=2, rowspan=2)
    ax02 = subplot2grid((3, 7), (0, 5), colspan=2, rowspan=2)
    ax03 = subplot2grid((3, 7), (0, 2), colspan=3, rowspan=2)
    ax04 = subplot2grid((3, 7), (2, 2), colspan=3)

    # ax01.set_title('Position vs Time')
    # ax02.set_title('Velocity vs Time')
    # ax03.set_title('Position and Velocity vs Time')
    if 'xlim' in kwargs:
        xlim = kwargs['xlim']
        ax03.set_xlim(*xlim)
        ax04.set_xlim(*xlim)
    if 'plim' in kwargs:
        plim = kwargs['plim']
        ax03.set_ylim(*plim)
        ax01.set_ylim(*plim)
        ax02.set_ylim(*plim)
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

    ax03.plot(xs[i], ps[i],linewidth=2)
    

    def get_x_otn(ax, x):
        x0, x1 = ax.get_xlim()
        return (x-x0)/(x1-x0)

    def plot_sub(ax, ind_x, col, where):
        data = []
        tss = []
        for j in range(i+1):
            tss.append(ts[j])
            data.append(ps[j][ind_x])
        ax.plot(tss, data, color=col,linewidth=2)
        ax.scatter(tss[-1], data[-1], color=col)
        ax03.scatter(xs[i][ind_x], ps[i][ind_x], color=col)
        if where == 'left':
            ax.axhline(y=data[-1], xmin=get_x_otn(ax, tss[-1]), xmax=1.05,c=col,linewidth=1,zorder=0, clip_on=False)
            ax03.axhline(y=data[-1], xmin=-0.05, xmax=get_x_otn(ax03, xs[i][ind_x]),c=col,linewidth=1,zorder=0, clip_on=False)
        else:
            ax.axhline(y=data[-1], xmin=-0.05, xmax=get_x_otn(ax, tss[-1]),c=col,linewidth=1,zorder=0, clip_on=False)
            ax03.axhline(y=data[-1], xmin=get_x_otn(ax03, xs[i][ind_x]), xmax=1.05,c=col,linewidth=1,zorder=0, clip_on=False)

    for ind_x, col in zip(kwargs['inds_left'], kwargs['colors_left']):
        plot_sub(ax01, ind_x, col, where='left')
    for ind_x, col in zip(kwargs['inds_right'], kwargs['colors_right']):
        plot_sub(ax02, ind_x, col, where='right')    

    if fname is not None:
        f0.savefig(os.path.dirname(os.path.realpath(__file__))+'\\imgs\\'+fname)
    return f0

# savefig(ts,xs,ps,50,xlim=(0,25), plim=(0,4), inds_left=[0,33], inds_right=[-1,66], colors_left=['r','orange'], colors_right=['g','lime'])
# plt.show()

for i in tqdm(range(1, 100)):    
    savefig(ts,xs,ps, i, f'{i}.png',xlim=(0,25), plim=(0,4), inds_left=[0,33], inds_right=[-1,66], colors_left=['r','orange'], colors_right=['g','lime'])

fpattern = os.path.dirname(os.path.realpath(__file__))+'\\imgs\\*.png'
def get_file_index(fn):
    basename = os.path.basename(fn)
    index = int(os.path.splitext(basename)[0])
    return index
files = glob.glob(fpattern)
files.sort(key=get_file_index)

images = []
for file_name in tqdm(files):
    images.append(imageio.imread(file_name))
imageio.mimsave(os.path.dirname(os.path.realpath(__file__))+'\\movie.gif', images, duration=0.05)

for f in tqdm(files):
    os.remove(f)

color = '#FF0000'