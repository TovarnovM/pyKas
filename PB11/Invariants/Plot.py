import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Polygon
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def plot_one(x, y, name, name_x, name_y):
    """
    Функция построения графиков
    :param x:
    :param y:
    :param name: Название графика
    :param name_x: Название оси X
    :param name_y: Название оси Y
    :return:
    """
    plt.figure()
    plt.plot(x, y, color='red')
    plt.title(name)
    plt.xlabel(name_x)
    plt.ylabel(name_y)
    plt.grid(True)
    plt.show()


def plot_two(x, y1, y2, name, legend1, legend2):
    """

    :param x:
    :param y1:
    :param y2:
    :param name:
    :return:
    """
    plt.figure()
    plt.title(name)
    plt.plot(x, y1, label=legend1, color='red')
    plt.plot(x, y2, label=legend2, color='blue')
    plt.legend()
    plt.grid(True)
    plt.show()


def plot_contours(X, x, y, Z, name, name_x, name_y):
    plt.figure()
    x, Y = np.meshgrid(x, y)
    cp = plt.contourf(X, Y, Z, 100)
    plt.colorbar(cp)
    plt.title(name)
    plt.xlabel(name_x)
    plt.ylabel(name_y)
    plt.show()


def plot_3d(X, x, y, Z, name, name_x, name_y):
    fig = plt.figure()
    ax = Axes3D(fig)
    x, Y = np.meshgrid(x, y)
    surf = ax.plot_surface(X, Y, Z)
    # fig.colorbar(surf, shrink=0.75, aspect=15)
    plt.title(name)
    ax.set_xlabel(name_x)
    ax.set_ylabel(name_y)
    plt.show()

    # ax.title(name)


def plot_muzzle(layers):
    ax = plt.axes()
    tube = layers[0].tube
    xs, ys = tube.d.x, tube.d.y / 2
    line = plt.Line2D(xs, ys, marker='.', label='ствол')
    ax.add_line(line)
    line2 = plt.Line2D([xs[0], xs[-1]], [0, 0], linestyle='-.', linewidth=1)
    ax.add_line(line2)
    line3 = plt.Line2D([xs[0], xs[0]], [0, ys[0]], marker='.')
    ax.add_line(line3)
    line4 = plt.Line2D([xs[-1], xs[-1]], [0, ys[-1]], marker='.')
    ax.add_line(line4)
    ax.legend()
    plt.xlim((xs[0]-0.1, xs[-1]+0.1))
    plt.ylim((-0.05, max(ys)+0.05))
    plt.grid(True)
    plt.show()


def plot_layers_values(layers, attr_name, ax):
    lines = []
    cmap = matplotlib.pyplot.get_cmap('Set1')
    norm = matplotlib.colors.Normalize(vmin=0.,vmax=len(layers))
    for i, l in enumerate(layers):
        ys = getattr(l, attr_name)
        xs = l.x_c
        line = plt.Line2D(xs, ys, marker='.', label=f"сетка {i}", color=cmap(i))
        lines.append(ax.add_line(line))
    plt.legend()
    return lines
