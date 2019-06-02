import random
from math import cos, pi, sin, sqrt, e, exp


def optifun_DeJong1(chromo):
    """
    Sphere model (DeJong1)
    хО(-5,12; 5,12)
    Один минимум равный 0 в точке, где xi=0.0.
    :param chromo: хромосома, у которой все гены DRange
    :return: fitness
    """
    summ = 0.0
    for value in chromo.values():
        summ += value ** 2
    return -summ


def optifun_DeJong2(chromo):
    """
    Rosenbrock's saddle (DeJong2)
    хО(-2,048; 2,048)
    Минимум равный 0 в точке, где xi=1.0. Функция имеет большое медленно убывающее плато.
    :param chromo: хромосома, у которой все гены DRange
    :return: fitness
    """
    x = list(chromo.values())
    summ = 0.0
    for i in range(len(x) - 1):
        summ += 100 * (x[i + 1] - x[i] ** 2) ** 2 + (x[i] - 1) ** 2
    return -summ


def optifun_DeJong3(chromo):
    """
    Step function (DeJong3)
    хО(-5,12; 5,12)
    |[xi]| - модуль целой части
    Минимум равный 0 в точке, где xi=0.0.
    :param chromo: хромосома, у которой все гены DRange
    :return: fitness
    """
    summ = 0.0
    for value in chromo.values():
        summ += abs(value)
    return -summ


def optifun_DeJong4(chromo):
    """
    Gaussian quartic (DeJong4)
    хО(1.5; 1.5)
    gauss(0,1) - функция, возвращающая случайную величину с нормальным распределением с мат. ожиданием в 0 и дисперсией равной 1.
    Минимум равный 0 в точке, где xi=0.0.
    :param chromo: хромосома, у которой все гены DRange
    :return: fitness
    """
    summ = 0.0
    for value in chromo.values():
        summ += value ** 4 + random.gauss(0, 1)
    return -summ


def optifun_Rastrigin(chromo):
    """
    Функция со сложным рельефом. Считается сложной для оптимизации.
    хО(-5,12; 5,12)
    Минимум равный 0 в точке, где xi=0.0.
    Локальный минимум в точке, где одна координата равна 1.0, а остальные равны 0.0
    :param chromo: хромосома, у которой все гены DRange
    :return: fitness
    """
    n = len(chromo.values())
    summ = 10.0 * n
    for x in chromo.values():
        summ += x ** 2 - 10 * cos(2 * pi * x)
    return -summ


def optifun_Schwefel(chromo):
    """
    Schwefel's (Sine Root) Function
    хО(-500; 500)
    Минимум равный 0 в точке, где xi=420.9687.
    Локальный минимум в точке, где одна координата равна -302.5232, а остальные равны 420.9687.
    Т.к. локальный минимум находится далеко от глобального, то есть опасность, что алгоритм "собъется с пути".
    Угол под знаком синуса получается в радианах.
    :param chromo: хромосома, у которой все гены DRange
    :return: fitness
    """
    n = len(chromo.values())
    summ = 418.9829 * n
    for x in chromo.values():
        summ += -x * sin(sqrt(abs(x)))
    return -summ


def optifun_Griewangk(chromo):
    """
    Griewangk's Function
    хО(-600; 600)
    Минимум равный 0 в точке, где xi=0.0.
    При n=10 есть ещё 4 локальных минимума равных 0,0074 приблизительно в точке (+-Pi, +-Pi*sqrt(2), 0, ..., 0).
    :param chromo: хромосома, у которой все гены DRange
    :return: fitness
    """
    summ = 1.0
    mult = 1.0
    for i, x in enumerate(chromo.values()):
        summ += x ** 2 / 4000
        mult *= cos(x / sqrt(i + 1))
    summ -= mult
    return -summ


def optifun_Ackley(chromo):
    """
    Ackley's Function
    хО(-30; 30)
    Минимум равный 0 в точке, где xi=0.0.
    :param chromo: хромосома, у которой все гены DRange
    :return: fitness
    """
    n = len(chromo.values())
    summ1, summ2 = 0.0, 0.0
    for x in chromo.values():
        summ1 += x ** 2
        summ2 += cos(2 * pi * x)
    summ = 20 + e - 20 * exp(-0.2 * sqrt(summ1 / n)) - exp(summ2 / n)
    return -summ


def get_funcs_diap():
    """
    Функция получения словаря со всеми необходимыми функциями и рекомендуемыми диапазонами
    :return: dict вида key = имя функции, value = tuple(ссылка на функцию, x0_1, x0_2, xi_min)
    """
    return {
        'optifun_DeJong1': (optifun_DeJong1, -5.12, 5.12, 0.0),
            'optifun_DeJong2': (optifun_DeJong2, -2.048, 2.048, 1.0),
            'optifun_DeJong3': (optifun_DeJong3, -5.12, 5.12, 0.0),
            'optifun_DeJong4': (optifun_DeJong4, -1.5, 1.5, 0.0),
            'optifun_Rastrigin': (optifun_Rastrigin, -5.12, 5.12, 0.0),
            'optifun_Schwefel': (optifun_Schwefel, -500, 500, 420.9687),
            'optifun_Griewangk': (optifun_Griewangk, -600, 600, 0.0),
            'optifun_Ackley': (optifun_Ackley, -10, 10, 0.0)
            }


def print_surf(name, func, x0, x1, n, xmax):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import cm

    # plt.xkcd()
    # Create heatmap
    xedges, yedges = np.linspace(x0, x1, n), np.linspace(x0, x1, n)
    extent = [xedges[0], xedges[-1], yedges[-1], yedges[0]]

    heatmap = np.zeros((n, n))
    for i, x in enumerate(xedges):
        for j, y in enumerate(yedges):
            # if i==j:
            #     heatmap[i, i] = 10
            #     continue
            heatmap[i, j] = func({1: x, 2: y})
    # Plot heatmap

    fig, ax = plt.subplots()
    plt.title(name + f" max({xmax}; {xmax}), maxValue = {np.max(heatmap)}")
    plt.ylabel('y')
    plt.xlabel('x')
    cax = ax.imshow(heatmap, extent=extent, cmap=cm.inferno)
    fig.colorbar(cax)
    ax.plot(xmax, xmax, '-o', color='r', linewidth=1)
    plt.show()


def main():
    # from tqdm import tqdm

    d = get_funcs_diap()
    n = 200
    for k in d:
        func, x0, x1, xmax = d[k]
        print_surf(k, func, x0, x1, n,xmax)

def main2():
    '''
    ======================
    3D surface (color map)
    ======================

    Demonstrates plotting a 3D surface colored with the coolwarm color map.
    The surface is made opaque by using antialiased=False.

    Also demonstrates using the LinearLocator and custom formatting for the
    z axis tick labels.
    '''

    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    import numpy as np


    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Make data.
    n = 1000
    x0= -2
    x1=2
    X = np.linspace(x0, x1, n)
    Y = np.linspace(x0, x1, n)
    X, Y = np.meshgrid(X, Y)
    Z = np.zeros_like(X)
    d = get_funcs_diap()
    cmap = cm.jet
    
    func, x00, x11, xmax = d['optifun_Ackley']
    for i, (x,y) in enumerate(zip(X, Y)):
        for j, (xx, yy) in enumerate(zip(x, y)):
            Z[i,j] = func({1: xx, 2: yy})+10

    # Plot the surface.
    surf = ax.plot_surface(X, Y, Z, cmap=cmap,
                        linewidth=0, antialiased=True)


    # Customize the z axis.
    ax.set_zlim(1, 12)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()


    # xedges, yedges = np.linspace(x0, x1, n), np.linspace(x0, x1, n)
    # extent = [xedges[0], xedges[-1], yedges[-1], yedges[0]]

    # heatmap = np.zeros((n, n))
    # for i, x in enumerate(xedges):
    #     for j, y in enumerate(yedges):
    #         # if i==j:
    #         #     heatmap[i, i] = 10
    #         #     continue
    #         heatmap[i, j] = func({1: x, 2: y})+10
    # # Plot heatmap

    # fig, ax = plt.subplots()
    # # plt.title(name + f" max({xmax}; {xmax}), maxValue = {np.max(heatmap)}")
    # plt.ylabel('y')
    # plt.xlabel('x')
    # cax = ax.imshow(heatmap, extent=extent, cmap=cmap)
    # fig.colorbar(cax)
    # # ax.plot(xmax, xmax, '-o', color='r', linewidth=1)
    # plt.show()

    plt.contour(X, Y, Z, levels=13, alpha=0.8)
    plt.show()



if __name__ == '__main__':
    main2()
    # from tqdm import tqdm
    # import time
    # for i in tqdm(range(10)):
    #     time.sleep(1)