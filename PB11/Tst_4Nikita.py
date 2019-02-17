import numpy as np
from scipy import interpolate
from math import pi, ceil, sqrt
import sys
import time
import random as rnd

from Invariants.Tube import Tube

from Invariants.Plot import *
from Pneumatics.PnInit import pn_create_layer

consts_gas = {
    'air': {
        'covolume': 0.0010838,
        'R': 287,
        'gamma': 1.4},
    'He': {
        'covolume': 0.005925,
        'R': 2078,
        'gamma': 1.66},
    'CO_2': {
        'covolume': 0.0009702,
        'R': 189,
        'gamma': 1.3}

}

def get_tube(par_dict):
    d_0 = par_dict['d_0']
    l_0 = par_dict['l_0']
    xs = [0, l_0]
    ds = [d_0, d_0]
    return Tube(xs, ds)

def get_gas_grid_and_rborder(par_dict, tube, x_left):
    # 'gas_1': {
    #     'p': 200E5,
    #     'l': 0.3,
    #     'T': 300,
    #     'consts': {
    #         'covolume': 0.005925,
    #         'R': 2078,
    #         'gamma': 1.66
    #     }
    # },
    cell_dx = 0.001
    l = par_dict["l"]
    n_cells = ceil(l / cell_dx)
    x_right = x_left + l

    init_const_gas = {'p': par_dict["p"],
                      'T': par_dict["T"]}
    grid_gas = {'name': 'gas',
                'n_cells': n_cells,
                'consts': par_dict["consts"],
                'type': 'gas',
                'init_const': init_const_gas,
                'xl': x_left,
                'xr': x_right,
                'vl': 0.0,
                'vr': 0.0
                }

    return grid_gas, x_right

def func_in(x, grid):
    """
    Функция инициализации ячейки
    :param x: координата центра ячейки
    :return: значения параметров плотности, скорости и давления в ячейке
    """
    if x:
        ro = grid['init_const']['p'] / (grid['consts']['R'] * grid['init_const']['T'] +
                                        grid['init_const']['p'] * grid['consts']['covolume'])
        # ro = grid['init_const']['ro']
        # ro = 250
        u = 0
        p = grid['init_const']['p']
        return ro, u, p

def newton(tau, f, m, v0, x0):
    a = f / m
    return x0 + v0 * tau + 0.5 * a * tau ** 2, v0 + a * tau

def root_find(a, b, c, t1, t2):
    d = b ** 2 - 4 * a * c
    if d < 0:
        return t1
    d = sqrt(d)
    x1 = (-b + d) / (2 * a)
    x2 = (-b - d) / (2 * a)
    if t1 <= x1 <= t2:
        if t1 <= x2 <= t2:
            return min(x1, x2)
        return x1
    elif t1 <= x2 <= t2:
        return x2
    return max(x1, x2, t1)

def get_borders(layers, tau, l_min, m1, m2):
    p0, p1, p2, pa = layers[0].p[-1], layers[1].p[0], layers[1].p[-1], 101325.0
    tube = layers[0].tube
    x1, v1 = layers[0].x[-1], layers[0].V[-1]
    x2, v2 = layers[1].x[-1], layers[1].V[-1]
    s1, s2 = tube.get_S([x1, x2])
    f1, f2 = (p0 - p1) * s1, (p2 - pa) * s2
    x1_, v1_ = newton(tau, f1, m1, v1, x1)
    x2_, v2_ = newton(tau, f2, m2, v2, x2)
    if x2 - x1 >= l_min:
        return x1_, v1_, x2_, v2_

    tau_0 = root_find(0.5 * (f2 / m2 - f1 / m1), v2 - v1, x2 - x1 - l_min, 0, tau)
    x1, v1 = newton(tau_0, f1, m1, v1, x1)
    x2, v2 = newton(tau_0, f2, m2, v2, x2)
    v_s = (v1 * m1 + v2 * m2) / (m1 + m2)
    x1_, v1_ = newton(tau - tau_0, f1 + f2, m1 + m2, v_s, x1)
    x2_, v2_ = newton(tau - tau_0, f1 + f2, m1 + m2, v_s, x2)
    return x1_, v1_, x2_, v2_


def get_bord_before_ign2(layers, tau, m_sum):
    p0, pa = layers[0].p[-1], 101325.0
    tube = layers[0].tube
    x1, v1 = layers[0].x[-1], layers[0].V[-1]
    s1 = float(tube.get_S(x1))
    f1 = (p0 - pa) * s1
    x1_, v1_ = newton(tau, f1, m_sum, v1, x1)
    return x1_, v1_

def main(par_dict, only1 = False):
    tube = get_tube(par_dict)

    solv_dict1, x_r1 = get_gas_grid_and_rborder(par_dict['gas_1'], tube, 0)
    solv_dict2, x_r2 = get_gas_grid_and_rborder(par_dict['gas_2'], tube, x_r1)


    lr1 = pn_create_layer(solv_dict1, tube)
    lr2 = pn_create_layer(solv_dict2, tube)

    layers = [lr1, lr2]
    start_time = time.time()

    l_max = par_dict['l_0']
    courant_number = par_dict['courant_number']

    m2 = par_dict['m_elem']
    l_min = x_r2 - x_r1
    m1 = par_dict['m_1']
    l_ignite = par_dict['l_ignite']


    ind = 0
    time_arr = [0]  # Список времени для графика
    V_arr = [0]  # Список скорости для графика
    V_podd = [0]
    x_s = layers[0].x[-1] if only1 else layers[-1].x[-1]
    x_arr = [x_s]  # Список координаты для графика
    p_arr1 = [layers[0].p[-1]]
    p_arr = [layers[-1].p[-1]]

    p_arr_3d = np.append(layers[0].p, layers[1].p)
    x_c_arr_3d = np.append(layers[0].x_c, layers[1].x_c)
    time_arr_3d = [0]

    Vmax_curr = 0
    pmax_curr = 0

    ro, u ,p  = lr1.func_in(x_r2, solv_dict2)
    v = tube.get_W_between(x_r1, x_r2)
    omega2 = ro*v

    m_sum = m1 + m2 + omega2 if not only1 else m2

    isignite = False
    while True:
        if l_max < x_s:
            break

        layers1 = []

        if not isignite:
            if (not only1) and layers[1].x[-1] > l_ignite:
                isignite = True
                layers[1].init_arr_q()
                layers[1].e = layers[1].get_energ(layers[1].p, layers[1].ro)

                layers[1].q = layers[1].init_arr_q()  # Список нумпи массивов q1, q2, q3, размерностями n

                layers[1].h = layers[1].get_arr_h()  # Список нумпи массивов h1, h2, h3, размерностями n

        if isignite:
            tau_arr = [l.time_step() for l in layers]
            tau = courant_number * min(tau_arr)  # Вычисление шага по времени

            x1, v1, x2, v2 = get_borders(layers, tau, l_min, m1, m2)

            layers1.append(layers[0].euler_step_new(tau, 0, 0, x1, v1))
            layers1.append(layers[1].euler_step_new(tau, x1, v1, x2, v2))
        else:
            tau_arr = [l.time_step() for l in layers[:-1]]
            tau = courant_number * min(tau_arr)  # Вычисление шага по времени
            x1, v1 = get_bord_before_ign2(layers, tau, m_sum)
            x2, v2 = x1 + l_min, v1
            layers1.append(layers[0].euler_step_new(tau, 0, 0, x1, v1))
            layers1.append(layers[1].move_to(tau, x1, v1, x1 + l_min))

        layers = layers1

        if layers[-1].V[-1] > Vmax_curr:
            Vmax_curr = layers[-1].V[-1]
        pmax_curr1 = max([max(layer.p) for layer in layers])
        if pmax_curr < pmax_curr1:
            pmax_curr = pmax_curr1

        x_s = layers[0].x[-1] if only1 else layers[-1].x[-1]
        time_arr.append(layers[-1].time)  # Добавление текущего шага по времени в список для графика
        V_arr.append(layers[-1].V[-1])  # Добавление текущей скорости поршня в список для графика
        V_podd.append(layers[0].V[-1])
        x_arr.append(x_s)  # Добавление текущей координаты поршня в список для графика
        p_arr.append(layers[-1].p[-1])
        p_arr1.append(layers[0].p[-1])

        ind += 1
        if ind % 300 == 0:
            ro, u, e = layers[0].get_param(layers[0].q)
            print(layers[-1].time, Vmax_curr, pmax_curr, x2, isignite)
            x_c_arr_3d = np.vstack([x_c_arr_3d, np.append(layers[0].x_c, layers[1].x_c)])
            p_arr_3d = np.vstack([p_arr_3d, np.append(layers[0].p, layers[1].p)])
            time_arr_3d.append(layers[-1].time)


    print("--- %s seconds ---" % (time.time() - start_time))

    print('Время вылета:', time_arr[-1], 'с')
    print('Скорость вылета:', V_arr[-1], 'м/с')
    print('Скорость вылета2:', layers[-1].u[-1], 'м/с')

    plot_two(time_arr, V_arr, V_podd, 'График скорости снаряда от времени', 'снаряд', 'поддон')
    plot_one(x_arr, V_arr, 'График скорости снаряда от времени', 'Координата', 'Скорость')
    plot_two(time_arr, p_arr, p_arr1, 'График давления на дно снаряда от времени', 'снаряд', 'сборка')
    plot_contours(x_c_arr_3d, x_c_arr_3d[0], time_arr_3d, p_arr_3d, 'Распределение давления', 'Координата',
                  'Время')
    plot_3d(x_c_arr_3d, x_c_arr_3d[0], time_arr_3d, p_arr_3d, 'Распределение давления на дно снаряда',
            'Координата', 'Время')

    return Vmax_curr, pmax_curr

if __name__ == '__main__':
    gas_type = 'air'
    p0 = 300E5
    par_dict = {
        'd_0': 0.05,
        'l_0': 3,
        'gas_1': {
            'p': p0,
            'l': 0.3,
            'T': 300,
            'consts': consts_gas[gas_type]
        },
        'gas_2': {
            'p': p0,
            'l': 0.1,
            'T': 300,
            'consts': consts_gas[gas_type]
        },

        'm_elem': 23E-3,
        'm_1': 10E-3,
        'l_ignite': 2,

        'courant_number': 0.13
    }




    par_dict2 = {
        'd_0': 0.05,
        'l_0': 3,
        'gas_1': {
            'p': p0,
            'l': par_dict['gas_1']['l'] + par_dict['gas_2']['l'] ,
            'T': 300,
            'consts': consts_gas[gas_type]
        },
        'gas_2': {
            'p': p0,
            'l': 0.01,
            'T': 300,
            'consts': consts_gas[gas_type]
        },

        'm_elem': 23E-3,
        'm_1': 30E-3,
        'l_ignite': 0.5,

        'courant_number': 0.13
    }

    main(par_dict)
    main(par_dict2, True)