import numpy as np
from scipy import interpolate
from math import pi, ceil, sqrt
import sys
import time
import random as rnd

from Invariants.Plot import *
from Invariants.Tube import Tube
from Pneumatics.PnInit import pn_create_layer
from OneVelocity.OvInit import ov_create_layer
from ElasticPiston.ElPistInit import el_pist_create_layer

powders_dict = {
    '11\\7': {'I_k': 0.68,
              'T_1': 2855.0,
              'Z_k': 1.536,
              'alpha_k': 1.015,
              'etta': 0.227,
              'f': 0.998,
              'k_1': 0.715,
              'k_2': 0.549,
              'k_f': 0.0003,
              'k_l': 0.0016,
              'lambda_1': 0.194,
              'lambda_2': -0.933,
              'name': '11\\7',
              'ro': 1.6},
    '12\\7': {'I_k': 0.81,
              'T_1': 2890.0,
              'Z_k': 1.522,
              'alpha_k': 1.013,
              'etta': 0.226,
              'f': 1.004,
              'k_1': 0.725,
              'k_2': 0.546,
              'k_f': 0.0003,
              'k_l': 0.0016,
              'lambda_1': 0.183,
              'lambda_2': -0.957,
              'name': '12\\7',
              'ro': 1.6},
    '14\\7': {'I_k': 0.96,
              'T_1': 2865.0,
              'Z_k': 1.535,
              'alpha_k': 1.015,
              'etta': 0.227,
              'f': 1.0,
              'k_1': 0.728,
              'k_2': 0.541,
              'k_f': 0.0003,
              'k_l': 0.0016,
              'lambda_1': 0.174,
              'lambda_2': -0.934,
              'name': '14\\7',
              'ro': 1.6},
    '15\\7': {'I_k': 0.93,
              'T_1': 2915.0,
              'Z_k': 1.541,
              'alpha_k': 1.01,
              'etta': 0.225,
              'f': 1.01,
              'k_1': 0.73,
              'k_2': 0.539,
              'k_f': 0.0003,
              'k_l': 0.0016,
              'lambda_1': 0.17,
              'lambda_2': -0.924,
              'name': '15\\7',
              'ro': 1.6},
    '4\\1': {'I_k': 0.3,
             'T_1': 2970.0,
             'Z_k': 1.0,
             'alpha_k': 1.001,
             'etta': 0.221,
             'f': 1.014,
             'k_1': 1.07,
             'k_2': 0.0,
             'k_f': 0.0003,
             'k_l': 0.0016,
             'lambda_1': -0.066,
             'lambda_2': 0.0,
             'name': '4\\1',
             'ro': 1.6},
    '4\\1 фл': {'I_k': 0.4,
                'T_1': 2874.0,
                'Z_k': 1.069,
                'alpha_k': 1.026,
                'etta': 0.233,
                'f': 1.006,
                'k_1': 0.752,
                'k_2': 0.961,
                'k_f': 0.0003,
                'k_l': 0.0016,
                'lambda_1': 0.239,
                'lambda_2': 0.374,
                'name': '4\\1 фл',
                'ro': 1.6},
    '4\\7': {'I_k': 0.32,
             'T_1': 3006.0,
             'Z_k': 1.488,
             'alpha_k': 1.008,
             'etta': 0.228,
             'f': 1.027,
             'k_1': 0.811,
             'k_2': 0.505,
             'k_f': 0.0003,
             'k_l': 0.0016,
             'lambda_1': 0.081,
             'lambda_2': -1.024,
             'name': '4\\7',
             'ro': 1.6},
    '4\\7 св': {'I_k': 0.33,
                'T_1': 2956.0,
                'Z_k': 1.488,
                'alpha_k': 1.014,
                'etta': 0.229,
                'f': 1.018,
                'k_1': 0.811,
                'k_2': 0.505,
                'k_f': 0.0003,
                'k_l': 0.0016,
                'lambda_1': 0.081,
                'lambda_2': -1.024,
                'name': '4\\7 св',
                'ro': 1.6},
    '5\\1': {'I_k': 0.25,
             'T_1': 2930.0,
             'Z_k': 1.0,
             'alpha_k': 1.005,
             'etta': 0.223,
             'f': 1.009,
             'k_1': 1.062,
             'k_2': 0.0,
             'k_f': 0.0003,
             'k_l': 0.0016,
             'lambda_1': -0.058,
             'lambda_2': 0.0,
             'name': '5\\1',
             'ro': 1.6},
    '5\\1 х-10': {'I_k': 0.5,
                  'T_1': 2440.0,
                  'Z_k': 1.0,
                  'alpha_k': 1.026,
                  'etta': 0.257,
                  'f': 0.882,
                  'k_1': 1.045,
                  'k_2': 0.0,
                  'k_f': 0.0007,
                  'k_l': 0.0022,
                  'lambda_1': -0.043,
                  'lambda_2': 0.0,
                  'name': '5\\1 х-10',
                  'ro': 1.6},
    '5\\1 х-20': {'I_k': 0.7,
                  'T_1': 1890.0,
                  'Z_k': 1.0,
                  'alpha_k': 0.962,
                  'etta': 0.286,
                  'f': 0.713,
                  'k_1': 1.045,
                  'k_2': 0.0,
                  'k_f': 0.0007,
                  'k_l': 0.0022,
                  'lambda_1': -0.043,
                  'lambda_2': 0.0,
                  'name': '5\\1 х-20',
                  'ro': 1.6},
    '6\\7 гр': {'I_k': 0.5,
                'T_1': 2923.0,
                'Z_k': 1.498,
                'alpha_k': 1.02,
                'etta': 0.231,
                'f': 1.015,
                'k_1': 0.794,
                'k_2': 0.513,
                'k_f': 0.0003,
                'k_l': 0.0016,
                'lambda_1': 0.098,
                'lambda_2': -1.004,
                'name': '6\\7 гр',
                'ro': 1.6},
    '7\\1': {'I_k': 0.5,
             'T_1': 2930.0,
             'Z_k': 1.0,
             'alpha_k': 1.01,
             'etta': 0.223,
             'f': 1.009,
             'k_1': 1.059,
             'k_2': 0.0,
             'k_f': 0.0003,
             'k_l': 0.0016,
             'lambda_1': -0.056,
             'lambda_2': 0.0,
             'name': '7\\1',
             'ro': 1.6},
    '7\\14': {'I_k': 0.4,
              'T_1': 2962.0,
              'Z_k': 1.46,
              'alpha_k': 1.013,
              'etta': 0.229,
              'f': 1.021,
              'k_1': 0.581,
              'k_2': 0.091,
              'k_f': 0.0003,
              'k_l': 0.0016,
              'lambda_1': 0.289,
              'lambda_2': 1.087,
              'name': '7\\14',
              'ro': 1.6},
    '7\\7': {'I_k': 0.5,
             'T_1': 2900.0,
             'Z_k': 1.607,
             'alpha_k': 1.01,
             'etta': 0.224,
             'f': 1.006,
             'k_1': 0.769,
             'k_2': 0.506,
             'k_f': 0.0003,
             'k_l': 0.0016,
             'lambda_1': 0.101,
             'lambda_2': -0.823,
             'name': '7\\7',
             'ro': 1.6},
    '8\\1 тр': {'I_k': 0.7,
                'T_1': 2850.0,
                'Z_k': 1.0,
                'alpha_k': 1.017,
                'etta': 0.227,
                'f': 0.998,
                'k_1': 1.002,
                'k_2': 0.0,
                'k_f': 0.0003,
                'k_l': 0.0016,
                'lambda_1': -0.002,
                'lambda_2': 0.0,
                'name': '8\\1 тр',
                'ro': 1.6},
    '8\\7': {'I_k': 0.5,
             'T_1': 2855.0,
             'Z_k': 0.504,
             'alpha_k': 1.015,
             'etta': 0.227,
             'f': 0.998,
             'k_1': 0.783,
             'k_2': 0.542,
             'k_f': 0.0003,
             'k_l': 0.0016,
             'lambda_1': 0.17,
             'lambda_2': -0.993,
             'name': '8\\7',
             'ro': 1.6},
    '9\\7': {'I_k': 0.6,
             'T_1': 2855.0,
             'Z_k': 1.526,
             'alpha_k': 1.015,
             'etta': 0.227,
             'f': 0.998,
             'k_1': 0.724,
             'k_2': 0.545,
             'k_f': 0.0003,
             'k_l': 0.0016,
             'lambda_1': 0.183,
             'lambda_2': -0.951,
             'name': '9\\7',
             'ro': 1.6}}

powder_names = list(powders_dict.keys())


def shtraf(p, inv_dict):
    p_shtraf = inv_dict['p_shtraf']
    if p < p_shtraf:
        return 0
    p_max = inv_dict['p_max']
    v_shtraf = inv_dict['v_shtraf']
    k = v_shtraf / (p_max - p_shtraf)
    b = -k * p_shtraf
    return k * p + b


#  ['4\\1', '4\\1 фл', '4\\7', '4\\7 св','5\\1', '5\\1 х-10', '5\\1 х-20','6\\7 гр','7\\1', '7\\7','7\\14','8\\1 тр','8\\7','9\\7','11\\7','12\\7','14\\7','15\\7']




def main_tube_w_tst():
    from matplotlib import pyplot as plt
    t = Tube([1, 3, 4], [2, 1, 2])

    x1, x2 = 0.5, 3.5
    w = t.get_W_between(x1, x2)
    x22 = t.get_x2(x1, w)
    print(f"w = {w}, x2 = {x22}")

    xs = [0, 1, 2, 3, 4, 5]
    ys = [t.w(x) for x in xs]
    plt.plot(xs, ys)
    plt.grid()
    plt.show()


def get_powder_grid(cr, tube, x_left, i):
    ro = cr[f"delta_powder_{i}"]
    cell_dx = 0.0015
    omega = cr[f"om_{i}"]
    w = omega / ro
    x_right = float(tube.get_x2(w=w, x1=x_left))
    n_cells = ceil((x_right - x_left) / cell_dx)
    const_powder = {'covolume': 0,
                    'R': 287,
                    'gamma': 1.228,
                    'param_powder': dict(powders_dict[cr[f"powder_type_{i}"]]),
                    'nu': 1}
    init_const_powder = {'ro': ro}
    grid_powder = {'name': 'powder',
                   'n_cells': n_cells,
                   'consts': const_powder,
                   'type': 'powder',
                   'init_const': init_const_powder,
                   'xl': x_left,
                   'xr': x_right,
                   'vl': 0.0,
                   'vr': 0.0}

    return grid_powder, x_right


def get_tube(inv_dict):
    d_0 = inv_dict['d_0']
    l_0 = inv_dict['l_0']
    l_05 = inv_dict['l_05']
    d_1 = inv_dict['d_1']
    l_1 = inv_dict['l_1']
    xs = [0, l_0, l_0 + l_05, l_0 + l_05 + l_1]
    ds = [d_0, d_0, d_1, d_1]
    return Tube(xs, ds)


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


def get_borders(layers, tau, l_min, m1, m2, pf_1):
    p0, p1, p2, pa = layers[0].p[-1], layers[1].p[0], layers[1].p[-1], 101325.0
    tube = layers[0].tube
    x1, v1 = layers[0].x[-1], layers[0].V[-1]
    x2, v2 = layers[1].x[-1], layers[1].V[-1]
    s1, s2 = tube.get_S([x1, x2])
    f1, f2 = (p0 - p1) * s1, (p2 - pa) * s2
    if abs(p0 - p1) > pf_1 or abs(v1) > 1E-2:
        x1_, v1_ = newton(tau, f1, m1, v1, x1)
    else:
        x1_, v1_ = x1, v1
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


def get_bord_before_ign2(layers, tau, m_sum, pf_1):
    p0, pa = layers[0].p[-1], 101325.0
    tube = layers[0].tube
    x1, v1 = layers[0].x[-1], layers[0].V[-1]
    s1 = float(tube.get_S(x1))
    f1 = (p0 - pa) * s1
    if abs(p0 - pa) > pf_1 or abs(v1) > 1E-2:
        x1_, v1_ = newton(tau, f1, m_sum, v1, x1)
    else:
        x1_, v1_ = x1, v1
    return x1_, v1_


def main(param_dict, with_calc=True, print_gr=True):
    try:
        wchr = param_dict['wchr']

        inv_dict = param_dict['inv_dict']
        tube = get_tube(inv_dict)

        solv_dict1, x_r1 = get_powder_grid(wchr['chromo'], tube, 0, 1)
        solv_dict2, x_r2 = get_powder_grid(wchr['chromo'], tube, x_r1, 2)

        lr1 = ov_create_layer(solv_dict1, tube)
        lr2 = ov_create_layer(solv_dict2, tube)

        layers = [lr1, lr2]
        start_time = time.time()

        l_max = inv_dict['l_0'] + inv_dict['l_05'] + inv_dict['l_1']
        time_max = inv_dict['time_max']
        time_max_calc = inv_dict['time_max_calc']
        time_v0 = inv_dict['time_v0']
        Vmax_curr = 0
        pmax_curr = 0

        p_max = inv_dict['p_max']
        pf_1 = wchr['chromo']['pf_1']
        pmax_curr = 0
        courant_number = inv_dict['courant_number']

        m2 = inv_dict['m_elem']
        omega2 = wchr['chromo']["om_2"]
        l_min = x_r2 - x_r1
        m1 = wchr['chromo']['m_1']
        l_ignite = wchr['chromo']['l_ignite']
        p_ign_perc = wchr['chromo']['p_ign_perc']

        ind = 0
        time_arr = [0]  # Список времени для графика
        V_arr = [0]  # Список скорости для графика
        V_podd = [0]
        x_arr = [layers[-1].x[-1]]  # Список координаты для графика
        p_arr1 = [layers[0].p[-1]]
        p_arr = [layers[-1].p[-1]]

        p_arr_3d = np.append(layers[0].p, layers[1].p)
        x_c_arr_3d = np.append(layers[0].x_c, layers[1].x_c)
        time_arr_3d = [0]

        isignite = False

        if not with_calc:
            v = (wchr['vec_param'] - 0.5)*5.12
            # res = -sum(v * v)
            # return res
            summ = 10*len(v)
            summ += sum(v*v-10*np.cos(2*pi*v))
            return -summ

        while True:
            if l_max < layers[-1].x[-1]:
                break
            # if (Vmax_curr - layers[-1].V[-1]) > 20:
            #     break
            if layers[0].time >= time_max:
                return -9999
            # if pmax_curr > p_max:
            #     return -9998
            if time.time() - start_time > time_max_calc:
                return -9997
            if (time_v0 < layers[0].time) and (layers[-1].V[-1] < 1E-4):
                return -9995

            layers1 = []

            if not isignite:
                if layers[1].x[-1] > l_ignite:
                    isignite = True
                    p_lefty = layers[0].p[-1]
                    p_ign = p_lefty * p_ign_perc
                    print(f"p_lefty = {p_lefty},   p_ign = {p_ign}")
                    layers[1].p[:] = p_ign if p_ign > 101325.0 else 101325.0
                    layers[1].init_arr_q()
                    layers[1].y = layers[1].powd.psi(layers[1].z)
                    layers[1].e = layers[1].get_energ(layers[1].p, layers[1].ro, layers[1].y)

                    layers[1].q = layers[1].init_arr_q()  # Список нумпи массивов q1, q2, q3, размерностями n

                    layers[1].h = layers[1].get_arr_h()  # Список нумпи массивов h1, h2, h3, размерностями n

            if isignite:
                tau_arr = [l.time_step() for l in layers]
                tau = courant_number * min(tau_arr)  # Вычисление шага по времени

                x1, v1, x2, v2 = get_borders(layers, tau, l_min, m1, m2, pf_1)

                layers1.append(layers[0].euler_step_new(tau, 0, 0, x1, v1, True))
                layers1.append(layers[1].euler_step_new(tau, x1, v1, x2, v2, isignite))
            else:
                tau_arr = [l.time_step() for l in layers[:-1]]
                tau = courant_number * min(tau_arr)  # Вычисление шага по времени
                x1, v1 = get_bord_before_ign2(layers, tau, m1 + m2 + omega2, pf_1)
                x2, v2 = x1 + l_min, v1
                layers1.append(layers[0].euler_step_new(tau, 0, 0, x1, v1, True))
                layers1.append(layers[1].move_to(tau, x1, v1, x1 + l_min))

            layers = layers1

            if layers[-1].V[-1] > Vmax_curr:
                Vmax_curr = layers[-1].V[-1]
            pmax_curr1 = max([max(layer.p) for layer in layers])
            if pmax_curr < pmax_curr1:
                pmax_curr = pmax_curr1

            if print_gr:
                time_arr.append(layers[-1].time)  # Добавление текущего шага по времени в список для графика
                V_arr.append(layers[-1].V[-1])  # Добавление текущей скорости поршня в список для графика
                V_podd.append(layers[0].V[-1])
                x_arr.append(layers[-1].x[-1])  # Добавление текущей координаты поршня в список для графика
                p_arr.append(layers[-1].p[-1])
                p_arr1.append(layers[0].p[-1])

                ind += 1
                if ind % 300 == 0:
                    ro, u, e, z = layers[0].get_param(layers[0].q)
                    print(layers[-1].time, Vmax_curr, pmax_curr, x2, isignite)
                    x_c_arr_3d = np.vstack([x_c_arr_3d, np.append(layers[0].x_c, layers[1].x_c)])
                    p_arr_3d = np.vstack([p_arr_3d, np.append(layers[0].p, layers[1].p)])
                    time_arr_3d.append(layers[-1].time)

        if print_gr:
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
    except Exception as e:
        print(e)
        return -9996


def fitness_foo(param_dict):
    Vmax_curr, pmax_curr = main(param_dict, with_calc=True, print_gr=False)
    fitness = Vmax_curr - shtraf(pmax_curr, param_dict['inv_dict'])
    return fitness


def fit_fun4tst(param_dict):
    time.sleep(rnd.random())
    return main(param_dict, with_calc=False, print_gr=True)


def main2():
    invar_dict = {
        'm_elem': 40E-3,
        'time_max': 0.5,
        # sec
        'time_max_calc': 21 * 60,
        'time_v0': 0.01,

        'd_0': 0.03,
        'd_1': 0.03,
        'l_0': 1.7,
        'l_05': 0.3,
        'l_1': 3.7,

        'p_shtraf': 4000E5,
        'p_max': 7000E5,
        'v_shtraf': 3000,

        'courant_number': 0.13
    }
    param_dict = {'wchr': {'chromo': {'delta_powder_1': 675.4713285462944,
  'delta_powder_2': 595.122185293749,
  'l_ignite': 0.6451790850223771,
  'm_1': 0.0791291959059027,
  'om_1': 0.08430337660291476,
  'om_2': 0.09629952057375057,
  'p_ign_perc': 0.5806239465502379,
  'pf_1': 12358973.721830292,
  'powder_type_1': '5\\1',
  'powder_type_2': '5\\1'},
 'dop_info': None,
 'fitness': 2988.0375023761658,
 'id': 0,
 'name': 'Любомирот',
 'vec_param': np.array([ 0.82559307,  0.95888356,  0.32849553,  0.61602883,  0.43867832,
         0.23780546,  0.2150597 ,  0.60077993]),
 'vec_struct': np.array([ 4.,  4.])},
                  'inv_dict': invar_dict}
    print(fit_fun4tst(param_dict))


main2()
