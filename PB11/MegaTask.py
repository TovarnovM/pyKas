import numpy as np
from scipy import interpolate
from math import pi, ceil
import sys
import time
import random as rnd

from Invariants.Tube import Tube

sys.path.append('E:\\PyBallistic9')
sys.path.append('D:\\PyBallistic9')

from Invariants.Plot import *
from Pneumatics.PnInit import pn_create_layer
from OneVelocity.OvInit import ov_create_layer
from ElasticPiston.ElPistInit import el_pist_create_layer

from functional import seq

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

m_elem = 23E-3
time_max = 0.5
# sec
time_max_calc = 21 * 60
time_v0 = 0.007

l_max_mini = 7


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


def get_powder_grid_and_rborder(cr, tube, x_left):
    ro = cr["delta_powder"]
    cell_dx = 0.005
    omega = cr["omega_powder"]
    w = omega / ro
    x_right = float(tube.get_x2(w=w, x1=x_left))
    n_cells = ceil((x_right - x_left) / cell_dx)
    powder = dict(powders_dict[cr["powder_type"]])
    const_powder = {'covolume': 0,
                    'R': 287,
                    'gamma': 1.228,
                    'param_powder': powder,
                    'nu': 1}
    init_const_powder = {'ro': ro,
                         't_init': cr["t_powder"]}
    grid_powder = {'name': 'powder',
                   'n_cells': n_cells,
                   'consts': const_powder,
                   'type': 'powder',
                   'init_const': init_const_powder}

    rborder = {'m': cr["mass_bord_right"],
               'p_f': cr["p_bord_force_right"],
               'x': x_right,
               'V': 0}
    return grid_powder, rborder


def get_piston_grid_and_rborder(cr, tube, x_left):
    ro = 921.0
    cell_dx = 0.0025
    omega = cr["omega_piston"]
    w = omega / ro
    x_right = float(tube.get_x2(w=w, x1=x_left))
    n_cells = ceil((x_right - x_left) / cell_dx)
    const_piston = {'covolume': 0,
                    'gamma': 1.63098,
                    'c0': 2380,
                    'ro0': 919.03,
                    'sigmas': 25.2e6,
                    'k0': 0.054,
                    'b1': 0.027,
                    'b2': 0.00675,
                    'mu': 0.001,
                    'taus': 1e6}
    grid_piston = {'name': 'piston',
                   'n_cells': n_cells,
                   'consts': const_piston,
                   'type': 'piston'}

    rborder = {'m': cr["mass_bord_right"],
               'p_f': cr["p_bord_force_right"],
               'x': x_right,
               'V': 0}
    return grid_piston, rborder


def get_gas_grid_and_rborder(cr, tube, x_left):
    cell_dx = 0.005
    l = cr["l_gas"]
    n_cells = ceil(l / cell_dx)
    x_right = x_left + l
    type_gas = cr['type_gas']

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
    init_const_gas = {'p': cr["p_gas"],
                      'T': 300}
    grid_gas = {'name': 'gas',
                'n_cells': n_cells,
                'consts': consts_gas[type_gas],
                'type': 'gas',
                'init_const': init_const_gas}
    rborder = {'m': cr["mass_bord_right"],
               'p_f': cr["p_bord_force_right"],
               'x': x_right,
               'V': 0}
    return grid_gas, rborder


def wchr_to_solver(wchr):
    cr = wchr['chromo']
    i = 0
    xs, ds = [cr[f"l_{i}"]], [cr[f"d_{i}"]]
    i += 1
    while f"d_{i}" in cr:
        li, di = cr[f"l_{i}"], cr[f"d_{i}"]
        xs.append(xs[-1] + li)
        ds.append(di)
        i += 1
    geom = list(zip(xs, ds))
    tube = Tube(xs, ds)
    fun_dict = {
        "powder": get_powder_grid_and_rborder,
        "piston": get_piston_grid_and_rborder,
        "gas": get_gas_grid_and_rborder}
    borders = [
        {'m': 10000,
         'p_f': 1e10,
         'x': 0,
         'V': 0}]
    grids = []
    x_left = 0
    for gname, scr in zip(cr["struct"], cr["struct_chromos"]):
        grid, border = fun_dict[gname](scr, tube, x_left)
        x_left = border['x']
        grids.append(grid)
        borders.append(border)

    borders[-1]['m'] = m_elem
    solver = {'borders': borders,
              'grids': grids,
              'geom': geom,
              'courant_number': 0.13}
    return solver


def shtraf(p):
    if p < 4000E5:
        return 0
    k = 2e-5
    b = -8e3
    return k * p + b


# VITALY=============================================


# VITALY=============================================

def fitness_foo(wchr):
    try:
        solver = wchr_to_solver(wchr)
        start_time = time.time()
        layers = []
        l_max = solver['geom'][-1][0] + l_max_mini
        for i in range(len(solver['grids'])):
            if solver['grids'][i]['type'] == 'gas':
                layers.append(pn_create_layer(i, solver))
            elif solver['grids'][i]['type'] == 'powder':
                layers.append(ov_create_layer(i, solver))
            elif solver['grids'][i]['type'] == 'piston':
                layers.append(el_pist_create_layer(i, solver))

        Vmax = 0
        pmax = 0

        while True:
            if l_max < layers[-1].x[-1]:
                break
            if (Vmax - layers[-1].V[-1]) > 20:
                break
            if layers[0].time >= time_max:
                return -9999
            if pmax > 7e+8:
                return -9998
            if time.time() - start_time > time_max_calc:
                return -9997
            if (time_v0 < layers[0].time) and (layers[-1].V[-1] < 1E-4):
                return -9995

            layers1 = []
            tau_arr = [l.time_step() for l in layers]
            tau = solver['courant_number'] * min(tau_arr)  # Вычисление шага по времени

            for i in range(len(layers)):
                if len(layers) == 1:
                    layers1.append(layers[i].euler_step(layers, i, tau, flag='one'))
                elif i == 0:
                    layers1.append(layers[i].euler_step(layers, i, tau, flag='left'))
                    # layers[i] = layers1[i]
                elif i != 0 and i != (len(layers) - 1):
                    layers1.append(layers[i].euler_step(layers, i, tau, flag='inter'))
                    # layers[i] = layers1[i]
                elif i == (len(layers) - 1):
                    layers1.append(layers[i].euler_step(layers, i, tau, flag='right'))

            layers = layers1

            if layers[-1].V[-1] > Vmax:
                Vmax = layers[-1].V[-1]
            pmax = max([max(layer.p) for layer in layers])

        fitness = Vmax - shtraf(pmax)
        return fitness
    except Exception as e:
        return -9996


def fit_fun4tst(wchr):
    try:
        solver = wchr_to_solver(wchr)
        vars = [-9999, -9998, -9997, -9996, -9995, rnd.random() * 2000]
        layers = []
        l_max = solver['geom'][-1][0] + l_max_mini
        for i in range(len(solver['grids'])):
            if solver['grids'][i]['type'] == 'gas':
                layers.append(pn_create_layer(i, solver))
            elif solver['grids'][i]['type'] == 'powder':
                layers.append(ov_create_layer(i, solver))
            elif solver['grids'][i]['type'] == 'piston':
                layers.append(el_pist_create_layer(i, solver))

        v = wchr['vec_param'] - 0.5
        res = -sum(v * v)
        time.sleep(1 + rnd.random() * 2)
        # return rnd.choice(vars)+1000
        return res
    except Exception as e:
        return -9990


wchromo4tst = {
    'chromo': {'l_0': 1.03783393439955, 'l_1': 0.2837692999481771, 'l_2': 0.49639150100092677,
               'l_3': 0.14671523779549972, 'l_4': 0.17373486141996183, 'd_0': 0.023,
               'd_1': 0.020563940489849913, 'd_2': 0.017346801263944625, 'd_3': 0.015971082093544947,
               'd_4': 0.015971082093544947, 'struct': ['piston', 'powder', 'gas', 'piston'],
               'struct_chromos': [
                   {'omega_piston': 0.13683056525252102, 'mass_bord_left': 0.08398383125632956,
                    'mass_bord_right': 0.15723713961471766, 'p_bord_force_left': 1317876.6372259278,
                    'p_bord_force_right': 3880795.415544614},
                   {'omega_powder': 0.07902686153317322, 't_powder': 0.0018857243778830708,
                    'powder_type': '12\\7', 'delta_powder': 641.0248292001437,
                    'mass_bord_left': 0.15723713961471766, 'mass_bord_right': 0.15278679017016647,
                    'p_bord_force_left': 3880795.415544614, 'p_bord_force_right': 4500205.209741069},
                   {'p_gas': 49056603.74257242, 'l_gas': 1.1518921399826072, 'type_gas': 'CO_2',
                    'mass_bord_left': 0.15278679017016647, 'mass_bord_right': 0.1448718923802378,
                    'p_bord_force_left': 4500205.209741069, 'p_bord_force_right': 4965965.0644078525},
                   {'omega_piston': 0.04048987485394602, 'mass_bord_left': 0.1448718923802378,
                    'mass_bord_right': 0.2918304748525285, 'p_bord_force_left': 4965965.0644078525,
                    'p_bord_force_right': 1496070.552398767}]}, 'dop_info': None, 'fitness': None, 'id': 0,
    'name': 'Горгидиодо', 'vec_param': np.array([0.49359681, 0.56667194, 0.99276854, 0.2920145, 0.34616205,
                                                 1., 0.81261081, 0.56513856, 0.45931401, 0.45931401,
                                                 0.89869666, 0.25511666, 0.50771427, 0.12301784, 0.38189853,
                                                 0.76696513, 0.94286219, 0.35256207, 0.50771427, 0.49236824,
                                                 0.38189853, 0.44446517, 0.98109426, 0.56507289, 0.49236824,
                                                 0.46507549, 0.44446517, 0.49151162, 0.15761442, 0.46507549,
                                                 0.97182922, 0.49151162, 0.14101723]),
    'vec_struct': np.array([4., 1., 0., 2., 1., 15., 2.])}

print(f'done + sleep {fit_fun4tst(wchromo4tst)} {shtraf(4e8)} {shtraf(5e8)}')
# if __name__ == "__main__":
#     main_tube_w_tst()


