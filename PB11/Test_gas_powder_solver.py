# Чистый порох

border0 = {'m': 1000000,
           'p_f': 1e8,
           'x': 0,
           'V': 0}

border1 = {'m': 0.015,
           'p_f': 20e6,
           'x': 0.3,
           'V': 0}

border2 = {'m': 0.036,
           'p_f': 100e6,
           'x': 1.8,
           'V': 0}

powder = {'I_k': 0.3,
          'T_1': 2736.0,
          'Z_k': 1.602,
          'alpha_k': 1.053,
          'etta': 0.24,
          'f': 0.988,
          'k_1': 0.653,
          'k_2': 0.650,
          'k_f': 0.0003,
          'k_l': 0.0016,
          'lambda_1': 0.247,
          'lambda_2': -0.791,
          'name': 'порох',
          'ro': 1.6}

const0 = {'covolume': 0.0010838,
          'R': 340,
          'gamma': 1.24,
          'param_powder': powder,
          'nu': 1}

init_const0 = {'ro': 900, 't_init': 0.00}

const1 = {'covolume': 0.0010838,
          'R': 287,
          'gamma': 1.4}

init_const1 = {'p': 20e6,
               'T': 300}

grid0 = {'name': 'powder',
         'n_cells': 100,
         'consts': const0,
         'type': 'powder',
         'init_const': init_const0}

grid1 = {'name': 'gas',
         'n_cells': 100,
         'consts': const1,
         'type': 'gas',
         'init_const': init_const1}

solver = {'borders': [border0, border1, border2],
          'grids': [grid0, grid1],
          'geom': [(0, 0.035), (0.27, 0.035), (0.3, 0.03), (1.4, 0.03), (1.8, 0.0145), (13.8, 0.0145)],
          'courant_number': 0.1}