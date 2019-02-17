# Чистый порох

border0 = {'m': 1000000,
           'p_f': 1e8,
           'x': 0,
           'V': 0}

border1 = {'m': 45.359,
           'p_f': 13.79e6,
           'x': 0.762,
           'V': 0}

powder = {'I_k': 0.312799,
          'T_1': 2795.0,
          'Z_k': 1.536,
          'alpha_k': 1.0838,
          'etta': 0.231,
          'f': 1.009,
          'k_1': 0.7185,
          'k_2': 0.5386,
          'k_f': 0.0003,
          'k_l': 0.0016,
          'lambda_1': 0.2049,
          'lambda_2': -0.8977,
          'name': 'порох',
          'ro': 1.575}

const0 = {'covolume': 0.0010838,
          'R': 340,
          'gamma': 1.27,
          'param_powder': powder,
          'nu': 0.9}

init_const0 = {'ro': 913.93, 't_init': 0.01}

grid0 = {'name': 'powder',
         'n_cells': 100,
         'consts': const0,
         'type': 'powder',
         'init_const': init_const0}

solver = {'borders': [border0, border1],
          'grids': [grid0],
          'geom': [(0, 0.132), (5.08, 0.132)],
          'courant_number': 0.3}