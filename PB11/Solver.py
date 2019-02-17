# Чистая пневматика

# border0 = {'m': 100,
#            'p_f': 1e8,
#            'x': 0,
#            'V': 0}
#
# border1 = {'m': 0.05,
#            'p_f': 101325,
#            'x': 0.5,
#            'V': 0}
#
# const0 = {'covolume': 0.0010838,
#           'R': 287,
#           'gamma': 1.4}
#
# init_const0 = {'p': 40e6,
#                'T': 300}
#
# grid0 = {'name': 'gas',
#          'n_cells': 200,
#          'consts': const0,
#          'type': 'gas',
#          'init_const': init_const0}
#
#
# solver = {'borders': [border0, border1],
#           'grids': [grid0],
#           'geom': [(0, 0.025), (2, 0.025)],
#           'courant_number': 0.5}


# # Чистый порох
#
# border0 = {'m': 100,
#            'p_f': 1e8,
#            'x': 0,
#            'V': 0}
#
# border1 = {'m': 45.359,
#            'p_f': 13.79e6,
#            'x': 0.762,
#            'V': 0}
#
# powder = {'I_k': 0.312799,
#           'T_1': 2795.0,
#           'Z_k': 1.536,
#           'alpha_k': 1.0838,
#           'etta': 0.231,
#           'f': 1.009,
#           'k_1': 0.7185,
#           'k_2': 0.5386,
#           'k_f': 0.0003,
#           'k_l': 0.0016,
#           'lambda_1': 0.2049,
#           'lambda_2': -0.8977,
#           'name': 'порох',
#           'ro': 1.575}
#
# const0 = {'covolume': 0.0010838,
#           'R': 287,
#           'gamma': 1.27,
#           'param_powder': powder,
#           'nu': 0.9}
#
# init_const0 = {'omega': 9.5255}
#
# grid0 = {'name': 'powder',
#          'n_cells': 200,
#          'consts': const0,
#          'type': 'powder',
#          'init_const': init_const0}
#
# solver = {'borders': [border0, border1],
#           'grids': [grid0],
#           'geom': [(0, 0.132), (5.08, 0.132)],
#           'courant_number': 0.3}


# Газ - порох

# border0 = {'m': 100,
#            'p_f': 1e8,
#            'x': 0,
#            'V': 0}
#
# border1 = {'m': 45,
#            'p_f': 13.79e6,
#            'x': 0.762,
#            'V': 0}
#
# border2 = {'m': 1,
#            'p_f': 101325,
#            'x': 1,
#            'V': 0}
#
# powder = {'I_k': 0.312799,
#           'T_1': 2795.0,
#           'Z_k': 1.536,
#           'alpha_k': 1.0838,
#           'etta': 0.231,
#           'f': 1.009,
#           'k_1': 0.7185,
#           'k_2': 0.5386,
#           'k_f': 0.0003,
#           'k_l': 0.0016,
#           'lambda_1': 0.2049,
#           'lambda_2': -0.8977,
#           'name': 'порох',
#           'ro': 1.575}
#
# const0 = {'covolume': 0.0010838,
#           'R': 287,
#           'gamma': 1.27,
#           'param_powder': powder,
#           'nu': 0.9}
#
# const1 = {'covolume': 0.0010838,
#           'R': 287,
#           'gamma': 1.4}
#
# init_const0 = {'omega': 9.5255}
#
# init_const1 = {'p': 40e6,
#                'T': 300}
#
# grid0 = {'name': 'powder',
#          'n_cells': 100,
#          'consts': const0,
#          'type': 'powder',
#          'init_const': init_const0}
#
# grid1 = {'name': 'gas',
#          'n_cells': 100,
#          'consts': const1,
#          'type': 'gas',
#          'init_const': init_const1}
#
# solver = {'borders': [border0, border1, border2],
#           'grids': [grid0, grid1],
#           'geom': [(0, 0.132), (5.08, 0.132)],
#           'courant_number': 0.5}

# Газ - эластичный пистон

# border0 = {'m': 100,
#            'p_f': 1e8,
#            'x': 0,
#            'V': 0}
#
# border1 = {'m': 0.1,
#            'p_f': 101325,
#            'x': 0.5,
#            'V': 0}
#
# border2 = {'m': 0.1,
#            'p_f': 101325,
#            'x': 0.6,
#            'V': 0}
#
# const0 = {'covolume': 0.0010838,
#           'R': 287,
#           'gamma': 1.4}
#
# const1 = {'covolume': 0,
#           'gamma': 1.63098,
#           'c0': 2380,
#           'ro0': 919.03,
#           'sigmas': 25.2e6,
#           'k0': 0.054,
#           'b1': 0.027,
#           'b2': 0.00675,
#           'mu': 0.001,
#           'taus': 1e6}
#
# init_const0 = {'p': 40e6,
#                'T': 300}
#
# grid0 = {'name': 'gas',
#          'n_cells': 100,
#          'consts': const0,
#          'type': 'gas',
#          'init_const': init_const0}
#
# grid1 = {'name': 'piston',
#          'n_cells': 50,
#          'consts': const1,
#          'type': 'piston'}
#
#
# solver = {'borders': [border0, border1, border2],
#           'grids': [grid0, grid1],
#           'geom': [(0, 0.025), (1, 0.025)],
#           'courant_number': 0.4}

# 'geom': [(0, 0.025), (1.2, 0.025), (1.5, 0.02), (2, 0.02)],

# Порох - эластичный пистон

# border0 = {'m': 10000,
#            'p_f': 1e8,
#            'x': 0,
#            'V': 0}
#
# border1 = {'m': 0.1,
#            'p_f': 101325,
#            'x': 0.2,
#            'V': 0}
#
# border2 = {'m': 0.005,
#            'p_f': 101325,
#            'x': 0.27,
#            'V': 0}
#
# powder = {'I_k': 0.32,
#           'T_1': 3006.0,
#           'Z_k': 1.488,
#           'alpha_k': 1.008,
#           'etta': 0.228,
#           'f': 1.027,
#           'k_1': 0.811,
#           'k_2': 0.505,
#           'k_f': 0.0003,
#           'k_l': 0.0016,
#           'lambda_1': 0.081,
#           'lambda_2': -1.024,
#           'name': 'порох',
#           'ro': 1.6}
#
# const0 = {'covolume': 0,
#           'R': 287,
#           'gamma': 1.228,
#           'param_powder': powder,
#           'nu': 1}
#
# const1 = {'covolume': 0,
#           'gamma': 1.63098,
#           'c0': 2380,
#           'ro0': 919.03,
#           'sigmas': 25.2e6,
#           'k0': 0.054,
#           'b1': 0.027,
#           'b2': 0.00675,
#           'mu': 0.001,
#           'taus': 1e6}
#
# init_const0 = {'ro': 840,
#                't_init': 0.0001}
#
# grid0 = {'name': 'powder',
#          'n_cells': 50,
#          'consts': const0,
#          'type': 'powder',
#          'init_const': init_const0}
#
# grid1 = {'name': 'piston',
#          'n_cells': 50,
#          'consts': const1,
#          'type': 'piston'}
#
#
# solver = {'borders': [border0, border1, border2],
#           'grids': [grid0, grid1],
#           'geom': [(0, 0.023), (1.149, 0.023), (1.229, 0.014), (1.5, 0.014)],
#           'courant_number': 0.13}

# Порох - поршень - газ

# border0 = {'m': 100,
#            'p_f': 1e8,
#            'x': 0,
#            'V': 0}
#
# border1 = {'m': 0.01,
#            'p_f': 101325,
#            'x': 0.2,
#            'V': 0}
#
# border2 = {'m': 0.1,
#            'p_f': 101325,
#            'x': 0.3,
#            'V': 0}
#
# border3 = {'m': 0.01,
#            'p_f': 1e8,
#            'x': 1,
#            'V': 0}
#
# powder = {'I_k': 0.3,
#           'T_1': 2736.0,
#           'Z_k': 1.602,
#           'alpha_k': 1.053,
#           'etta': 0.236,
#           'f': 0.988,
#           'k_1': 0.653,
#           'k_2': 0.65,
#           'k_f': 0.0003,
#           'k_l': 0.0016,
#           'lambda_1': 0.247,
#           'lambda_2': -0.791,
#           'name': 'порох',
#           'ro': 1.5}
#
# const0 = {'covolume': 0.0010838,
#           'R': 287,
#           'gamma': 1.27,
#           'param_powder': powder,
#           'nu': 1}
#
# const1 = {'covolume': 0,
#           'gamma': 1.63098,
#           'c0': 2380,
#           'ro0': 919.03,
#           'sigmas': 25.2e6,
#           'k0': 0.054,
#           'b1': 0.027,
#           'b2': 0.00675,
#           'mu': 0.001,
#           'taus': 1e6}
#
# const2 = {'covolume': 0.0010838,
#           'R': 287,
#           'gamma': 1.4}
#
# init_const0 = {'omega': 9.5255}
#
# init_const2 = {'p': 200e6,
#                'ro': 250}
#
# grid0 = {'name': 'powder',
#          'n_cells': 50,
#          'consts': const0,
#          'type': 'powder',
#          'init_const': init_const0}
#
# grid1 = {'name': 'piston',
#          'n_cells': 90,
#          'consts': const1,
#          'type': 'piston'}
#
# grid2 = {'name': 'gas',
#          'n_cells': 50,
#          'consts': const2,
#          'type': 'gas',
#          'init_const': init_const2}
#
# solver = {'borders': [border0, border1, border2, border3],
#           'grids': [grid0, grid1, grid2],
#           'geom': [(0, 0.023), (0.92, 0.023), (1, 0.014), (6.5, 0.014)],
#           'courant_number': 0.3}


border0 = {'m': 0.04252,
           'p_f': 101325,
           'x': 0,
           'V': 600}

border1 = {'m': 0.01756,
           'p_f': 101325,
           'x': 0.063,
           'V': 600}

const0 = {'covolume': 0,
          'gamma': 1.63098,
          'c0': 2380,
          'ro0': 919.03,
          'sigmas': 25.2e6,
          'k0': 0.054,
          'b1': 0.027,
          'b2': 0.00675,
          'mu': 0.001,
          'taus': 1e6}


grid0 = {'name': 'piston',
         'n_cells': 100,
         'consts': const0,
         'type': 'piston'}


solver = {'borders': [border0, border1],
          'grids': [grid0],
          'geom': [(0, 0.023), (0.08, 0.023), (0.2, 0.018), (0.3, 0.018)],
          'courant_number': 0.1}
