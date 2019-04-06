# Чистая пневматика

border0 = {'m': 100,
           'p_f': 1e8,
           'x': 0,
           'V': 0}

border1 = {'m': 0.05,
           'p_f': 1e8,
           'x': 0.5,
           'V': 0}

const0 = {'covolume': 0.0010838,
          'R': 287,
          'gamma': 1.4}

init_const0 = {'p': 40e6,
               'T': 300}

grid0 = {'name': 'gas',
         'n_cells': 100,
         'consts': const0,
         'type': 'gas',
         'init_const': init_const0}


solver = {'borders': [border0, border1],
          'grids': [grid0],
          'geom': [(0, 0.025), (2, 0.025)],
          'courant_number': 0.1}