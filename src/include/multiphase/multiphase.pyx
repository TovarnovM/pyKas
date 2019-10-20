# distutils: language=c++
# cython: language_level=3, boundscheck=False, nonecheck=False, cdivision=True, initializedcheck=False

from gaslayer cimport GasLayer, GasEOS, GasFluxCalculator, GridStrecher
from tube cimport Tube, InterpXY

import cython

import numpy as np
cimport numpy as np



cdef class GasPhase(GasLayer):
    def __init__(self, n_cells, Tube tube, GasEOS gasEOS, GasFluxCalculator flux_calculator, GridStrecher grid_strecher):
        super().__init__(n_cells, tube, gasEOS, flux_calculator, grid_strecher, 4)
        self.alphas = np.ones(n_cells, dtype=np.double)
        
    cpdef void copy_params_to_gas_phase(self, GasPhase to_me):
        self.copy_params_to(to_me)
        to_me.alphas[:] = self.alphas

    cpdef GasLayer copy(self):
        cdef GasPhase res = GasPhase(self.n_cells, self.tube, self.gasEOS, self.flux_calculator, self.grid_strecher)
        self.copy_params_to_gas_phase(res)
        return res
    
    


cdef class SolidPhase(GasPhase):
    pass

class MultiPhaseLayer(object):
    @classmethod
    def from_dict(cls, tube: Tube, x_1: float, layer_dict, calc_settings):
        """основной метод для создания полидисперсного слоя

        
        Arguments:
            layer_dict {dict} -- словарь с начальными данными
            {
                'init_conds':{
                    'p_0': 100e5,   # начальное давление газа
                    'x_right': 0.7, # x_left, W_0
                }
                'phases': {
                    'phase1':{
                        'type': 'gas',
                        'name': 'air',
                        'gamma': 1.4,
                        'kappa': 0.0010838,
                        'R': 287,
                        'T_0': 300,     # K
                        'W_0': 0.0002,  # м^3 объем газа
                        'u_0': 0,       # начальная скорость
                    },
                    'phase2':{
                        'type': 'gas',
                        'name': 'powdergas',
                        'gamma': 1.228,
                        'kappa': 0.0010838,
                        'R': 287,
                        'T_0': 300,     # K
                        'W_0': 0.0002,  # м^3 объем газа
                        'u_0': 0,       # начальная скорость
                    },
                    'phase3':{
                        'type': 'solid',
                        'name': 'powder',
                        'powder': {
                            'name': '4\\7',
                            'f': 1.027,
                            'etta': 0.228,
                            'alpha_k': 1.008,
                            'T_1': 3006.0,
                            'ro': 1.6,
                            'I_k': 0.32,
                            'Z_k': 1.488,
                            'k_1': 0.811,
                            'lambda_1': 0.081,
                            'k_2': 0.505,
                            'lambda_2': -1.024,
                            'k_f': 0.0003,
                            'k_l': 0.0016
                        },
                        'omega': 35,  # грамм
                        'delta': 700, # г/cm^3
                        't_ign': 0.00001, # начало горения
                        'u_0': 0,     #начальная скорость
                    }

                },
                'connections':[
                    ('phase3', 'phase2')
                ]
            }

            'calc_settings': {
                'cell_dx': 0.0025,
                'n_cells_min': 13,
                'n_cells_max': 300,
                'GridStrecher_kwargs': {}
            }
        """
        def get_phase_layer(d):
            if d['type'] == 'gas':
                return GasPhase.get_standart(tube, x_1, d, calc_settings)
            if d['type'] == 'solid':
                return SolidPhase.get_standart(tube, x_1, d, calc_settings)
            raise Exception(f'Неправильный тип у фазы {d}')
        layers = {name: get_phase_layer(d) for name, d in layer_dict['phases'].items()}



    def __init__(self, layers):
        self.layers = layers

    def copy(self):
        res = MultiPhaseLayer({name:lr.copy() for name, lr in self.layers.items()})
        return res