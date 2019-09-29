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
    
    


cdef class SolidPhase(GasLayer):
    pass

class MultiPhaseLayer(object):
    @classmethod
    def from_dict(cls, d):
        """основной метод для создания полидисперсного слоя

        
        Arguments:
            d {dict} -- словарь с начальными данными
            {
                'phases': {
                    'phase1':{
                        'type': 'gas',
                        'name': 'air',
                        'gamma': 1.4,
                        'kappa': 0.0010838,
                        'R': 287,
                        'T_0': 300,     # K
                        'W_0': 0.0002,  # м^3 объем газа
                        'p_0': 100e5,   # начальное давление газа
                        'u_0': 0,       # начальная скорость
                    },
                    'phase3':{
                        'type': 'powder',
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
                            'p_0': 5e6, # начальное давление
                            't_ign': 0.00001, # начало горения
                            'u_0': 0,     #начальная скорость
                    }

                },
                'calc_settings': {
                    'cell_dx': 0.0025,
                    'n_cells_min': 13,
                    'n_cells_max': 300,
                    'GridStrecher_kwargs': {}
                }
            }
        """


    def __init__(self, layers):
        self.layers = layers

    def copy(self):
        res = MultiPhaseLayer({name:lr.copy() for name, lr in self.layers.items()})
        return res