from Pneumatics.PnLayer import PnLayer
from Pneumatics.PnBorder import *


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


def pn_create_layer(solver, tube):
    return PnLayer(solver, func_in, get_flux_left, get_flux_right, tube)
