from ElasticPiston.ElPistLayer import ElPistLayer
from ElasticPiston.ElPistBorder import *


def func_in(x, grid):
    """
    Функция инициализации ячейки
    :param x: координата центра ячейки
    :return: значения параметров плотности, скорости и давления в ячейке
    """
    if x:
        # ro = grid['consts']['ro0']
        ro = 921
        u = 0
        p = 101325
        return ro, u, p


def el_pist_create_layer(solver, tube):
    return ElPistLayer(solver, func_in, get_flux_left, get_flux_right, tube)
