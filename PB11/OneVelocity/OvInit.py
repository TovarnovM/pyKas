from OneVelocity.OvLayer import OvLayer
from OneVelocity.OvBorder import *


def func_in(x, grid, W_km):
    """
    Функция инициализации ячейки
    :param x: координата центра ячейки
    :return: значения параметров плотности, скорости и давления в ячейке
    """
    if x:
        # ro = grid['init_const']['omega'] / W_km
        ro = grid['init_const']['ro']
        u = 0
        p = 101325
        z = 0
        return ro, u, p, z


def ov_create_layer(solver_grid, tube):
    return OvLayer(solver_grid, func_in, get_flux_left, get_flux_right, tube)
