import numpy as np
import math
from Pneumatics import PnAUSMplus


def AUSMp(l):
    """
    Схема AUSM плюс для получения потоков f на границах ячеек
    :param l: слой
    :return:  список нумпи массивов f - потоков через границы (не считая левой и правой границы сетки)
    """
    f = PnAUSMplus.AUSMp(l)

    u = l.q[1] / l.q[0]
    u_right = np.append(0.5 * (u + np.roll(u, -1))[:-1], l.V[-1])
    u_left = np.roll(u_right, 1)
    u_left[0] = l.V[0]

    ha = 2 * (u_right - u_left) / l.dx - u * l.tube.get_stuff(l.x) / l.tube.get_S(l.x_c)
    hh = (0.5 / math.sqrt(3)) * np.abs(ha)

    tauxx = np.zeros_like(ha, dtype=np.float64)

    cond1 = np.abs(hh) > 1e-11
    tauxx[cond1] = (2 / 3) * (l.mu + 0.5 * l.taus / hh[cond1]) * ha[cond1]

    f[1] -= tauxx
    f[2] -= tauxx * u_right

    return f


def get_flux_bord(l, i, func_bord):
    """
    Схема AUSM плюс для получения потоков f на границах сетки
    :param l: слой
    :param i: индекс 0 или -1 для левой и правой границы соответственно
    :param func_bord: функция граничных условий
    :return: кортеж потоков f через границы
    """
    f1, f2, f3 = PnAUSMplus.get_flux_bord(l, i, func_bord)

    u = l.q[1][i] / l.q[0][i]

    if i == -1:
        ha = 4 * (l.V[i] - u) / l.dx[i]
    else:
        ha = 4 * (u - l.V[0]) / l.dx[0]

    hh = (0.5 / math.sqrt(3)) * abs(ha)

    if hh != 0:
        tauxx = (2 / 3) * (l.mu + 0.5 * l.taus / hh) * ha
    else:
        tauxx = 0

    f2 -= tauxx
    f3 -= tauxx * l.V[i]

    return f1, f2, f3


def get_f(l):
    """
    Функция для получения конечных списков массивов потоков f
    :param l: слой
    :return: два списка нумпи массивов потоков f через левые и правые границы ячеек, включая границы сетки
    """
    f_right = AUSMp(l)
    f_border_r = l.flux_right(l)
    for i in range(len(f_right)):
        f_right[i][-1] = f_border_r[i]
    f_left = [np.roll(ar, 1) for ar in f_right]
    f_border_l = l.flux_left(l)
    for i in range(len(f_left)):
        f_left[i][0] = f_border_l[i]
    return f_left, f_right
