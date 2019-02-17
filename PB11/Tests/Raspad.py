"""
Тестовая задача о распаде разрыва.
Сравнить графики с Зализняком (глава 6).
"""


from Invariants.Plot import plot_one
from Pneumatics.PnLayer import PnLayer
from Pneumatics.PnAUSMplus import get_flux_bord
import numpy as np


def func_in(x):
    if x < 0:
        ro = 1
        u = 0
        p = 1
    else:
        ro = 0.125
        u = 0
        p = 0.1
    return ro, u, p


def get_x_v_left(l, tau, p_left):
    """
    Функция получения координаты и скорости левой границы
    :param l: слой
    :param tau: шаг по времени
    :return: координату и скорость левой границы
    """
    V_l = l.V[0]
    x_l = l.x[0] + tau * V_l
    return x_l, V_l


def get_x_v_right(l, tau, p_right):
    """
    Функция получения координаты и скорости правой границы
    :param l: слой
    :param tau: шаг по времени
    :return: координату и скорость правой границы
    """
    V_r = l.V[-1]
    x_r = l.x[-1] + tau * V_r
    return x_r, V_r


def func_bord_left(r1, u1, p1, V):
    """

    :param r1:
    :param u1:
    :param p1:
    :param V:
    :return:
    """
    r2 = r1
    u2 = -u1 + 2 * V
    p2 = p1
    return r2, u2, p2


def func_bord_right(r1, u1, p1, V):
    """

    :param r1:
    :param u1:
    :param p1:
    :param V:
    :return:
    """
    r2 = r1
    u2 = -u1 + 2 * V
    p2 = p1
    return r2, u2, p2


def get_flux_left(l):
    """

    :param l: слой
    :return: Фукнцию для определения потока через левую границу сетки
    """
    return get_flux_bord(l, 0, func_bord_left)


def get_flux_right(l):
    """

    :param l: слой
    :return: Фукнцию для определения потока через правую границу сетки
    """
    return get_flux_bord(l, -1, func_bord_right)


def main():
    # Дано:
    air1 = {'коволюм': 0,
            'показатель адиабаты': 1.4,
            'координата слева': -1,
            'координата справа': 1,
            'скорость лев гран': 0,
            'скорость прав гран': 0}

    geom = {'координата слева': -10,
            'координата справа': 10,
            'диаметр слева': 0.1,
            'диаметр справа': 0.1,
            'количество ячеек': 100}

    gas = {'воздух1': air1}

    r = 0.1  # Число Куранта
    m = 0

    # Создание ячейки и слоя
    layer = PnLayer(gas, 'воздух1', geom, get_x_v_left, get_x_v_right, get_flux_left, get_flux_right, func_in, m)

    time_arr = [0]  # Список времени для графика

    while time_arr[-1] <= 0.5:
        tau = r * layer.time_step()  # Вычисление шага по времени
        layer1 = layer.euler_step(layer, tau)  # Шаг по времени
        layer = layer1
        time_arr.append(layer.time)  # Добавление текущего шага по времени в список для графика

    xr = np.roll(layer.x, -1)
    xc = ((xr + layer.x) / 2)[:-1]

    plot_one(xc, layer.q[0], 'Распределение плотности', 'Координата', 'Плотность')
    plot_one(xc, layer.q[1] / layer.q[0], 'Распределение скорости', 'Координата', 'Скорость')
    plot_one(xc, layer.p, 'Распределение давления', 'Координата', 'Давление')


if __name__ == '__main__':
    main()
