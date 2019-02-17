from OneVelocity.OvLayer import OvLayer
from OneVelocity.OvAUSMplus import get_flux_bord as flux_bord_powd
from Pneumatics.PnLayer import PnLayer
from Pneumatics.PnAUSMplus import get_flux_bord as flux_bord_gas
from Invariants.Plot import plot_one


def func_in_powd(x):
    """
    Функция инициализации ячейки
    :param x: координата центра ячейки
    :return: значения параметров плотности, скорости и давления в ячейке
    """
    if x:
        ro = 900
        u = 0
        p = 1e6
        z = 0
        return ro, u, p, z


def func_in_gas(x):
    """
    Функция инициализации ячейки
    :param x: координата центра ячейки
    :return: значения параметров плотности, скорости и давления в ячейке
    """
    if x:
        ro = 700
        u = 0
        p = 2.5e5
        return ro, u, p


def get_x_v_left_powd(l, tau, p_left):
    """
    Функция получения координаты и скорости левой границы
    :param l: слой
    :param tau: шаг по времени
    :return: координату и скорость левой границы
    """
    V_l = 0
    x_l = l.x[0] + tau * l.V[0]
    return x_l, V_l


def get_x_v_right_powd(l, tau, p_right):
    """
    Функция получения координаты и скорости правой границы
    :param p_right: давление справа от границы
    :param l: слой
    :param tau: шаг по времени
    :return: координату и скорость правой границы
    """
    V_r = l.V[-1] + l.get_a(-1, p_right) * tau
    x_r = l.x[-1] + tau * l.V[-1] + 0.5 * l.get_a(-1, p_right) * tau * tau
    return x_r, V_r


def get_x_v_left_gas(l, tau, p_left):
    """
    Функция получения координаты и скорости левой границы
    :param p_left: давление слева от границы
    :param l: слой
    :param tau: шаг по времени
    :return: координату и скорость левой границы
    """
    V_l = l.V[-1] + l.get_a(-1, p_left) * tau
    x_l = l.x[0] + tau * l.V[0] + 0.5 * l.get_a(0, p_left) * tau * tau
    return x_l, V_l


def get_x_v_right_gas(l, tau, p_right):
    """
    Функция получения координаты и скорости правой границы
    :param p_right: давление справа от границы
    :param l: слой
    :param tau: шаг по времени
    :return: координату и скорость правой границы
    """
    V_r = l.V[-1] + l.get_a(-1, p_right) * tau
    x_r = l.x[-1] + tau * l.V[-1] + 0.5 * l.get_a(-1, p_right) * tau * tau
    return x_r, V_r


def func_bord_left_powd(r1, u1, p1, z1, V):
    """

    :param z1:
    :param r1:
    :param u1:
    :param p1:
    :param V:
    :return:
    """
    r2 = r1
    u2 = -u1 + 2 * V
    p2 = p1
    z2 = z1
    return r2, u2, p2, z2


def func_bord_right_powd(r1, u1, p1, z1, V):
    """

    :param z1:
    :param r1:
    :param u1:
    :param p1:
    :param V:
    :return:
    """
    r2 = r1
    u2 = -u1 + 2 * V
    p2 = p1
    z2 = z1
    return r2, u2, p2, z2


def get_flux_left_powd(l):
    """

    :param l: слой
    :return: Фукнцию для определения потока через левую границу сетки
    """
    return flux_bord_powd(l, 0, func_bord_left_powd)


def get_flux_right_powd(l):
    """

    :param l: слой
    :return: Фукнцию для определения потока через правую границу сетки
    """
    return flux_bord_powd(l, -1, func_bord_right_powd)


def func_bord_left_gas(r1, u1, p1, V):
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


def func_bord_right_gas(r1, u1, p1, V):
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


def get_flux_left_gas(l):
    """

    :param l: слой
    :return: Фукнцию для определения потока через левую границу сетки
    """
    return flux_bord_gas(l, 0, func_bord_left_gas)


def get_flux_right_gas(l):
    """

    :param l: слой
    :return: Фукнцию для определения потока через правую границу сетки
    """
    return flux_bord_gas(l, -1, func_bord_right_gas)


def main():
    # Дано:
    powd = {'коволюм': 0.0010838,
            'показатель адиабаты': 1.4,
            'координата слева': 0,
            'координата справа': 0.25,
            'скорость лев гран': 0,
            'скорость прав гран': 0,
            'nu': 1,
            'давление форсирования': 0}

    air = {'коволюм': 0.0010838,
           'показатель адиабаты': 1.4,
           'координата слева': 0.25,
           'координата справа': 0.5,
           'скорость лев гран': 0,
           'скорость прав гран': 0,
           'nu': 1}

    geom = {'координата слева': 0,
            'координата справа': 2,
            'диаметр слева': 0.025,
            'диаметр справа': 0.025,
            'количество ячеек': 100}

    gas = {'гпс': powd,
           'газ': air}

    a = {'I_k': 1.257,
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

    dict_powder = {'порох': a}

    r = 0.3  # Число Куранта
    m = 1  # Масса снаряда
    m1 = 1
    eks = 0.001  # Эксцентреситет - отклонение координаты поршня от длины ствола

    # Создание слоя
    layer_powd = OvLayer(gas, 'гпс', geom, dict_powder, 'порох', get_x_v_left_powd, get_x_v_right_powd,
                         get_flux_left_powd, get_flux_right_powd, func_in_powd, m1)

    layer_gas = PnLayer(gas, 'газ', geom, get_x_v_left_gas, get_x_v_right_gas, get_flux_left_gas, get_flux_right_gas,
                        func_in_gas, m)

    time_arr = [0]  # Список времени для графика
    V_arr = [0]  # Список скорости для графика
    x_arr = [layer_gas.x[-1]]  # Список координаты для графика

    while (geom['координата справа'] - layer_gas.x[-1]) >= eks:
        tau = r * min(layer_powd.time_step(), layer_gas.time_step())  # Вычисление шага по времени

        layer_powd1 = layer_powd.euler_step(layer_powd, tau, 0, layer_gas.p[0])  # Шаг по времени
        layer_powd = layer_powd1

        layer_gas1 = layer_gas.euler_step(layer_gas, tau, layer_powd.p[-1], 101325)  # Шаг по времени
        layer_gas = layer_gas1

        time_arr.append(layer_powd.time)  # Добавление текущего шага по времени в список для графика
        V_arr.append(layer_powd.V[-1])  # Добавление текущей скорости поршня в список для графика
        x_arr.append(layer_powd.x[-1])  # Добавление текущей координаты поршня в список для графика

    print('Время вылета:', time_arr[-1], 'с')
    print('Скорость вылета:', V_arr[-1], 'м/с')

    plot_one(time_arr, V_arr, 'График скорости поршня от времени', 'Время', 'Скорость')
    plot_one(x_arr, V_arr, 'График скорости поршня от координаты', 'Координата', 'Скорость')
    plot_one(time_arr, x_arr, 'График координаты поршня от времени', 'Время', 'Координата')


if __name__ == '__main__':
    main()
