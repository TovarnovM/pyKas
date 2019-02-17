"""
Тестовая задача для сравнения с аналитическим решением.
"""


from Pneumatics.PnLayer import PnLayer
import Invariants.Plot as plot
from Pneumatics.PnAUSMplus import get_flux_bord


def func_in(x):
    """
    Функция инициализации ячейки
    :param x: координата центра ячейки
    :return: значения параметров плотности, скорости и давления в ячейке
    """
    if x:
        ro = 40e6 / (287 * 300)
        u = 0
        p = 40e6
        return ro, u, p


def get_x_v_left(l, tau, p_left):
    """
    Функция получения координаты и скорости левой границы
    :param l: слой
    :param tau: шаг по времени
    :return: координату и скорость левой границы
    """
    V_l = l.V[0]
    x_l = l.x[0] + tau * l.V[0]
    return x_l, V_l


def get_x_v_right(l, tau, p_right):
    """
    Функция получения координаты и скорости правой границы
    :param l: слой
    :param tau: шаг по времени
    :return: координату и скорость правой границы
    """
    V_r = l.V[-1] + l.get_a(-1, p_right) * tau
    x_r = l.x[-1] + tau * l.V[-1] + 0.5 * l.get_a(-1, p_right) * tau * tau
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


def get_Vanal(l, c0, S, p0):
    """
    Аналитическая функция получения скорости
    :param l: слой
    :param c0: скорость звука в начальный момент времени
    :param S: площадь
    :param p0: давление в начальных момент времени
    :return:
    """
    return ((2 * c0) / l.const.g[9]) * \
           (1 - ((S * p0 * l.time * (l.const.g[0] + 1)) / (2 * c0 * l.m) + 1) ** (-l.const.g[6]))


def main():
    # Дано:
    air1 = {'коволюм': 0,
            'показатель адиабаты': 1.4,
            'координата слева': 0,
            'координата справа': 0.5,
            'скорость лев гран': 0,
            'скорость прав гран': 0}

    geom = {'координата слева': 0,
            'координата справа': 2,
            'диаметр слева': 0.025,
            'диаметр справа': 0.025,
            'количество ячеек': 100}

    gas = {'воздух1': air1}

    r = 0.1  # Число Куранта
    m = 0.5  # Масса снаряда
    eks = 0.001  # Эксцентреситет - отклонение координаты поршня от длины ствола

    # Создание слоя
    layer = PnLayer(gas, 'воздух1', geom, get_x_v_left, get_x_v_right, get_flux_left, get_flux_right, func_in, m)

    c0 = layer.get_Csound(layer.ro[0], layer.p[0])
    p0 = layer.p[0]
    S = layer.S[0]

    time_arr = [0]  # Список времени для графика
    V_arr = [0]  # Список скорости для графика
    x_arr = [layer.x[-1]]  # Список координаты для графика
    V_arr_anal = [0]
    V_r = [0]

    time_k = (gas['воздух1']['координата справа'] - gas['воздух1']['координата слева']) / c0
    x_v = 0

    while (geom['координата справа'] - layer.x[-1]) >= eks:
        tau = r * layer.time_step()  # Вычисление шага по времени
        layer1 = layer.euler_step(layer, tau, 0)  # Шаг по времени
        layer = layer1

        time_arr.append(layer.time)  # Добавление текущего шага по времени в список для графика
        V_arr.append(layer.V[-1])  # Добавление текущей скорости поршня в список для графика
        x_arr.append(layer.x[-1])  # Добавление текущей координаты поршня в список для графика

        if layer.time >= time_k:
            x_v += c0 * tau
        if x_v <= layer.x[-1]:
            anal = get_Vanal(layer, c0, S, p0)
            V_arr_anal.append(anal)
            V_r.append((layer.V[-1] - anal) / abs(anal) * 100)
        else:
            V_arr_anal.append(V_arr_anal[-1])
            V_r.append(V_r[-1])

    plot.plot_two(time_arr, V_arr, V_arr_anal, 'Сравнение численного и аналитического решения',
                  'Численный метод', 'Аналитический метод')
    plot.plot_one(time_arr, V_r, 'Погрешность', 'Время', 'Проценты')


if __name__ == '__main__':
    main()
