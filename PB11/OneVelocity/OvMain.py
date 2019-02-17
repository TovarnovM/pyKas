"""
Тестовая задача односкоростной газопороховой смеси.
Сравнить с термодинамикой Никита Валерьевича.
"""


from OneVelocity.OvLayer import OvLayer
from Invariants.Plot import plot_one
from OneVelocity.OvAUSMplus import get_flux_bord


def func_in(x):
    """
    Функция инициализации ячейки
    :param x: координата центра ячейки
    :return: значения параметров плотности, скорости и давления в ячейке
    """
    if x:
        ro = 913.9345
        u = 0
        p = 0.1e6
        z = 0
        return ro, u, p, z


def get_x_v_left(l, tau, p=0):
    """
    Функция получения координаты и скорости левой границы
    :param p: давление с другой стороны границы
    :param l: слой
    :param tau: шаг по времени
    :return: координату и скорость левой границы
    """
    V_l = 0
    x_l = l.x[0] + tau * l.V[0]
    return x_l, V_l


def get_x_v_right(l, tau, p=0):
    """
    Функция получения координаты и скорости правой границы
    :param p: давление с другой стороны границы
    :param l: слой
    :param tau: шаг по времени
    :return: координату и скорость правой границы
    """
    V_r = l.V[-1] + l.get_a(-1, p) * tau
    x_r = l.x[-1] + tau * l.V[-1] + 0.5 * l.get_a(-1, p) * tau * tau
    return x_r, V_r


def func_bord_left(r1, u1, p1, z1, V):
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


def func_bord_right(r1, u1, p1, z1, V):
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
    gps = {'коволюм': 0.0010838,
           'показатель адиабаты': 1.27,
           'координата слева': 0,
           'координата справа': 0.762,
           'скорость лев гран': 0,
           'скорость прав гран': 0,
           'nu': 0.9,
           'давление форсирования': 13.79e6}

    geom = {'координата слева': 0,
            'координата справа': 5.08,
            'диаметр слева': 0.132,
            'диаметр справа': 0.132,
            'количество ячеек': 200}

    gas = {'гпс': gps}

    a = {'I_k': 0.312799,
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
    m = 45.359  # Масса снаряда
    eks = 0.001  # Эксцентреситет - отклонение координаты поршня от длины ствола

    # Создание слоя
    layer = OvLayer(gas, 'гпс', geom, dict_powder, 'порох', get_x_v_left, get_x_v_right, get_flux_left,
                    get_flux_right, func_in, m)

    time_arr = [0]  # Список времени для графика
    V_arr = [0]  # Список скорости для графика
    x_arr = [layer.x[-1]]  # Список координаты для графика
    p_arr = [layer.p[-1]]   # Список давления на дно снаряда для графика
    z_arr = [0]     # Список z для графика
    psi_arr = [layer.powd.psi(layer.z[-1])]     # Список psi для графика

    while (geom['координата справа'] - layer.x[-1]) >= eks:
        tau = r * layer.time_step()  # Вычисление шага по времени
        layer1 = layer.euler_step(layer, tau, 0, 101325)  # Шаг по времени
        layer = layer1
        time_arr.append(layer.time)  # Добавление текущего шага по времени в список для графика
        V_arr.append(layer.V[-1])  # Добавление текущей скорости поршня в список для графика
        x_arr.append(layer.x[-1])  # Добавление текущей координаты поршня в список для графика
        p_arr.append(layer.p[-1])   # Добавление текущего давления на дно снаряда я в список для графика
        z_arr.append(layer.q[3][-1] / layer.q[0][-1])   # Добавление текущего z в список для графика
        psi_arr.append(layer.powd.psi(layer.q[3][-1] / layer.q[0][-1]))  # Добавление текущего psi в список для графика

    print('Время вылета:', time_arr[-1], 'с')
    print('Скорость вылета:', V_arr[-1], 'м/с')

    plot_one(time_arr, V_arr, 'График скорости поршня от времени', 'Время', 'Скорость')
    plot_one(time_arr, p_arr, 'График давления на дно снаряда от времени', 'Время', 'Давление')
    plot_one(time_arr, x_arr, 'График координаты поршня от времени', 'Время', 'Координата')
    plot_one(time_arr, z_arr, 'График z от времени', 'Время', 'z')
    plot_one(time_arr, psi_arr, 'График psi от времени', 'Время', 'psi')


if __name__ == '__main__':
    main()
