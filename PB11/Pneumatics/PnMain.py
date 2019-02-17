from Pneumatics.PnLayer import PnLayer
from Invariants.Plot import plot_one
from Pneumatics.PnAUSMplus import get_flux_bord


def func_in(x):
    """
    Функция инициализации ячейки
    :param x: координата центра ячейки
    :return: значения параметров плотности, скорости и давления в ячейке
    """
    if x:
        ro = 40e6 / ((287 * 300) + 0.0010838 * 40e6)
        # ro = 40e6 / (287 * 300)
        u = 0
        p = 40e6
        return ro, u, p


def get_x_v_left(l, tau, p_left):
    """
    Функция получения координаты и скорости левой границы
    :param p_left: давление слева от левой границы
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
    :param p_right: давление справа от правой границы
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


def main():
    # Дано:
    # коволюм 0.0010838
    air1 = {'коволюм': 0.0010838,
            'показатель адиабаты': 1.4,
            'координата слева': 0,
            'координата справа': 0.5,
            'скорость лев гран': 0,
            'скорость прав гран': 0}

    geom = {'координата слева': 0,
            'координата справа': 2,
            'диаметр слева': 0.025,
            'диаметр справа': 0.025,
            'количество ячеек': 200}

    gas = {'воздух1': air1}

    r = 0.5                         # Число Куранта
    m = 0.05                        # Масса снаряда
    eks = 0.001                     # Эксцентреситет - отклонение координаты поршня от длины ствола
    p_atm = 101325

    # Создание слоя
    layer = PnLayer(gas, 'воздух1', geom, get_x_v_left, get_x_v_right, get_flux_left, get_flux_right, func_in, m)

    time_arr = [0]              # Список времени для графика
    V_arr = [0]                 # Список скорости для графика
    x_arr = [layer.x[-1]]           # Список координаты для графика

    while (geom['координата справа'] - layer.x[-1]) >= eks:
        tau = r * layer.time_step()     # Вычисление шага по времени
        layer1 = layer.euler_step(layer, tau, p_right=p_atm)    # Шаг по времени
        layer = layer1
        time_arr.append(layer.time)         # Добавление текущего шага по времени в список для графика
        V_arr.append(layer.V[-1])           # Добавление текущей скорости поршня в список для графика
        x_arr.append(layer.x[-1])           # Добавление текущей координаты поршня в список для графика

    print('Время вылета:', time_arr[-1], 'с')
    print('Скорость вылета:', V_arr[-1], 'м/с')

    plot_one(time_arr, V_arr, 'График скорости поршня от времени', 'Время', 'Скорость')
    plot_one(x_arr, V_arr, 'График скорости поршня от координаты', 'Координата', 'Скорость')
    plot_one(time_arr, x_arr, 'График координаты поршня от времени', 'Время', 'Координата')


if __name__ == '__main__':
    main()
