from ElasticPiston.ElAUSMplus import get_flux_bord


def get_x_v_left(l, tau, p_left):
    """
    Функция получения координаты и скорости левой границы
    :param p_left: давление слева от левой границы
    :param l: слой
    :param tau: шаг по времени
    :return: координату и скорость левой границы
    """
    V_l = round(l.V[0] + l.get_a(0, p_left) * tau, 8)
    x_l = round(l.x[0] + tau * l.V[0] + 0.5 * l.get_a(0, p_left) * tau * tau, 8)
    return x_l, V_l


def get_x_v_right(l, tau, p_right):
    """
    Функция получения координаты и скорости правой границы
    :param p_right: давление справа от правой границы
    :param l: слой
    :param tau: шаг по времени
    :return: координату и скорость правой границы
    """
    V_r = round(l.V[-1] + l.get_a(-1, p_right) * tau, 8)
    x_r = round(l.x[-1] + tau * l.V[-1] + 0.5 * l.get_a(-1, p_right) * tau * tau, 8)
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
