import numpy as np


def AUSMp(l):
    """
    Схема AUSM плюс для получения потоков f на границах ячеек
    :param l: слой
    :return:  список нумпи массивов f - потоков через границы (не считая левой и правой границы сетки)
    """
    r1, u1, e1 = l.get_param(l.q)
    p1 = l.p
    H1 = e1 + 0.5 * np.square(u1) + p1 / r1
    c1 = l.get_Csound(r1, p1)

    r2 = np.roll(r1, -1)
    u2 = np.roll(u1, -1)
    p2 = np.roll(p1, -1)
    H2 = np.roll(H1, -1)
    c2 = np.roll(c1, -1)

    cs = 0.5 * (c1 + c2)

    V_right = np.roll(l.V, -1)[:-1]

    Mr1 = (u1 - V_right) / cs
    Mr2 = (u2 - V_right) / cs

    cond1 = (np.abs(Mr1)) >= 1
    cond2 = (np.abs(Mr2)) >= 1

    M4p = np.zeros_like(Mr1, dtype=np.float64)
    M4p[cond1] = 0.5 * (Mr1[cond1] + np.abs(Mr1[cond1]))
    M4p[cond1 == False] = 0.25 * np.square(Mr1[cond1 == False] + 1) * (
        1 + 2 * 0.25 * np.square(Mr1[cond1 == False] - 1))

    P5p = np.zeros_like(Mr1, dtype=np.float64)
    P5p[cond1] = 0.5 * (Mr1[cond1] + np.abs(Mr1[cond1])) / Mr1[cond1]
    P5p[cond1 == False] = 0.25 * np.square(Mr1[cond1 == False] + 1) * \
                          ((2 - Mr1[cond1 == False]) + 3 * 0.25 * Mr1[cond1 == False] * np.square(
                              Mr1[cond1 == False] - 1))

    M4m = np.zeros_like(Mr1, dtype=np.float64)
    M4m[cond2] = 0.5 * (Mr2[cond2] - np.abs(Mr2[cond2]))
    M4m[cond2 == False] = -0.25 * np.square(Mr2[cond2 == False] - 1) * (
        1 + 2 * 0.25 * np.square(Mr2[cond2 == False] + 1))

    P5m = np.zeros_like(Mr1, dtype=np.float64)
    P5m[cond2] = 0.5 * (Mr2[cond2] - np.abs(Mr2[cond2])) / Mr2[cond2]
    P5m[cond2 == False] = 0.25 * np.square(Mr2[cond2 == False] - 1) * \
                          ((2 + Mr2[cond2 == False]) + 3 * 0.25 * Mr2[cond2 == False] * np.square(
                              Mr2[cond2 == False] + 1))

    Mrf = M4p + M4m
    pf = P5p * p1 + P5m * p2

    f1 = 0.5 * (cs * Mrf * (r1 + r2) - cs * np.abs(Mrf) * (r2 - r1))
    f2 = 0.5 * (cs * Mrf * (r1 * u1 + r2 * u2) - cs * np.abs(Mrf) * (r2 * u2 - r1 * u1)) + pf
    f3 = 0.5 * (cs * Mrf * (r1 * H1 + r2 * H2) - cs * np.abs(Mrf) * (r2 * H2 - r1 * H1)) + pf * V_right

    return [f1, f2, f3]


def get_flux_bord(l, i, func_bord):
    """
    Схема AUSM плюс для получения потоков f на границах сетки
    :param l: слой
    :param i: индекс 0 или -1 для левой и правой границы соответственно
    :param func_bord: функция граничных условий
    :return: кортеж потоков f через границы
    """
    if i == -1:
        q = [l.q[0][i], l.q[1][i], l.q[2][i]]
        r1, u1, e1 = l.get_param(q)
        p1 = l.p[i]
        H1 = e1 + 0.5 * (u1 ** 2) + p1 / r1
        c1 = l.get_Csound(r1, p1)

        V = l.V[i]

        r2, u2, p2 = func_bord(r1, u1, p1, V)

        e2 = l.get_energ(p2, r2)
        H2 = e2 + 0.5 * (u2 ** 2) + p2 / r2
        c2 = l.get_Csound(r2, p2)

    else:
        q = [l.q[0][i], l.q[1][i], l.q[2][i]]
        r2, u2, e2 = l.get_param(q)
        p2 = l.p[i]
        H2 = e2 + 0.5 * (u2 ** 2) + p2 / r2
        c2 = l.get_Csound(r2, p2)

        V = l.V[i]

        r1, u1, p1 = func_bord(r2, u2, p2, V)

        e1 = l.get_energ(p1, r1)
        H1 = e1 + 0.5 * (u1 ** 2) + p1 / r1
        c1 = l.get_Csound(r1, p1)

    cs = 0.5 * (c1 + c2)

    Mr1 = (u1 - V) / cs
    Mr2 = (u2 - V) / cs

    if abs(Mr1) >= 1:
        M4p = 0.5 * (Mr1 + abs(Mr1))
        P5p = 0.5 * (Mr1 + abs(Mr1)) / Mr1
    else:
        M4p = 0.25 * ((Mr1 + 1) ** 2) * (1 + 2 * 0.25 * ((Mr1 - 1) ** 2))
        P5p = 0.25 * ((Mr1 + 1) ** 2) * ((2 - Mr1) + 3 * 0.25 * Mr1 * ((Mr1 - 1) ** 2))

    if abs(Mr2) >= 1:
        M4m = 0.5 * (Mr2 - abs(Mr2))
        P5m = 0.5 * (Mr2 - abs(Mr2)) / Mr2
    else:
        M4m = -0.25 * ((Mr2 - 1) ** 2) * (1 + 2 * 0.25 * ((Mr2 + 1) ** 2))
        P5m = 0.25 * ((Mr2 - 1) ** 2) * ((2 + Mr2) - 3 * 0.25 * Mr2 * ((Mr2 + 1) ** 2))

    Mrf = M4p + M4m
    pf = P5p * p1 + P5m * p2

    f1 = 0.5 * (cs * Mrf * (r1 + r2) - cs * abs(Mrf) * (r2 - r1))
    f2 = 0.5 * (cs * Mrf * (r1 * u1 + r2 * u2) - cs * abs(Mrf) * (r2 * u2 - r1 * u1)) + pf
    f3 = 0.5 * (cs * Mrf * (r1 * H1 + r2 * H2) - cs * abs(Mrf) * (r2 * H2 - r1 * H1)) + pf * V

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
