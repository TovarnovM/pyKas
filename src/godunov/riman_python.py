from math import *

import numpy as np


def get_e_13_1(p, ro, p_0, c_0, kappa):
    """Уравнение состояния
    
    Arguments:
        p {float} -- давление газа
        ro {float} -- плотность газа
        p_0 {float} -- параметр в уравнении состояния
        c_0 {float} -- еще один параметр там же
        kappa {float} -- коэфф. адиабаты
    
    Returns:
        float -- внутренняя энергия
    """

    return (p+kappa*p_0)/((kappa-1)*ro)-c_0*c_0/(kappa-1)

def get_p_0(ro_0, c_0, kappa):
    """Нахождение параметра в уравнении состояния p_0 по его другой "постановке" через плотность(?)
    
    Arguments:
        ro_0 {float} -- параметр в уравнеии состояния (плотность kind of)
        c_0 {float} -- параметр в уравнеии состояния
        kappa {float} -- коэфф. адиабаты
    
    Returns:
        float -- p_0 - параметр в уравнеии состояния
    """

    return ro_0*c_0*c_0/kappa

def get_R_13_2(ro, p, P, p_0, kappa):
    """Адиабата Гюгонио (формула 13.2 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        ro {float} -- плотность газа перед фронтом УВ
        p {float} -- давление газа перед фронтом УВ
        P {float} -- давление газа за фронтом УВ
        p_0 {float} -- параметр в уравнеии состояния
        kappa {float} -- коэфф. адиабаты
    
    Returns:
        float -- R - плотность газа за фронтом УВ
    """

    chsl = (kappa+1)*(P+p_0)+(kappa-1)*(p+p_0)
    znam = (kappa-1)*(P+p_0)+(kappa+1)*(p+p_0)
    return ro*chsl/znam

def get_a_13_4_shock(ro, p, P, p_0, kappa):
    """Получить массовую скорость для ударной волны 
    (формула 13.4 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        ro {float} -- плотность газа перед фронтом УВ
        p {float} -- давление газа перед фронтом УВ
        P {float} -- давление газа за фронтом УВ
        p_0 {float} -- параметр в уравнеии состояния
        kappa {float} -- коэфф. адиабаты
    
    Returns:
        float -- a - массовую скорость для ударной волны
    """

    return sqrt(ro*(0.5*(kappa+1)*(P+p_0) + 0.5*(kappa-1)*(p+p_0)))

def get_c_13_8(p, ro, p_0, kappa):
    """Получить скорость звука в газе
    (формула 13.8 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        p {float} -- давление газа
        ro {float} -- плотность газа
        p_0 {float} -- параметр в уравнеии состояния
        kappa {float} -- коэфф. адиабаты
    
    Returns:
        float -- c - скорость звука в газе
    """

    return sqrt(kappa*(p+p_0)/ro)

def get_a_13_11_discharge(ro, p, P, p_0, kappa):
    """Получить массовую скорость для волны разряжения
    (формула 13.11 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        ro {float} -- плотность газа перед фронтом волны разряжения
        p {float} -- давление газа перед фронтом волны разряжения
        P {float} -- давление газа за фронтом волны разряжения
        p_0 {float} -- параметр в уравнеии состояния
        kappa {float} -- коэфф. адиабаты
    
    Returns:
        float -- a - массовую скорость для волны разряжения
    """
    c = get_c_13_8(p, ro, p_0, kappa)
    chsl = 1 - (P+p_0)/(p+p_0)
    znam = 1 - pow((P+p_0)/(p+p_0), 0.5*(kappa-1)/kappa)
    if abs(znam)<1e-13:
        return 0
    return 0.5*(kappa-1)*ro*c/kappa*chsl/znam

def get_f_13_16(P, p_k, ro_k, p_0, kappa):
    """
    (формула 13.16 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        P {float} -- давление в месте контактного разрыва
        p_k {float} -- давление невозмущенного газа
        ro_k {float} -- плотность невозмущенного газа
        p_0 {float} -- параметр в уравнеии состояния
        kappa {float} -- коэфф. адиабаты
    
    Returns:
        float -- формула 13.16
    """

    pi_k = (P+p_0)/(p_k+p_0)
    c_k = get_c_13_8(p=p_k, ro=ro_k, p_0=p_0, kappa=kappa)
    if P >= p_k:
        znam = ro_k*c_k*sqrt((kappa+1)*pi_k*0.5/kappa+(kappa-1)*0.5/kappa)
        return (P-p_k)/znam
    else:
        return 2/(kappa-1)*c_k*(pow(pi_k, 0.5*(kappa-1)/kappa)-1)

def get_df_13_17(P, p_k, ro_k, p_0, kappa):
    """
    Производная формулы 13.16
    (формула 13.17 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        P {float} -- давление в месте контактного разрыва
        p_k {float} -- давление невозмущенного газа
        ro_k {float} -- плотность невозмущенного газа
        p_0 {float} -- параметр в уравнеии состояния
        kappa {float} -- коэфф. адиабаты
    
    Returns:
        float -- формула 13.17
    """
    pi_k = (P+p_0)/(p_k+p_0)
    c_k = get_c_13_8(p=p_k, ro=ro_k, p_0=p_0, kappa=kappa)
    if P >= p_k:
        chsl = (kappa+1)*pi_k+(3*kappa-1)
        znam = 4*kappa*ro_k*c_k*sqrt(pow((kappa+1)*pi_k*0.5/kappa+(kappa-1)*0.5/kappa, 3))
        return chsl/znam
    else:
        chsl = c_k*pow(pi_k, 0.5*(kappa-1)/kappa)
        znam = kappa*(P+p_0)
        return chsl/znam

def get_ddf_13_18(P, p_k, ro_k, p_0, kappa):
    """
    Вторая производная формулы 13.16
    (формула 13.18 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        P {float} -- давление в месте контактного разрыва
        p_k {float} -- давление невозмущенного газа
        ro_k {float} -- плотность невозмущенного газа
        p_0 {float} -- параметр в уравнеии состояния
        kappa {float} -- коэфф. адиабаты
    
    Returns:
        float -- формула 13.18
    """
    pi_k = (P+p_0)/(p_k+p_0)
    c_k = get_c_13_8(p=p_k, ro=ro_k, p_0=p_0, kappa=kappa)
    if P >= p_k:
        chsl = (kappa+1)*((kappa+1)*pi_k+(7*kappa-1))
        znam = 16*kappa*ro_k*ro_k*pow(c_k, 3)*sqrt(pow((kappa+1)*pi_k*0.5/kappa+(kappa-1)*0.5/kappa, 5))
        return -chsl/znam
    else:
        chsl = (kappa+1)*c_k*pow(pi_k, 0.5*(kappa-1)/kappa)
        znam = 2*kappa*kappa*(P+p_0)*(P+p_0)
        return -chsl/znam

def get_F_13_15(P, p_1, ro_1, p_2, ro_2, p_0, kappa):
    """Уравнение для давления в зоне контактного разрыва
    (формула 13.15 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        P {float} -- давление в месте контактного разрыва
        p_1 {float} -- давление невозмущенного газа слева от контактного разрыва (p_1 < p_2)
        ro_1 {float} -- плотность невозмущенного газа слева от контактного разрыва
        p_2 {float} -- давление невозмущенного газа справа от контактного разрыва (p_1 < p_2)
        ro_2 {float} -- плотность невозмущенного газа справа от контактного разрыва
        p_0 {float} -- параметр в уравнеии состояния
        kappa {float} -- коэфф. адиабаты
    
    Returns:
        float -- u_1-u_2  -  разница скоростей газовых потоков
    """

    return get_f_13_16(P=P, p_k=p_1, ro_k=ro_1, p_0=p_0, kappa=kappa)+\
           get_f_13_16(P=P, p_k=p_2, ro_k=ro_2, p_0=p_0, kappa=kappa)

def get_dF_13_15(P, p_1, ro_1, p_2, ro_2, p_0, kappa):
    """Производная уравнения 13.15 для давления в зоне контактного разрыва
    (формула 13.17 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        P {float} -- давление в месте контактного разрыва
        p_1 {float} -- давление невозмущенного газа слева от контактного разрыва (p_1 < p_2)
        ro_1 {float} -- плотность невозмущенного газа слева от контактного разрыва
        p_2 {float} -- давление невозмущенного газа справа от контактного разрыва (p_1 < p_2)
        ro_2 {float} -- плотность невозмущенного газа справа от контактного разрыва
        p_0 {float} -- параметр в уравнеии состояния
        kappa {float} -- коэфф. адиабаты
    
    Returns:
        float -- (u_1-u_2)'  - производная разницы скоростей газовых потоков
    """
    return get_df_13_17(P=P, p_k=p_1, ro_k=ro_1, p_0=p_0, kappa=kappa)+\
           get_df_13_17(P=P, p_k=p_2, ro_k=ro_2, p_0=p_0, kappa=kappa)

def get_ddF_13_15(P, p_1, ro_1, p_2, ro_2, p_0, kappa):
    """Вторая производная уравнения 13.15 для давления в зоне контактного разрыва
    (формула 13.18 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        P {float} -- давление в месте контактного разрыва
        p_1 {float} -- давление невозмущенного газа слева от контактного разрыва (p_1 < p_2)
        ro_1 {float} -- плотность невозмущенного газа слева от контактного разрыва
        p_2 {float} -- давление невозмущенного газа справа от контактного разрыва (p_1 < p_2)
        ro_2 {float} -- плотность невозмущенного газа справа от контактного разрыва
        p_0 {float} -- параметр в уравнеии состояния
        kappa {float} -- коэфф. адиабаты
    
    Returns:
        float -- (u_1-u_2)''  - вторая производная разницы скоростей газовых потоков
    """
    return get_ddf_13_18(P=P, p_k=p_1, ro_k=ro_1, p_0=p_0, kappa=kappa)+\
           get_ddf_13_18(P=P, p_k=p_2, ro_k=ro_2, p_0=p_0, kappa=kappa)

def get_Us_13_22(p_1, p_2, ro_1, c_1, c_2, p_0, kappa):
    """Значения функции get_F_13_15 (значения разницы скоростей газовых потоков) в точках P = -p_0, p_1, p_2, 
    для определения конфигурации, возникающих при распаде разрыва.
    (формула 13.22 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")

    Arguments:
        p_1 {float} -- давление невозмущенного газа слева от контактного разрыва (p_1 < p_2)
        p_2 {float} -- давление невозмущенного газа справа от контактного разрыва (p_1 < p_2)
        ro_1 {float} --  плотность невозмущенного газа слева от контактного разрыва
        c_1 {float} -- скорость звука в невозмущенном газе слева от контактного разрыва
        c_2 {float} -- скорость звука в невозмущенном газе справа от контактного разрыва
        p_0 {float} -- параметр в уравнеии состояния
        kappa {float} -- коэфф. адиабаты
    
    Returns:
        tuple(float, float, float) -- (U_уд, U_раз, U_вак) - критические значения разницы скоростей газовых потоков, 
                для определения конфигурации, возникающих при распаде разрыва:
                1) если u_1 - u_2 > U_уд,     -- то в обе стороны идут ударные волны;
                2) U_раз < u_1 - u_2 < U_уд   -- то вправо уданрная волна, влево волна разряжения;
                3) U_вак < u_1 - u_2 < U_раз  -- возникают 2 волны разряжения
                4) U_вак > u_1 - u_2          -- возникает область вакуума  
    """

    U_ud = (p_2-p_1)/sqrt(ro_1*((kappa+1)*(p_2+p_0)/2 + (kappa-1)*(p_1+p_0)/2))
    U_raz = -2*c_2/(kappa-1) * (1 - pow((p_1+p_0)/(p_2+p_0), 0.5*(kappa-1)/kappa))
    U_vak = -2*c_1/(kappa-1) - 2*c_2/(kappa-1)
    return U_ud, U_raz, U_vak

def mega_foo(p_1, ro_1, u_1, p_2, ro_2, u_2, p_0, kappa, eps_F=1e-2, n_iter_max=100, **kwargs):
    """Функция определения параметров газового потока в задаче Римана о распаде разрыва
    
    Arguments:
        p_1 {float} -- давление слева
        ro_1 {float} -- плотность слева
        u_1 {float} -- скорость потока слева
        p_2 {float} -- давление справа
        ro_2 {float} -- плотность справа
        u_2 {float} -- скорость потока справа
        p_0 {float} -- параметр в уравнении состояния
        kappa {float} -- параметр в уравнеии состояния
    
    Keyword Arguments:
        eps_F {float} -- точность определения решения (default: {1e-2})
        n_iter_max {int} -- максимальное количество итерайий при нахождении точного решения (default: {100})
    
    Returns:
        tuple -- success, UD_left, UD_right, D_1, D_star_1, U, D_star_2, D_2, R_1, R_2, P
                 [0] success   {bool} -- успешен ли расчет (False при создании области вакуума)
                 [1] UD_left   {bool} -- есть ли ударная волна слева
                 [2] UD_right  {bool} -- есть ли ударная волна справа
                 [3] D_1      {float} -- скорость левой границы
                 [4] D_star_1 {float} -- скорость хвоста ВР слева (если ВР слева, если УД -- D_star_1=D_1)
                 [5] U        {float} -- скорость контактного разрыва
                 [6] D_star_2 {float} -- скорость хвоста ВР справа (если ВР справа, если УД -- D_star_2=D_2)
                 [7] D_2      {float} -- скорость правой границы
                 [8] R_1      {float} -- плотность слева контактного разрыва
                 [9] R_2      {float} -- плотность справа контактного разрыва
                 [10] P        {float} -- давление в зоне контактного разрыва
    """

    reverse = False
    if (p_1 > p_2):
        reverse = True
        p_1, p_2 = p_2, p_1
        ro_1, ro_2 = ro_2, ro_1
        u_1, u_2 = -u_2, -u_1
    c_1 = get_c_13_8(p=p_1, ro=ro_1, p_0=p_0, kappa=kappa)
    c_2 = get_c_13_8(p=p_2, ro=ro_2, p_0=p_0, kappa=kappa)
    U_ud, U_raz, U_vak = get_Us_13_22(p_1, p_2, ro_1, c_1, c_2, p_0, kappa)
    if u_1-u_2 < U_vak:
        print('Вакуум !!')
        return False, 0,0,0,0,0,0,0,0,0,0
    
    #13.26

    P = (p_1*ro_2*c_2+p_2*ro_1*c_1 + (u_1-u_2)*ro_1*c_1*ro_2*c_2)/(ro_1*c_1+ro_2*c_2)
    if P < p_1:
        P = (p_1+p_0)*pow((u_1-u_2-U_vak)/(U_raz-U_vak), 2*kappa/(kappa-1)) - p_0
    for i in range(n_iter_max):
        F = get_F_13_15(P, p_1, ro_1, p_2, ro_2, p_0, kappa) - (u_1-u_2)
        if abs(F) < eps_F:
            break
        dF = get_dF_13_15(P, p_1, ro_1, p_2, ro_2, p_0, kappa)
        ddF =get_ddF_13_15(P, p_1, ro_1, p_2, ro_2, p_0, kappa)

        P = P - (F/dF)*(1+0.5*F*ddF/(dF*dF))
    else:
        print(f'Не достигнута точность={eps_F}    за {n_iter_max} итераций')

    
    # 2 уданрые волны
    if u_1-u_2 > U_ud:
        a_1 = get_a_13_4_shock(ro_1, p_1, P, p_0, kappa)
        a_2 = get_a_13_4_shock(ro_2, p_2, P, p_0, kappa)
    # УД слева, ВР справа
    elif U_raz < u_1-u_2 <= U_ud:
        a_1 = get_a_13_4_shock(ro_1, p_1, P, p_0, kappa)
        a_2 = get_a_13_11_discharge(ro_2, p_2, P, p_0, kappa)
    # 2 ВР
    else:
        a_1 = get_a_13_11_discharge(ro_1, p_1, P, p_0, kappa)
        a_2 = get_a_13_11_discharge(ro_2, p_2, P, p_0, kappa)

    # 13.14 скорость контактного разрыва
    # U = (a_1*u_1+a_2*u_2+p_1-p_2)/(a_1+a_2)
    U=0.5*(u_1 - get_f_13_16(P,p_1,ro_1,p_0,kappa) + u_2 + get_f_13_16(P,p_2,ro_2,p_0,kappa))
    # УД слева
    UD_left = u_1-u_2 > U_raz
    if UD_left:
        D_1 = u_1 - a_1/ro_1
        R_1 = get_R_13_2(ro_1, p_1, P, p_0, kappa)
        D_star_1 = D_1
    else:
        D_1 = u_1 - c_1
        c_star_1 = c_1 + 0.5*(kappa-1)*(u_1 - U)
        D_star_1 = U - c_star_1
        R_1 = kappa*(P+p_0)/pow(c_star_1, 2)

    # УД справа
    UD_right = u_1-u_2 > U_ud
    if UD_right:
        D_2 = u_2 + a_2/ro_2
        R_2 = get_R_13_2(ro_2, p_2, P, p_0, kappa)
        D_star_2 = D_2
    else:
        D_2 = u_2 + c_2
        c_star_2 = c_2 - 0.5*(kappa-1)*(u_2 - U)
        D_star_2 = U + c_star_2
        R_2 = kappa*(P+p_0)/pow(c_star_2, 2)
    if not reverse:
        return True, UD_left, UD_right, D_1, D_star_1, U, D_star_2, D_2, R_1, R_2, P
    else:
        return True, UD_right, UD_left, -D_2, -D_star_2, -U, -D_star_1, -D_1, R_2, R_1, P


def plot_rays(show=True, **init_cond):
    import matplotlib.pyplot as plt
    suc, UD_left, UD_right, D_1, D_star_1, U, D_star_2, D_2, R_1, R_2, P = mega_foo(**init_cond)
    print(UD_left, UD_right, D_1, D_star_1, U, D_star_2, D_2, R_1, R_2, P)

    line_l = init_cond.get('line_l', 1)
    lw=init_cond.get('lw', 2)
    x0 = init_cond.get('x0', 0)
    y0 = init_cond.get('y0', 0)
    y1 = y0 + line_l/sqrt(1+D_1**2)
    x1 = x0 + D_1 * y1
    if UD_left:
        plt.plot([x0,x1], [y0, y1], color='r', lw=lw, label='УД')
    else:
        plt.plot([x0,x1], [y0, y1], color='b', lw=lw, label='ВР')
        y1 = y0 + line_l/sqrt(1+D_star_1**2)
        x1 = x0 + D_star_1 * y1
        plt.plot([x0,x1], [y0, y1], color='b', lw=lw, label='ВР хвост', ls=':')
    
    y1 = y0 + line_l/sqrt(1+U**2)
    x1 = x0 + U * y1
    plt.plot([x0,x1], [y0, y1], color='green', lw=lw, label='КР', ls='--')

    y1 = y0 + line_l/sqrt(1+D_2**2)
    x1 = x0 + D_2 * y1
    if UD_right:
        plt.plot([x0,x1], [y0, y1], color='r', lw=lw, label='УД')
    else:
        plt.plot([x0,x1], [y0, y1], color='b', lw=lw, label='ВР')
        y1 = y0 + line_l/sqrt(1+D_star_2**2)
        x1 = x0 + D_star_2 * y1
        plt.plot([x0,x1], [y0, y1], color='b', lw=lw, label='ВР хвост', ls=':')
    plt.grid()
    plt.legend()
    if show:
        plt.show()

def get_distrs_to_time(t, **init_cond):
    suc, UD_left, UD_right, D_1, D_star_1, U, D_star_2, D_2, R_1, R_2, P = mega_foo(**init_cond)
    x0 = init_cond.get('x0', 0)
    n = init_cond.get('n', 0)
    c_0 = init_cond.get('c_0', 0)
    p_1, ro_1, u_1, p_2, ro_2, u_2, p_0, kappa = \
        init_cond['p_1'], init_cond['ro_1'], init_cond['u_1'], \
        init_cond['p_2'], init_cond['ro_2'], init_cond['u_2'], init_cond['p_0'], init_cond['kappa']
    width = t*D_2 - t*D_1
    x1 = t*D_1 - 0.15*width
    x2 = t*D_2 + 0.15*width
    xs = np.linspace(x1,x2,n)
    ros, ps, us, es, ms = [],[],[],[],[]
    for x in xs:
        if x <= t*D_1:
            ros.append(ro_1)
            ps.append(p_1)
            us.append(u_1)
            es.append(get_e_13_1(p_1, ro_1, p_0, c_0, kappa))
            c1 = get_c_13_8(p_1, ro_1, p_0, kappa)
            ms.append(u_1/c1)
            continue
        if (not UD_left) and (D_1*t < x <= D_star_1*t):
            k = (u_1 - U)/(D_1*t - D_star_1*t)
            b = u_1 - k*D_1*t
            u = k*x+b
            c1 = get_c_13_8(p_1, ro_1, p_0, kappa)
            # u = 2/(kappa+1)*(c1+x/t)
            ro = ro_1*pow(1-0.5*(kappa-1)*(u-u_1)/c1, 2/(kappa-1))
            p = p_1*pow(1-0.5*(kappa-1)*(u-u_1)/c1, 2*kappa/(kappa-1))
            e = get_e_13_1(p, ro, p_0, c_0, kappa)
            ros.append(ro)
            ps.append(p)
            us.append(u)
            es.append(e)
            c1=get_c_13_8(p, ro, p_0, kappa)
            ms.append(u/c1)
            continue
        if D_star_1*t < x <= U*t:
            ros.append(R_1)
            ps.append(P)
            us.append(U)
            es.append(get_e_13_1(P, R_1, p_0, c_0, kappa))
            c1=get_c_13_8(P, R_1, p_0, kappa)
            ms.append(U/c1)
            continue
        if U*t < x <= D_star_2*t:
            ros.append(R_2)
            ps.append(P)
            us.append(U)
            es.append(get_e_13_1(P, R_2, p_0, c_0, kappa))
            c1=get_c_13_8(P, R_2, p_0, kappa)
            ms.append(U/c1)
            continue
        if (not UD_right) and (D_star_2*t < x <= D_2*t):
            k = (u_2 - U)/(D_2*t - D_star_2*t)
            b = u_2 - k*D_2*t
            u = k*x+b
            c1 = get_c_13_8(p_2, ro_2, p_0, kappa)
            ro = ro_2*pow(1-0.5*(kappa-1)*(u_2-u)/c1, 2/(kappa-1))
            p = p_2*pow(1-0.5*(kappa-1)*(u_2-u)/c1, 2*kappa/(kappa-1))
            e = get_e_13_1(p, ro, p_0, c_0, kappa)
            ros.append(ro)
            ps.append(p)
            us.append(u)
            es.append(e)
            c1=get_c_13_8(p, ro, p_0, kappa)
            ms.append(u/c1)
            continue
        if x > t*D_2:
            ros.append(ro_2)
            ps.append(p_2)
            us.append(u_2)
            es.append(get_e_13_1(p_2, ro_2, p_0, c_0, kappa))
            c1=get_c_13_8(p_2, ro_2, p_0, kappa)
            ms.append(u_2/c1)
            continue
    xs = xs + x0
    return {
        'xs': xs,
        'ros': np.array(ros),
        'ps': np.array(ps),
        'us': np.array(us),
        'es': np.array(es),
        'ms': np.array(ms)
    }


def plot_distrs(**init_cond):
    import matplotlib.pyplot as plt

    def plot_one_t(t, col, lw):
        res = get_distrs_to_time(t, **init_cond)
        plt.subplot(232)
        plt.plot(res['xs'], res['ros'], c=col, lw=lw, label=f't = {t}')
        plt.grid(True)
        plt.legend()
        plt.ylabel(r'$\rho$')
        plt.subplot(233)
        plt.plot(res['xs'], res['ps'], c=col, lw=lw, label=f't = {t}')
        plt.grid(True)
        plt.legend()
        plt.ylabel(r'$p$')
        plt.subplot(234)
        plt.plot(res['xs'], res['us'], c=col, lw=lw, label=f't = {t}')
        plt.grid(True)
        plt.legend()
        plt.ylabel(r'$u$')
        plt.xlabel(r'$x$')
        plt.subplot(235)
        plt.plot(res['xs'], res['ms'], c=col, lw=lw, label=f't = {t}')
        plt.grid(True)
        plt.legend()
        plt.ylabel(r'$M$')
        plt.xlabel(r'$x$')
        plt.subplot(236)
        plt.plot(res['xs'], res['es'], c=col, lw=lw, label=f't = {t}')
        plt.grid(True)
        plt.legend()
        plt.ylabel(r'$e$')
        plt.xlabel(r'$x$')
    for t, col in zip(init_cond['ts'], plt.get_cmap('Set1').colors):
        plot_one_t(t, col, init_cond.get('lw', 2))
    plt.subplot(231)
    plot_rays(show=False, **init_cond)
    plt.show()

def get_init_conds_4_tsts():
    # из статьи ОДНОМЕРНЫЕ ЗАДАЧИ ГАЗОВОЙ ДИНАМИКИ И ИХ РЕШЕНИЕ  ПРИ ПОМОЩИ РАЗНОСТНЫХ СХЕМ ВЫСОКОЙ  РАЗРЕШАЮЩЕЙ СПОСОБНОСТИ 
    #           LEFT          RIGHT              
    #      ro       u       p       ro      u       p       t
    v = [
         ( 1,       0,      1,      0.125,  0,      0.1,    0.15  ),
         ( 0.445,   0.698,  3.528,  0.5,    0,      0.571,  0.15  ),       
         ( 1,       0,      1,      0.02,   3.55,   1,      0.15  ),
         ( 3.857,   0.920,  10.333, 1,      3.55,   1,      0.09  ),
         ( 1,       0,      0.5,    0.5,    0,      0.5,    0.42  ),
         ( 1,       0.5,    0.5,    0.5,    0.5,    0.5,    0.43  ),
         ( 1,       -1,     1,      0.9275, -1.0781,0.9,    0.18  ),
         ( 1,       0,      1000,   1,      0,      0.01,   0.012  ),
         ( 10,      2000,   500,    20,     0,      500,    0.012  ),
         ( 1,       -2,     0.4,    1,      2,      0.4,    0.15  )
        ]
    for ro_1, u_1, p_1, ro_2, u_2, p_2, t in v:
        yield {
            'p_1' : p_1,          # давление слева
            'ro_1' : ro_1,         # плотность слева   
            'u_1' : u_1,          # скорость слева      
            'p_2' : p_2,        # давление справа     
            'ro_2': ro_2,      # плотность справа    
            'u_2' : u_2,          # скорость справа        
            'p_0' : 0,          # параметр в ур-ии состояния        
            'kappa' : 1.4,      # параметр в ур-ии состояния          
            'c_0' : 0,          # параметр в ур-ии состояния     
            'eps_F':1e-6,       # точность определения решения           
            'n_iter_max':100,   # максимальное количество итераций при определении решения              
            'x0' : 0.5,         # положение разрыва в момент t=0        
            'ts': [0.15],        # времена, в которых нужно построить графики распределения параметров         
            'n': 10000          # кол-во точек, в которых ищутся параметры волны разрежения         
        } 

  

if __name__ == "__main__":
    # init_cond = {
    #     'p_1' : 1,          # давление слева
    #     'ro_1' : 1,         # плотность слева   
    #     'u_1' : 0,          # скорость слева      
    #     'p_2' : 0.1,        # давление справа     
    #     'ro_2': 0.125,      # плотность справа    
    #     'u_2' : 0,          # скорость справа        
    #     'p_0' : 0,          # параметр в ур-ии состояния        
    #     'kappa' : 1.4,      # параметр в ур-ии состояния          
    #     'c_0' : 0,          # параметр в ур-ии состояния     
    #     'eps_F':1e-3,       # точность определения решения           
    #     'n_iter_max':100,   # максимальное количество итераций при определении решения              
    #     'x0' : 0.5,         # положение разрыва в момент t=0        
    #     'ts': [0.2],        # времена, в которых нужно построить графики распределения параметров         
    #     'n': 10000          # кол-во точек, в которых ищутся параметры волны разрежения         
    # }

    # plot_rays(**init_cond)
    for i, ic in enumerate(get_init_conds_4_tsts()):
        try:
            plot_distrs(**ic)
        except Exception as e:
            print(i+1, e)
