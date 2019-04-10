# distutils: language=c++
# cython: language_level=3, boundscheck=False, nonecheck=False, cdivision=True

import cython
from libc.math cimport pi, sqrt, abs, copysign, pow
from libcpp cimport bool as bool_t

cpdef inline double get_e_13_1(double p, double ro, double p_0, double c_0, double gamma):
    """Уравнение состояния
    (формула 13.1 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")

    Arguments:
        p {float} -- давление газа
        ro {float} -- плотность газа
        p_0 {float} -- параметр в уравнении состояния
        c_0 {float} -- еще один параметр там же
        gamma {float} -- коэфф. адиабаты
    
    Returns:
        float -- внутренняя энергия
    """
    # return p*(1/ro-0.001)/(gamma-1)
    return (p+gamma*p_0)/((gamma-1)*ro)-c_0*c_0/(gamma-1)

cpdef inline double get_p_0(double ro_0, double c_0, double gamma):
    """Нахождение параметра в уравнении состояния p_0 по его другой "постановке" через плотность(?)
    
    Arguments:
        ro_0 {float} -- параметр в уравнеии состояния (плотность kind of)
        c_0 {float} -- параметр в уравнеии состояния
        gamma {float} -- коэфф. адиабаты
    
    Returns:
        float -- p_0 - параметр в уравнеии состояния
    """

    return ro_0*c_0*c_0/gamma

cpdef inline double get_R_13_2(double ro, double p, double P, double p_0, double gamma):
    """Адиабата Гюгонио (формула 13.2 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        ro {float} -- плотность газа перед фронтом УВ
        p {float} -- давление газа перед фронтом УВ
        P {float} -- давление газа за фронтом УВ
        p_0 {float} -- параметр в уравнеии состояния
        gamma {float} -- коэфф. адиабаты
    
    Returns:
        float -- R - плотность газа за фронтом УВ
    """

    cdef double chsl = (gamma+1)*(P+p_0)+(gamma-1)*(p+p_0)
    cdef double znam = (gamma-1)*(P+p_0)+(gamma+1)*(p+p_0)
    return ro*chsl/znam

cpdef inline double get_a_13_4_shock(double ro, double p, double P, double p_0, double gamma):
    """Получить массовую скорость для ударной волны 
    (формула 13.4 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        ro {float} -- плотность газа перед фронтом УВ
        p {float} -- давление газа перед фронтом УВ
        P {float} -- давление газа за фронтом УВ
        p_0 {float} -- параметр в уравнеии состояния
        gamma {float} -- коэфф. адиабаты
    
    Returns:
        float -- a - массовую скорость для ударной волны
    """

    return sqrt(ro*(0.5*(gamma+1)*(P+p_0) + 0.5*(gamma-1)*(p+p_0)))

cpdef inline double get_c_13_8(double p, double ro, double p_0, double gamma):
    """Получить скорость звука в газе
    (формула 13.8 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        p {float} -- давление газа
        ro {float} -- плотность газа
        p_0 {float} -- параметр в уравнеии состояния
        gamma {float} -- коэфф. адиабаты
    
    Returns:
        float -- c - скорость звука в газе
    """
    # return sqrt(p / ((1/gamma) * ro * (1 - 0.001*ro)))
    return sqrt(gamma*(p+p_0)/ro)

cpdef double get_a_13_11_discharge(double ro, double p, double P, double c, double p_0, double gamma):
    """Получить массовую скорость для волны разряжения
    (формула 13.11 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        ro {float} -- плотность газа перед фронтом волны разряжения
        p {float} -- давление газа перед фронтом волны разряжения
        P {float} -- давление газа за фронтом волны разряжения
        c {float} -- скорость звука перед фронтом волны разряжения
        p_0 {float} -- параметр в уравнеии состояния
        gamma {float} -- коэфф. адиабаты
    
    Returns:
        float -- a - массовую скорость для волны разряжения
    """
    # c = get_c_13_8(p, ro, p_0, gamma)
    cdef double chsl = 1.0 - (P+p_0)/(p+p_0)
    cdef double znam = 1.0 - pow((P+p_0)/(p+p_0), 0.5*(gamma-1)/gamma)
    if abs(znam)<1e-13:
        return 0
    return 0.5*(gamma-1)*ro*c/gamma*chsl/znam

cpdef double get_f_13_16(double P, double p_k, double ro_k, double c_k, double p_0, double gamma):
    """
    (формула 13.16 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        P {float} -- давление в месте контактного разрыва
        p_k {float} -- давление невозмущенного газа
        ro_k {float} -- плотность невозмущенного газа
        c_k {float} -- скорость звука в невозмущенном газе
        p_0 {float} -- параметр в уравнеии состояния
        gamma {float} -- коэфф. адиабаты
    
    Returns:
        float -- формула 13.16
    """

    cdef double pi_k = (P+p_0)/(p_k+p_0)
    # c_k = get_c_13_8(p=p_k, ro=ro_k, p_0=p_0, gamma=gamma)
    cdef double znam
    if P >= p_k:
        znam = ro_k*c_k*sqrt((gamma+1)*pi_k*0.5/gamma+(gamma-1)*0.5/gamma)
        return (P-p_k)/znam
    else:
        return 2/(gamma-1)*c_k*(pow(pi_k, 0.5*(gamma-1)/gamma)-1)

cpdef double get_df_13_17(double P, double p_k, double ro_k, double c_k, double p_0, double gamma):
    """
    Производная формулы 13.16
    (формула 13.17 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        P {float} -- давление в месте контактного разрыва
        p_k {float} -- давление невозмущенного газа
        ro_k {float} -- плотность невозмущенного газа
        c_k {float} -- скорость звука в невозмущенном газе
        p_0 {float} -- параметр в уравнеии состояния
        gamma {float} -- коэфф. адиабаты
    
    Returns:
        float -- формула 13.17
    """
    cdef double pi_k = (P+p_0)/(p_k+p_0)
    # c_k = get_c_13_8(p=p_k, ro=ro_k, p_0=p_0, gamma=gamma)
    cdef double chsl, znam
    if P >= p_k:
        chsl = (gamma+1)*pi_k+(3*gamma-1)
        znam = 4*gamma*ro_k*c_k*sqrt(pow((gamma+1)*pi_k*0.5/gamma+(gamma-1)*0.5/gamma, 3))
        return chsl/znam
    else:
        chsl = c_k*pow(pi_k, 0.5*(gamma-1)/gamma)
        znam = gamma*(P+p_0)
        return chsl/znam

cpdef double get_ddf_13_18(double P, double p_k, double ro_k, double c_k, double p_0, double gamma):
    """
    Вторая производная формулы 13.16
    (формула 13.18 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        P {float} -- давление в месте контактного разрыва
        p_k {float} -- давление невозмущенного газа
        ro_k {float} -- плотность невозмущенного газа
        c_k {float} -- скорость звука в невозмущенном газе
        p_0 {float} -- параметр в уравнеии состояния
        gamma {float} -- коэфф. адиабаты
    
    Returns:
        float -- формула 13.18
    """
    cdef double pi_k = (P+p_0)/(p_k+p_0)
    cdef double chsl, znam
    # c_k = get_c_13_8(p=p_k, ro=ro_k, p_0=p_0, gamma=gamma)
    if P >= p_k:
        chsl = (gamma+1)*((gamma+1)*pi_k+(7*gamma-1))
        znam = 16*gamma*ro_k*ro_k*pow(c_k, 3)*sqrt(pow((gamma+1)*pi_k*0.5/gamma+(gamma-1)*0.5/gamma, 5))
        return -chsl/znam
    else:
        chsl = (gamma+1)*c_k*pow(pi_k, 0.5*(gamma-1)/gamma)
        znam = 2*gamma*gamma*(P+p_0)*(P+p_0)
        return -chsl/znam

cpdef inline double get_F_13_15(double P, double p_1, double ro_1, double c_1, double p_2, double ro_2, double c_2, double p_0, double gamma):
    """Уравнение для давления в зоне контактного разрыва
    (формула 13.15 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        P {float} -- давление в месте контактного разрыва
        p_1 {float} -- давление невозмущенного газа слева от контактного разрыва (p_1 < p_2)
        ro_1 {float} -- плотность невозмущенного газа слева от контактного разрыва
        c_1 {float} -- скорость звука в газе слева от контактного разрыва
        p_2 {float} -- давление невозмущенного газа справа от контактного разрыва (p_1 < p_2)
        ro_2 {float} -- плотность невозмущенного газа справа от контактного разрыва
        c_2 {float} -- скорость звука в газе справа от контактного разрыва
        p_0 {float} -- параметр в уравнеии состояния
        gamma {float} -- коэфф. адиабаты
    
    Returns:
        float -- u_1-u_2  -  разница скоростей газовых потоков
    """

    return get_f_13_16(P=P, p_k=p_1, ro_k=ro_1, c_k=c_1, p_0=p_0, gamma=gamma)+\
           get_f_13_16(P=P, p_k=p_2, ro_k=ro_2, c_k=c_2, p_0=p_0, gamma=gamma)

cpdef inline double get_dF_13_15(double P, double p_1, double ro_1, double c_1, double p_2, double ro_2, double c_2, double p_0, double gamma):
    """Производная уравнения 13.15 для давления в зоне контактного разрыва
    (формула 13.17 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        P {float} -- давление в месте контактного разрыва
        p_1 {float} -- давление невозмущенного газа слева от контактного разрыва (p_1 < p_2)
        ro_1 {float} -- плотность невозмущенного газа слева от контактного разрыва
        c_1 {float} -- скорость звука в газе слева от контактного разрыва
        p_2 {float} -- давление невозмущенного газа справа от контактного разрыва (p_1 < p_2)
        ro_2 {float} -- плотность невозмущенного газа справа от контактного разрыва
        c_2 {float} -- скорость звука в газе справа от контактного разрыва
        p_0 {float} -- параметр в уравнеии состояния
        gamma {float} -- коэфф. адиабаты
    
    Returns:
        float -- (u_1-u_2)'  - производная разницы скоростей газовых потоков
    """
    return get_df_13_17(P=P, p_k=p_1, ro_k=ro_1, c_k=c_1, p_0=p_0, gamma=gamma)+\
           get_df_13_17(P=P, p_k=p_2, ro_k=ro_2, c_k=c_2, p_0=p_0, gamma=gamma)

cpdef inline double get_ddF_13_15(double P, double p_1, double ro_1, double c_1, double p_2, double ro_2, double c_2, double p_0, double gamma):
    """Вторая производная уравнения 13.15 для давления в зоне контактного разрыва
    (формула 13.18 из монографии С.К. Годунова "Численное решение многомерных задач газовой динамики")
    
    Arguments:
        P {float} -- давление в месте контактного разрыва
        p_1 {float} -- давление невозмущенного газа слева от контактного разрыва (p_1 < p_2)
        ro_1 {float} -- плотность невозмущенного газа слева от контактного разрыва
        c_1 {float} -- скорость звука в газе слева от контактного разрыва
        p_2 {float} -- давление невозмущенного газа справа от контактного разрыва (p_1 < p_2)
        ro_2 {float} -- плотность невозмущенного газа справа от контактного разрыва
        c_2 {float} -- скорость звука в газе справа от контактного разрыва
        p_0 {float} -- параметр в уравнеии состояния
        gamma {float} -- коэфф. адиабаты
    
    Returns:
        float -- (u_1-u_2)''  - вторая производная разницы скоростей газовых потоков
    """
    return get_ddf_13_18(P=P, p_k=p_1, ro_k=ro_1, c_k=c_1, p_0=p_0, gamma=gamma)+\
           get_ddf_13_18(P=P, p_k=p_2, ro_k=ro_2, c_k=c_2, p_0=p_0, gamma=gamma)

cpdef (double, double, double) get_Us_13_22(double p_1, double p_2, double ro_1, double c_1, double c_2, double p_0, double gamma):
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
        gamma {float} -- коэфф. адиабаты
    
    Returns:
        tuple(float, float, float) -- (U_уд, U_раз, U_вак) - критические значения разницы скоростей газовых потоков, 
                для определения конфигурации, возникающих при распаде разрыва:
                1) если u_1 - u_2 > U_уд,     -- то в обе стороны идут ударные волны;
                2) U_раз < u_1 - u_2 < U_уд   -- то вправо уданрная волна, влево волна разряжения;
                3) U_вак < u_1 - u_2 < U_раз  -- возникают 2 волны разряжения
                4) U_вак > u_1 - u_2          -- возникает область вакуума  
    """

    cdef double U_ud = (p_2-p_1)/sqrt(ro_1*((gamma+1)*(p_2+p_0)/2 + (gamma-1)*(p_1+p_0)/2))
    cdef double U_raz = -2*c_2/(gamma-1) * (1 - pow((p_1+p_0)/(p_2+p_0), 0.5*(gamma-1)/gamma))
    cdef double U_vak = -2*c_1/(gamma-1) - 2*c_2/(gamma-1)
    return U_ud, U_raz, U_vak

cpdef (double, double, double, double, double, double, double, double, double, double, double) \
    mega_foo(double p_1, double ro_1, double u_1, double c_1, \
             double p_2, double ro_2, double u_2, double c_2, \
             double p_0, double gamma, double eps_F=1e-5, int n_iter_max=100):
    """Функция определения параметров газового потока в задаче Римана о распаде разрыва
    
    Arguments:
        p_1 {float} -- давление слева
        ro_1 {float} -- плотность слева
        u_1 {float} -- скорость потока слева
        c_1 {float} -- скорость звука в невозмущенном газе слева от контактного разрыва
        p_2 {float} -- давление справа
        ro_2 {float} -- плотность справа
        u_2 {float} -- скорость потока справа
        c_2 {float} -- скорость звука в невозмущенном газе справа от контактного разрыва
        p_0 {float} -- параметр в уравнении состояния
        gamma {float} -- параметр в уравнеии состояния (коэфф. адиабаты)
    
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

    cdef bool_t reverse = False
    if (p_1 > p_2):
        reverse = True
        p_1, p_2 = p_2, p_1
        ro_1, ro_2 = ro_2, ro_1
        u_1, u_2 = -u_2, -u_1
        c_1, c_2 = c_2, c_1
    # c_1 = get_c_13_8(p=p_1, ro=ro_1, p_0=p_0, gamma=gamma)
    # c_2 = get_c_13_8(p=p_2, ro=ro_2, p_0=p_0, gamma=gamma)
    cdef double U_ud, U_raz, U_vak
    U_ud, U_raz, U_vak = get_Us_13_22(p_1, p_2, ro_1, c_1, c_2, p_0, gamma)
    if u_1-u_2 < U_vak:
        # print('Вакуум !!')
        return False, 0,0,0,0,0,0,0,0,0,0
    
    #13.26
    cdef double F, dF, ddF
    cdef double P = (p_1*ro_2*c_2+p_2*ro_1*c_1 + (u_1-u_2)*ro_1*c_1*ro_2*c_2)/(ro_1*c_1+ro_2*c_2)
    if P < p_1:
        P = (p_1+p_0)*pow((u_1-u_2-U_vak)/(U_raz-U_vak), 2*gamma/(gamma-1)) - p_0
    for i in range(n_iter_max):
        F = get_F_13_15(P, p_1, ro_1, c_1, p_2, ro_2, c_2, p_0, gamma) - (u_1-u_2)
        if abs(F) < eps_F:
            break
        dF = get_dF_13_15(P, p_1, ro_1, c_1, p_2, ro_2, c_2, p_0, gamma)
        ddF =get_ddF_13_15(P, p_1, ro_1, c_1, p_2, ro_2, c_2, p_0, gamma)

        P = P - (F/dF)*(1+0.5*F*ddF/(dF*dF))
    # else:
    #     print(f'Не достигнута точность={eps_F}    за {n_iter_max} итераций')

    cdef double a_1, a_2
    # 2 уданрые волны
    if u_1-u_2 > U_ud:
        a_1 = get_a_13_4_shock(ro_1, p_1, P, p_0, gamma)
        a_2 = get_a_13_4_shock(ro_2, p_2, P, p_0, gamma)
    # УД слева, ВР справа
    elif U_raz < u_1-u_2 <= U_ud:
        a_1 = get_a_13_4_shock(ro_1, p_1, P, p_0, gamma)
        a_2 = get_a_13_11_discharge(ro_2, p_2, P, c_2, p_0, gamma)
    # 2 ВР
    else:
        a_1 = get_a_13_11_discharge(ro_1, p_1, P, c_1, p_0, gamma)
        a_2 = get_a_13_11_discharge(ro_2, p_2, P, c_2, p_0, gamma)

    # 13.14 скорость контактного разрыва
    # U = (a_1*u_1+a_2*u_2+p_1-p_2)/(a_1+a_2)
    cdef double U=0.5*(u_1 - get_f_13_16(P,p_1,ro_1,c_1,p_0,gamma) + u_2 + get_f_13_16(P,p_2,ro_2,c_2,p_0,gamma))
    
    cdef double D_1, R_1, D_star_1, D_2, D_star_2, R_2
    # УД слева
    cdef double UD_left = u_1-u_2 > U_raz
    if UD_left:
        D_1 = u_1 - a_1/ro_1
        R_1 = get_R_13_2(ro_1, p_1, P, p_0, gamma)
        D_star_1 = D_1
    else:
        D_1 = u_1 - c_1
        c_star_1 = c_1 + 0.5*(gamma-1)*(u_1 - U)
        D_star_1 = U - c_star_1
        R_1 = gamma*(P+p_0)/pow(c_star_1, 2)

    # УД справа
    cdef double UD_right = u_1-u_2 > U_ud
    if UD_right:
        D_2 = u_2 + a_2/ro_2
        R_2 = get_R_13_2(ro_2, p_2, P, p_0, gamma)
        D_star_2 = D_2
    else:
        D_2 = u_2 + c_2
        c_star_2 = c_2 - 0.5*(gamma-1)*(u_2 - U)
        D_star_2 = U + c_star_2
        R_2 = gamma*(P+p_0)/pow(c_star_2, 2)
    if not reverse:
        return True, UD_left, UD_right, D_1, D_star_1, U, D_star_2, D_2, R_1, R_2, P
    else:
        return True, UD_right, UD_left, -D_2, -D_star_2, -U, -D_star_1, -D_1, R_2, R_1, P

cpdef (double, double, double) get_ray_URP(
    double ray_W, double UD_left, double UD_right, double D_1, double D_star_1, double U, 
    double D_star_2, double D_2, double R_1, double R_2, double P, 
    double p_1, double ro_1, double u_1, double c_1, 
    double p_2, double ro_2, double u_2, double c_2, double gamma):
    """ОООооо даааа =) функция получения вектора характеристик газа, испытавшего распад разрыва, на подвижной границе.
    В начальный момент времени граница находится в точке разрыва. Скорость границы (rayW) постоянна
    
    Arguments:
        ray_W    {[type]} -- скорость подвижной границы
        UD_left    {bool} -- есть ли ударная волна слева
        UD_right   {bool} -- есть ли ударная волна справа
        D_1       {float} -- скорость левой границы
        D_star_1  {float} -- скорость хвоста ВР слева (если ВР слева, если УД -- D_star_1=D_1)
        U         {float} -- скорость контактного разрыва
        D_star_2  {float} -- скорость хвоста ВР справа (если ВР справа, если УД -- D_star_2=D_2)
        D_2       {float} -- скорость правой границы
        R_1       {float} -- плотность слева контактного разрыва
        R_2       {float} -- плотность справа контактного разрыва
        P         {float} -- давление в зоне контактного разрыва
        p_1       {float} -- давление невозмущенного потока слева (левее D_1)
        ro_1      {float} -- плотность невозмущенного потока слева (левее D_1)
        u_1       {float} -- скорость невозмущенного потока слева (левее D_1)
        c_1       {float} -- скорость звука в невозмущенном газе слева от контактного разрыва
        p_2       {float} -- давление невозмущенного потока справа (правее D_2)
        ro_2      {float} -- плотность невозмущенного потока справа (правее D_2)
        u_2       {float} -- скорость невозмущенного потока справа (правее D_2)
        c_2       {float} -- скорость звука в невозмущенном газе справа от контактного разрыва
        gamma     {float} -- параметр в уравнеии состояния (коэфф. адиабаты)
    
    Returns:
        tuple(float, float, float) -- вектор характеристик потока. 
               U,     R,     P      (скорость, плотность, давление, внутр. энергия)
    """
    cdef double resU,resR,resP,k,b
    resU,resR,resP = 0,0,0
    if D_star_1 < ray_W <= U:
        resR = R_1
        resP = P
        resU = U
    elif U < ray_W <= D_star_2:
        resR = R_2
        resP = P
        resU = U
    elif (not UD_left) and (D_1 < ray_W <= D_star_1):
        k = (u_1 - U)/(D_1 - D_star_1)
        b = u_1 - k*D_1
        resU = k*ray_W+b
        resR = ro_1*pow(1-0.5*(gamma-1)*(resU-u_1)/c_1, 2/(gamma-1))
        resP = p_1*pow(1-0.5*(gamma-1)*(resU-u_1)/c_1, 2*gamma/(gamma-1))
    elif (not UD_right) and (D_star_2 < ray_W <= D_2):
        k = (u_2 - U)/(D_2 - D_star_2)
        b = u_2 - k*D_2
        resU = k*ray_W+b
        resR = ro_2*pow(1-0.5*(gamma-1)*(u_2-resU)/c_2, 2/(gamma-1))
        resP = p_2*pow(1-0.5*(gamma-1)*(u_2-resU)/c_2, 2*gamma/(gamma-1))
    elif ray_W > D_2:
        resR=ro_2
        resP=p_2
        resU=u_2
    elif ray_W <= D_1:
        resR=ro_1
        resP=p_1
        resU=u_1
    return resU, resR, resP