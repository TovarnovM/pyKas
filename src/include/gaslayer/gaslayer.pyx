# distutils: language=c++
# cython: language_level=3, boundscheck=False, nonecheck=False, cdivision=True, initializedcheck=False

from tube cimport Tube, InterpXY
import cython
from libc.math cimport pi, sqrt, copysign, exp, fabs
from godunov cimport get_e_13_1, get_p_13_1, get_c_13_8, mega_foo_cython, \
                     get_ray_URP, MegaFooResult, border_wall_URPDD_result, Border_URPDD_Result, \
                     get_p_0, border_wall_fill_rr, mega_foo_fill_rr
import numpy as np
cimport numpy as np


cpdef double foo():
    """функция для тестов
    
    Returns:
        [type] -- [description]
    """
    cdef Tube t = Tube([1,2,3],[3,3,3])
    cdef double d = t.get_d(4)
    return d

cpdef inline double abs(double x) nogil:
    if x < 0:
        return -x
    return x

cpdef inline double min2(double a, double b) nogil:
    if a < b:
        return a
    else:
        return b

cpdef inline double max2(double a, double b) nogil:
    if a > b:
        return a
    else:
        return b

cpdef inline void roue_to_q_(
    double ro,
    double u,
    double e,
    double[:] q) nogil:

    q[0] = ro
    q[1] = ro * u
    q[2] = ro * (e + 0.5 * u * u)

cpdef inline (double, double, double) roue_to_q(
    double ro,
    double u,
    double e) nogil:

    return ro, ro * u, ro * (e + 0.5 * u * u)


cpdef inline (double, double, double) q_to_roue(double[:] q) nogil:
    cdef double ro = q[0]
    cdef double u = q[1] / q[0]
    cdef double e = q[2] / q[0] - 0.5 * u * u
    return ro, u, e

cpdef inline (double, double, double) q123_to_roue(double q1, double q2, double q3) nogil:
    cdef double ro = q1
    cdef double u = q2 / q1
    cdef double e = q3 / q1 - 0.5 * u * u
    return ro, u, e


cpdef inline double roe_to_p(
    double ro,
    double e,
    double gamma,
    double b) nogil:

    return (gamma-1)*e*ro/(1-b*ro)



cpdef inline double rop_to_e(
    double ro,
    double p,
    double gamma,
    double b) nogil:

    return p*(1/ro-b)/(gamma-1)


cpdef inline double rop_to_csound(
    double ro,
    double p,
    double gamma,
    double b) nogil:

    return sqrt(p / ((1/gamma) * ro * (1 - b*ro)))



cpdef (double, double, double) AUSM_gas_(
    double p1, 
    double ro1, 
    double u1, 
    double e1, 
    double c1,
    double p2, 
    double ro2, 
    double u2, 
    double e2, 
    double c2,
    double vbi) nogil:

    cdef double r1=ro1
    cdef double r2=ro2
    cdef double H1 = e1 + 0.5*u1*u1 + p1/r1
    cdef double H2 = e2 + 0.5*u2*u2 + p2/r2

    cdef double cs = 0.5*(c1+c2)
    cdef double Mr1 = (u1-vbi)/cs
    cdef double Mr2 = (u2-vbi)/cs

    # ! Vacuum solution (?)
    # uvac = 2.0*g(4)*cs - du
    # if (uvac <= 0) then
    #     write(*,*) "Vacuum generated by given data"
    #     return
    # end if

    cdef double M4p, P5p, M4m, P5m
    if abs(Mr1) >= 1.0 :
        M4p = 0.5*(Mr1+abs(Mr1))
        P5p = 0.5*(Mr1+abs(Mr1))/Mr1    
    else:
        M4p = 0.25*((Mr1+1.0)*(Mr1+1.0))*(1.0+2.0*0.25*(Mr1-1.0)*(Mr1-1.0))
        P5p = 0.25*((Mr1+1.0)*(Mr1+1.0))*((2.0-Mr1)+3.0*Mr1*0.25*(Mr1-1.0)*(Mr1-1.0))   

    if abs(Mr2) >= 1.0:
        M4m = 0.5*(Mr2-abs(Mr2))
        P5m = 0.5*(Mr2-abs(Mr2))/Mr2
    else:
        M4m = -0.25*((Mr2-1.0)*(Mr2-1.0))*(1.0+2.0*0.25*(Mr2+1.0)*(Mr2+1.0))
        P5m = 0.25*((Mr2-1.0)*(Mr2-1.0))*((2.0+Mr2)-3.0*Mr2*0.25*(Mr2+1.0)*(Mr2+1.0))
    
    cdef double phi1 = 1.0
    cdef double phi2 = 1.0
    cdef double Mrf = M4p + M4m
    cdef double pf = P5p*phi1*p1 + P5m*phi2*p2

    cdef double flux1 = 0.5*(cs*Mrf*(phi1*r1+phi2*r2)-cs*abs(Mrf)*(phi2*r2-phi1*r1))
    cdef double flux2 = copysign(0.5*(cs*Mrf*(phi1*r1*u1+phi2*r2*u2)-cs*abs(Mrf)*(phi2*r2*u2-phi1*r1*u1)) + pf, Mrf)
    cdef double flux3 = 0.5*(cs*Mrf*(phi1*r1*H1+phi2*r2*H2)-cs*abs(Mrf)*(phi2*r2*H2-phi1*r1*H1)) + pf*vbi

    return flux1, flux2, flux3


cpdef inline double get_extrapol_2e_v15(double v2, double v3) nogil:
    """cpdef inline double get_extrapol_2e_v15(double v2, double v3)
    функция получения экстраполированного значения 2 порядка (рис. 3 статьи Богданова-Миллера)

    v2 - значение в 2 ячейке
    v3 - значение в 3 ячейке
    return - экстраполированное значение на границе 1 и  2 ячейки
    """

    return 1.5 * v2 - 0.5 * v3

cpdef inline double get_extraopl_3e_v15(double v2, double v3, double v4) nogil:
    """cpdef inline double get_extraopl_3e_v15(double v2, double v3, double v4)
    функция получения экстраполированного значения 3 порядка (рис. 3 статьи Богданова-Миллера)

    v2 - значение в 2 ячейке
    v3 - значение в 3 ячейке
    v4 - значение в 4 ячейке
    return - экстраполированное значение на границе 1 и  2 ячейки
    """
    
    return 1.875*v2 - 1.25*v3 + 0.375*v4

cpdef inline double get_interp_3i_v15(double v1, double v2, double v3) nogil:
    """cpdef inline double get_interp_3i_v15(double v1, double v2, double v3)
    функция получения экстраполированного значения 3 порядка (рис. 3 статьи Богданова-Миллера)

    v1 - значение в 1 ячейке
    v2 - значение в 2 ячейке
    v3 - значение в 3 ячейке
    return - интерполированное значение на границе 1 и  2 ячейки
    """
    
    return 0.375*v1 + 0.75*v2 - 0.125*v3

cpdef double get_final_limited_value(double v_init, double v1, double v2, double delta_L=0.025) nogil:
    """cpdef double get_final_limited_value(double v_init, double v1, double v2, double delta_L=0.025)
    Функция получения финального занчения интер-экстраполированного значения (рис. 4 статьи Богданова-Миллера)

    v_init - изначальное занчение парамтера (ось Х)
    v1  - одно из соседних занчений
    v2  - другое из соседних занчений
    delta_L - значение для сглаживания
    """
    cdef double v_min = min2(v1, v2)
    cdef double v_max = max2(v1, v2)
    cdef double L = v_max - v_min
    cdef double delta = delta_L * L
    cdef double k0,k1,k2
    if v_init >= v_min + delta and v_init <= v_max - delta:
        return v_init
    if v_init >= v_max + delta:
        return v_max
    if v_init <= v_min - delta:
        return v_min
    if v_init > v_min-delta and v_init < v_min+delta:
        k0 = 0.25*delta + 0.5*v_min + 0.25*v_min*v_min/delta
        k1 = 0.5*(delta-v_min)/delta
        k2 = 0.25/delta
        return k0 + k1*v_init + k2*v_init*v_init
    else:
        k0 = -0.25*delta + v_max*0.5 - 0.25*v_max*v_max/delta
        k1 = 0.5*(delta+v_max)/delta
        k2 = -0.25/delta
        return k0 + k1*v_init + k2*v_init*v_init


cdef class GasEOS(object):
    """Класс для описания уравонений состояний газов и т.д. Поддерживает различные представления
    Абеля, и из книги Годунова

    Предоставляет методы определения внутр. энергии, давления и скорости звука
    self.kind == 1   - уравнение состояния Абеля
    self.kind == 2   - из книги Годунова
    """
    def __init__(self, gamma, kappa=0.0, p_0=0.0, c_0=0.0, kind=1):
        self.gamma = gamma
        self.kappa = kappa
        self.p_0 = p_0
        self.c_0 = c_0
        self.kind = kind

    def __repr__(self):
        return f'GasEOS(gamma={self.gamma}, kappa={self.kappa}, p_0={self.p0}, c_0={self.c_0}, kind={self.kind})'

    def __str__(self):
        return repr(self)

    cpdef double get_e(self, double ro, double p):
        if self.kind == 1:
            return rop_to_e(ro, p, self.gamma, self.kappa)
        else:
            return get_e_13_1(p, ro, self.p_0, self.c_0, self.gamma)

    cpdef double get_p(self, double ro, double e):
        if self.kind == 1:
            return roe_to_p(ro, e, self.gamma, self.kappa)
        else:
            return get_p_13_1(e, ro, self.p_0, self.c_0, self.gamma)

    cpdef double get_csound(self, double ro, double p):
        if self.kind == 1:
            return rop_to_csound(ro, p, self.gamma, self.kappa)
        else:
            return get_c_13_8(p, ro, self.p_0, self.gamma)


cdef class GasFluxCalculator(object):
    """Класс для расчета потоков через 2 соседние ячейки. Так же предоставляет методы для граничных условий

    Поддерживает расчет потоков по точному Римановскому солверу.
    self.flux_type == 1
    """
    
    def __init__(self, 
                flux_type=1, 
                left_border_type=1, 
                right_border_type=1, 
                x_order=1,
                n_iter_max = 17,
                epsF=1e-6,
                delta_L=0.025,
                alpha_1=8,
                alpha_2=4):
        """Конструктор класса
        
        Keyword Arguments:
            flux_type {int} -- тип солвера 1-Точный Римановский солвер (default: {1})
            left_border_type {int} -- Тип правой стенки 1-непроницаемая стенка (default: {1})
            right_border_type {int} -- -- (default: {1})
            x_order - как брать значения для распада разрыва (1-просто значения в соседних ячейках, 2-как в статье богданова-миллера)
        """

        self.set_flux_type(flux_type)
        self.left_border_type = left_border_type
        self.right_border_type = right_border_type
        self.x_order = x_order
        self.epsF=epsF
        self.n_iter_max = n_iter_max
        self.delta_L = delta_L
        self.alpha_1 = alpha_1
        self.alpha_2 = alpha_2

    cpdef void set_flux_type(self, int new_flux_type):
        if new_flux_type == 1:
            self.flux_type = new_flux_type
            self.rr_vals_len = 8
            self.rr_bint_len = 3
        else:
            print("Неправильный тип расчета потока!!! Будет использоваться 1")
            self.set_flux_type(1)

    cpdef void create_rr_arrs(self, GasLayer layer):
        layer.rr_vals = np.zeros((layer.n_cells+1, self.rr_vals_len), dtype=np.double)
        layer.rr_bint = np.zeros((layer.n_cells+1, self.rr_bint_len), dtype=np.int)


    cpdef void fill_rr_vals(self, GasLayer layer):
        if self.flux_type == 1:
            self.fill_rr_vals_borders_fill_fluxesURP(layer)
            self.fill_rr_vals_Riman(layer)

    cpdef void fill_rr_vals_borders_fill_fluxesURP(self, GasLayer layer):
        cdef double p_1, ro_1, u_1, c_1, p_2, ro_2, u_2, c_2,rU, rR, rP,ray_W
        cdef size_t i
        if self.flux_type == 1:
            if self.left_border_type == 1:  
                p_1 = layer.ps[0]
                ro_1 = layer.ros[0] 
                u_1 = layer.us[0] 
                c_1 = layer.cs[0] 
                ray_W = layer.Vs_borders[0]
                rU, rR, rP = border_wall_fill_rr(left_border=True, vbi=ray_W, p=p_1, ro=ro_1, u=u_1, c=c_1, \
                    p_0=layer.gasEOS.p_0, gamma=layer.gasEOS.gamma, 
                    rr_vals=layer.rr_vals[0], rr_bint=layer.rr_bint[0], 
                    eps_F=self.epsF, n_iter_max=self.n_iter_max)
                layer.fluxes[0, 0] = rU
                layer.fluxes[1, 0] = rR
                layer.fluxes[2, 0] = rP

                layer.D_left[0] = layer.rr_vals[0, 4]
                layer.U_kr[0] = layer.rr_vals[0, 2]

            if self.right_border_type == 1:
                i=layer.n_cells-1
                p_1 = layer.ps[i]
                ro_1 = layer.ros[i] 
                u_1 = layer.us[i] 
                c_1 = layer.cs[i] 
                ray_W = layer.Vs_borders[i+1]

                rU, rR, rP = border_wall_fill_rr(left_border=False, vbi=ray_W, p=p_1, ro=ro_1, u=u_1, c=c_1, \
                    p_0=layer.gasEOS.p_0, gamma=layer.gasEOS.gamma, 
                    rr_vals=layer.rr_vals[i+1], rr_bint=layer.rr_bint[i+1], 
                    eps_F=self.epsF, n_iter_max=self.n_iter_max)
                layer.fluxes[0, i+1] = rU
                layer.fluxes[1, i+1] = rR
                layer.fluxes[2, i+1] = rP

                layer.D_right[i] = layer.rr_vals[i+1, 0]
                layer.U_kr[i+1] = layer.rr_vals[i+1, 2]

    cpdef void fill_rr_vals_Riman(self, GasLayer layer):
        if layer.rr_vals.shape[1] != self.rr_vals_len or layer.rr_bint.shape[1] != self.rr_bint_len:
            self.create_rr_arrs(layer)
        cdef size_t i
        cdef double p_1, p_2, ro_1, ro_2, u_1, u_2, c_1, c_2, betta1, betta2, betta3
        if self.x_order == 1:
            for i in range(1, layer.n_cells):
                mega_foo_fill_rr(
                    p_1=layer.ps[i-1], ro_1=layer.ros[i-1], u_1=layer.us[i-1], c_1=layer.cs[i-1], \
                    p_2=layer.ps[i],   ro_2=layer.ros[i],   u_2=layer.us[i],   c_2=layer.cs[i], \
                    p_0=layer.gasEOS.p_0, gamma=layer.gasEOS.gamma, 
                    rr_vals=layer.rr_vals[i], rr_bint=layer.rr_bint[i], 
                    eps_F=self.epsF, n_iter_max=self.n_iter_max)

                layer.D_left[i] = layer.rr_vals[i, 4]
                layer.D_right[i-1] = layer.rr_vals[i, 0]
                layer.U_kr[i] = layer.rr_vals[i, 2]
        else:
            for i in range(1, layer.n_cells):
                c_1 = layer.cs[i-1]
                c_2 = layer.cs[i]
                ro_1 = layer.ros[i-1]
                ro_2 = layer.ros[i]
                betta1 = fabs(layer.us[i-1] - layer.us[i])/min2(c_1, c_2)
                betta2 = fabs(layer.ps[i-1] - layer.ps[i])/min2(ro_1*c_1*c_1, ro_2*c_2*c_2)
                betta3 = fabs(ro_1*c_1 - ro_2*c_2)/min2(ro_1*c_1, ro_2*c_2)
                layer.bettas[i] = max2(betta1, max2(betta2, betta3))

            for i in range(1, layer.n_cells):
                p_1 = self.get_v_left_forRiman(layer.ps, layer.bettas, i)
                p_2 = self.get_v_right_forRiman(layer.ps, layer.bettas, i)

                ro_1 = self.get_v_left_forRiman(layer.ros, layer.bettas, i)
                ro_2 = self.get_v_right_forRiman(layer.ros, layer.bettas, i)

                u_1 = self.get_v_left_forRiman(layer.us, layer.bettas, i)
                u_2 = self.get_v_right_forRiman(layer.us, layer.bettas, i)
                
                c_1 = self.get_v_left_forRiman(layer.cs, layer.bettas, i)
                c_2 = self.get_v_right_forRiman(layer.cs, layer.bettas, i)
                mega_foo_fill_rr(
                    p_1=p_1, ro_1=ro_1, u_1=u_1, c_1=c_1, \
                    p_2=p_2, ro_2=ro_2, u_2=u_2, c_2=c_2, \
                    p_0=layer.gasEOS.p_0, gamma=layer.gasEOS.gamma, 
                    rr_vals=layer.rr_vals[i], rr_bint=layer.rr_bint[i], 
                    eps_F=self.epsF, n_iter_max=self.n_iter_max)

                layer.D_left[i] = layer.rr_vals[i, 4]
                layer.D_right[i-1] = layer.rr_vals[i, 0]
                layer.U_kr[i] = layer.rr_vals[i, 2]

    cpdef double get_v_right_forRiman(self, double[:] vs, double[:] bettas, int i):
        """
        cpdef (double, double) get_v1_v2_forRiman(self, double[:] vs, size_t i)

        метод получения значений для задачи распада разрыва по методу Бгданова-Миллера) (рис. 3, 4)

        vs - массив со значениями
        i - индекс границы, где нужно найти значения слева и справа

        return v_right - значения для распада разрыва справа от границы i
               v_right - интер-экстраполированное значение vs[i] 
        """
        cdef double v1, v2, v3, v4, betta12, betta23, betta34, betta_max
        cdef int imax = vs.shape[0]-1
        # значения справа
        v1 = vs[i-1]
        betta12 = bettas[i]
        v2 = vs[i]
        if i+1>imax:
            v3 = vs[imax]
            betta23 = bettas[imax]
        else:
            v3 = vs[i+1]
            betta23 = bettas[i+1]
        if i+2>imax:
            v4 = vs[imax]
            betta34 = bettas[imax]
        else:
            v4 = vs[i+2]
            betta34 = bettas[i+2]
        betta_max = max2(betta12, max2(betta23, betta34))
        return self.get_bogdanov_v(v1, v2, v3, v4, betta_max)


    cpdef double get_v_left_forRiman(self, double[:] vs, double[:] bettas, int i):
        """
        cpdef double get_v_left_forRiman(self, double[:] vs, int i)

        метод получения значений для задачи распада разрыва по методу Бгданова-Миллера) (рис. 3, 4)

        vs - массив со значениями
        i - индекс границы, где нужно найти значения слева и справа

        return v_left - значение для распада разрыва слева от границы i.
               v_left - интер-экстраполированное значение vs[i-1]
        """
        # значения справа
        cdef double v1, v2, v3, v4, v_left, betta12, betta23, betta34, betta_max
        v1 = vs[i]
        betta12 = bettas[i]
        v2 = vs[i-1]
        if i-2<0:
            v3 = vs[0]
            betta23 = bettas[1]
        else:
            v3 = vs[i-2]
            betta23 = bettas[i-1]
        if i-3<0:
            v4 = vs[0]
            betta34 = bettas[1]
        else:
            v4 = vs[i-3]
            betta34 = bettas[i-2]
        betta_max = max2(betta12, max2(betta23, betta34))
        return self.get_bogdanov_v(v1, v2, v3, v4, betta_max)
        

    cpdef double get_bogdanov_v(self, double v1, double v2, double v3, double v4, double betta_max):
        if betta_max >= self.alpha_1:
            return v2
        cdef double res = self.get_interextra_value(v1, v2, v3, v4)
        res = get_final_limited_value(res, v1, v2, self.delta_L)
        if betta_max < self.alpha_2:
            return res
        cdef double t = (betta_max-self.alpha_2)/(self.alpha_1-self.alpha_2)
        return t*res + (1-t)*v2

    cpdef double get_interextra_value(self, double v1, double v2, double v3, double v4):
        """
        cpdef double get_interextra_value(self, double v1, double v2, double v3, double v4)

        Получить интерполнированное значение для задачи распада разрыва по методу Богданова-Миллера (рис. 3)
        между ячейками 1 и 2 для ячейки 2
        v1 v2 v3 v4  - значения в ячейках
        return - интер-экстраполированное значение (справа) для v2 ячейки
        """

        cdef double e2, e3, i3, e_min, e_max
        e2 = get_extrapol_2e_v15(v2, v3)
        e3 = get_extraopl_3e_v15(v2, v3, v4)
        i3 = get_interp_3i_v15(v1, v2, v3)
        e_min = min2(e2, e3)
        e_max = max2(e2, e3)
        if i3 >= e_min and i3 <= e_max:
            return i3
        if i3 > e_max:
            return e_max
        if i3 < e_min:
            return e_min

    cpdef void fill_fluxesURP(self, GasLayer layer):
        if self.flux_type == 1:
            self.fill_fluxesURP_Godunov(layer)

    cpdef void fill_fluxesURP_Godunov(self, GasLayer layer):
        cdef size_t i
        for i in range(1, layer.n_cells):
            layer.fluxes[0, i], layer.fluxes[1, i], layer.fluxes[2, i] = get_ray_URP(
                ray_W=layer.Vs_borders[i], UD_left=layer.rr_bint[i, 1], UD_right=layer.rr_bint[i, 2], \
                D_1=layer.rr_vals[i, 0], D_star_1=layer.rr_vals[i, 1], U=layer.rr_vals[i, 2], \
                D_star_2=layer.rr_vals[i, 3], D_2=layer.rr_vals[i, 4], R_1=layer.rr_vals[i, 5], R_2=layer.rr_vals[i, 6], \
                P=layer.rr_vals[i, 7], \
                p_1=layer.ps[i-1], ro_1=layer.ros[i-1], u_1=layer.us[i-1], c_1=layer.cs[i-1], \
                p_2=layer.ps[i], ro_2=layer.ros[i], u_2=layer.us[i], c_2=layer.cs[i], gamma=layer.gasEOS.gamma)


cdef class GridStrecher(object):
    """Класс дял управления координатами и скоростями одномерных сеток

    self.strech_type == 1 - обычная равномерная сетка
                     == 2 - адаптивная сетка (еще не совсем гуд)

    st2_window_part, st2_adapt_prop - параметры для адаптивной сетки
    """

    def __init__(self, strech_type=1, st2_window_part=0.01, st2_adapt_prop=1, D_mnj=0.99):
        self.strech_type = strech_type
        self.st2_window_part = st2_window_part
        self.st2_adapt_prop = st2_adapt_prop
        self.bufarr = np.zeros(7)
        self.bufarr_border = np.zeros(8)
        self.interp_smooth = InterpXY([1,2,3], [1,2,3])
        self.interp_adapt = InterpXY([1,2,3,4], [1,2,3,4])
        self.D_mnj =D_mnj

    cpdef bint evaluate(self, double tau, GasLayer layer):
        cdef size_t i
        for i in range(layer.xs_borders.shape[0]):
            layer.xs_borders[i] = layer.xs_borders[i] + tau * layer.Vs_borders[i]
            if i > 0 and layer.xs_borders[i]-layer.xs_borders[i-1] <= 1e-12:
                return False
        self.fill_xs_cells(layer.xs_borders, layer.xs_cells)
        return True

    cpdef void init_regular(self, double v1, double v2, double[:] vs):
        cdef int n = vs.shape[0]
        cdef double dv = (v2-v1)/(n-1)
        cdef size_t i
        for i in range(vs.shape[0]):
            vs[i] = v1 + i*dv

    cpdef bint sync_layers(self, GasLayer layer0, GasLayer layer1, double tau, double v_left, double v_right):
        layer1.time = layer0.time + tau
        # self.fill_Vs_borders_proportional(v_left, v_right, layer0.xs_borders, layer0.Vs_borders)
        self.fill_Vs_borders_proportional(v_left, v_right, layer1.xs_borders, layer1.Vs_borders)
        cdef bint suc = self.evaluate(tau, layer1)
        if not suc:
            return False
        if self.strech_type == 2:
            # TODO Адаптивная сетка аццтой
            self.strech_layer_adaptive(layer1, layer_prev=layer0, tau=tau)
        cdef size_t i
        for i in range(layer0.Vs_borders.shape[0]):
            layer0.Vs_borders[i] = (layer1.xs_borders[i] - layer0.xs_borders[i])/tau
        layer1.Vs_borders[:] = layer0.Vs_borders
        return True

    cpdef void strech_layer_adaptive(self, GasLayer layer, GasLayer layer_prev, double tau):
        if self.bufarr.shape[0] != layer.xs_cells.shape[0]:
            self.bufarr = np.zeros(layer.xs_cells.shape[0])
        if self.bufarr_border.shape[0] != layer.xs_borders.shape[0]:
            self.bufarr_border = np.zeros(layer.xs_borders.shape[0])
        self.smooth_arr(layer.xs_cells, layer.ros, self.bufarr, self.st2_window_part)
        self.adaptine_borders(layer.xs_borders, self.bufarr, self.bufarr_border)
        cdef size_t i
        cdef double vi, D_1, D_2, mnj, v_min, v_max, xtmp, tmp1, t
        mnj = self.D_mnj
        for i in range(1, layer.xs_borders.shape[0]-1):
            vi = (self.bufarr_border[i] - layer_prev.xs_borders[i])/tau
            D_1 = layer_prev.rr_vals[i, 0]
            D_2 = layer_prev.rr_vals[i, 4]

            t = 0.9
            xtmp = (1-t)*layer_prev.xs_borders[i-1] + t*layer_prev.xs_borders[i]
            v_min = (xtmp - layer_prev.xs_borders[i])/tau

            t = 0.9
            xtmp = (1-t)*layer_prev.xs_borders[i+1] + t*layer_prev.xs_borders[i]
            v_max = (xtmp - layer_prev.xs_borders[i])/tau

            v_min = max2(v_min, D_1)
            v_max = min2(v_max, D_2)

            # v_min = D_1
            # v_max = D_2

            if vi < mnj*v_min:
                vi = mnj*v_min
            elif vi > mnj*v_max:
                vi = mnj*v_max
            self.bufarr_border[i] = layer_prev.xs_borders[i] + vi*tau
        layer.xs_borders[:] =  self.bufarr_border
        self.fill_xs_cells(layer.xs_borders, layer.xs_cells)


    cpdef void fill_Vs_borders_proportional(self, double v_left, double v_right, double[:] xs_borders, double[:] Vs_borders):
        cdef size_t i
        cdef double dV = v_right - v_left
        cdef double mnj = dV/(xs_borders[xs_borders.shape[0]-1]-xs_borders[0])
        Vs_borders[0] = v_left
        for i in range(1, Vs_borders.shape[0]-1):
            Vs_borders[i] = v_left + mnj * (xs_borders[i]-xs_borders[0])
        Vs_borders[Vs_borders.shape[0]-1] = v_right

    cpdef void fill_euler_vel0_regular(self, double tau, double[:] xs0, double[:] xs1, double[:] vel0_fill):
        cdef size_t i
        for i in range(vel0_fill.shape[0]):
            vel0_fill[i] = (xs1[i] - xs0[i])/tau

    cpdef void fill_xs_cells(self, double[:] xs_borders, double[:] xs_fill):
        cdef size_t i
        for i in range(xs_fill.shape[0]):
            xs_fill[i] = (xs_borders[i+1] + xs_borders[i])/2

    cpdef void smooth_arr(self, double[:] xs, double[:] vs, double[:] vs_smoothed, double window_part=0.1):
        if self.interp_smooth.length != vs.shape[0]:
            self.interp_smooth.set_length(vs.shape[0])
        
        vs_smoothed[:] = vs
        return

        cdef size_t i
        self.interp_smooth.xs[:] = xs 
        self.interp_smooth.ys[0] = 0 
        for i in range(1, xs.shape[0]):
            self.interp_smooth.ys[i] = self.interp_smooth.ys[i-1] + (xs[i] - xs[i-1])*0.5*(vs[i]+vs[i-1])     
        self.interp_smooth.sync_ks_bs()

        cdef double xsmax = xs[xs.shape[0]-1]
        cdef double xsmin = xs[0]
        cdef double window = window_part * (xs[xs.shape[0]-1] - xs[0])
        cdef double left, right
        for i in range(vs_smoothed.shape[0]):
            left = max2(xs[i] - window, xsmin)
            right = min2(xs[i] + window, xsmax)
            vs_smoothed[i] = (self.interp_smooth.get_v(right) - self.interp_smooth.get_v(left))/(right-left)
        

    cpdef void adaptine_borders(self, double[:] xs_borders, double[:] vs, double[:] xs_adapt):
        if self.interp_adapt.length != xs_borders.shape[0]:
            self.interp_adapt.set_length(xs_borders.shape[0])
        cdef size_t i
        cdef double xl, xr, yl, yr, li
        self.interp_adapt.xs[0] = 0
        self.interp_adapt.ys[:] = xs_borders

        
        cdef double v_min = vs[0]
        cdef double v_max = vs[0]
        for i in range(1, vs.shape[0]):
            if vs[i] > v_max:
                v_max = vs[i]
            elif vs[i]<v_min:
                v_min = vs[i]

        cdef double delta_x = xs_borders[xs_borders.shape[0]-1]-xs_borders[0]
        cdef double mnj = delta_x/(v_max-v_min)

        for i in range(1, xs_borders.shape[0]-1):
            xl = xs_borders[i-1]
            yl = vs[i-1]*mnj
            xr = xs_borders[i]
            yr = vs[i]*mnj
            li = sqrt((xr-xl)**2+(yr-yl)**2)
            self.interp_adapt.xs[i] = self.interp_adapt.xs[i-1] + li
        i=xs_borders.shape[0]-1
        self.interp_adapt.xs[i] = self.interp_adapt.xs[i-1] + (xs_borders[i]-xs_borders[i-1])
        self.interp_adapt.sync_ks_bs()

        li = (self.interp_adapt.xs[i] - self.interp_adapt.xs[0])/(xs_borders.shape[0]-1)
        for i in range(xs_adapt.shape[0]):
            xs_adapt[i] = self.interp_adapt.get_v(i*li)


cdef class GasLayer(object):
    def __init__(self, n_cells, Tube tube, GasEOS gasEOS, GasFluxCalculator flux_calculator, GridStrecher grid_strecher, int n_qs=3):
        self.n_cells = n_cells
        self.n_qs = n_qs
        # self.n_pars = n_pars
        self.tube = tube
        self.gasEOS = gasEOS
        self.time = 0
        self.flux_calculator = flux_calculator
        self.grid_strecher = grid_strecher
        
        self.xs_cells = np.zeros(n_cells, dtype=np.double)
        self.xs_borders = np.zeros(n_cells+1, dtype=np.double)
        self.Vs_borders = np.zeros(n_cells+1, dtype=np.double)

        self.S = np.zeros(n_cells+1, dtype=np.double)

        self.ds = np.zeros(n_cells, dtype=np.double)
        self.W = np.zeros(n_cells, dtype=np.double)

        self.ros = np.zeros(n_cells, dtype=np.double)
        self.ps = np.zeros(n_cells, dtype=np.double)
        self.us = np.zeros(n_cells, dtype=np.double)
        self.es = np.zeros(n_cells, dtype=np.double)
        self.cs = np.zeros(n_cells, dtype=np.double)

        # self.pars = np.zeros((n_pars, n_cells), dtype=np.double) # 0 - ros, 1 - ps, 2 - us 

        self.taus = np.zeros(n_cells, dtype=np.double)
        self.D_left = np.zeros(n_cells, dtype=np.double)
        self.D_right = np.zeros(n_cells, dtype=np.double)
        self.U_kr = np.zeros(n_cells+1, dtype=np.double)

        self.fluxes = np.zeros((n_qs, n_cells+1), dtype=np.double)

        self.qs = np.zeros((n_qs, n_cells), dtype=np.double)

        self.hs = np.zeros((n_qs, n_cells), dtype=np.double)
        # self.mega_foo_results = np.zeros((n_cells+1, ), dtype=np.double)
        self.flux_calculator.create_rr_arrs(self)

        self.bettas = np.zeros(n_cells+1, dtype=np.double)

    cpdef void copy_params_to(self, GasLayer to_me):
        to_me.time = self.time

        to_me.xs_cells[:] = self.xs_cells
        to_me.xs_borders[:] = self.xs_borders
        to_me.Vs_borders[:] = self.Vs_borders

        to_me.S[:] = self.S 

        to_me.ds[:] = self.ds
        to_me.W[:] = self.W

        to_me.ps[:] = self.ps
        to_me.ros[:] = self.ros
        to_me.us[:] = self.us
        to_me.es[:] = self.es
        to_me.cs[:] = self.cs

        # to_me.pars[:] = self.pars

        to_me.taus[:] = self.taus
        to_me.D_left[:] = self.D_left
        to_me.D_right[:] = self.D_right
        to_me.U_kr[:] = self.U_kr

        to_me.fluxes[:] = self.fluxes
        to_me.qs[:] = self.qs
        to_me.hs[:] = self.hs   

        to_me.rr_bint[:] = self.rr_bint
        to_me.rr_vals[:] = self.rr_vals   
        to_me.bettas[:] = self.bettas  

    cpdef GasLayer copy(self):
        cdef GasLayer res = GasLayer(self.n_cells, self.tube, self.gasEOS, self.flux_calculator, self.grid_strecher, self.n_qs)
        self.copy_params_to(res)
        return res


    cpdef void init_ropue_fromfoo(self, foo_ropu, bint init_q=True,  bint init_SsdW=True):
        self.grid_strecher.fill_xs_cells(self.xs_borders, self.xs_cells)
        cdef size_t i
        cdef double x
        cdef GasEOS eos = self.gasEOS
        cdef double ro, p, u, e
        for i in range(self.n_cells):
            ro, p, u = foo_ropu(self.xs_cells[i], self)
            self.ros[i], self.ps[i], self.us[i] = ro, p, u
            self.es[i] = eos.get_e(ro, p)
            self.cs[i] = eos.get_csound(ro, p)
        if init_q:
            self.init_q()
        if init_SsdW:
            self.init_SsdW()
        self.init_taus_acustic()
    
    cpdef void init_SsdW(self):
        self.tube.fill_S(self.xs_borders, self.S)
        self.tube.fill_dsdx(self.xs_borders, self.ds)
        self.tube.fill_W(self.xs_borders, self.W)
    
    cpdef void init_q(self):
        cdef size_t i
        cdef size_t n = self.ros.shape[0]
        for i in range(n):
            q1, q2, q3 = roue_to_q(self.ros[i], self.us[i], self.es[i])
            self.qs[0, i] = q1
            self.qs[1, i] = q2 
            self.qs[2, i] = q3

    cpdef void init_ropue(self):
        cdef size_t i
        cdef size_t n = self.ros.shape[0]
        cdef double ro, p, u, e
        for i in range(n):
            ro, u, e = q123_to_roue(self.qs[0, i], self.qs[1, i], self.qs[2, i])
            p = self.gasEOS.get_p(ro, e)
            self.ros[i] = ro
            self.ps[i] = p
            self.us[i] = u
            self.es[i] = e
            self.cs[i] = self.gasEOS.get_csound(ro, p)
        
    cpdef void init_h(self):
        cdef size_t i
        for i in range(self.hs.shape[1]):
            self.hs[0, i] = 0
            self.hs[1, i] = self.ps[i] * self.ds[i] 
            self.hs[2, i] = 0

    cpdef double get_tau_min(self):
        cdef size_t i
        cdef double tau_min = self.taus[0]
        for i in range(1, self.taus.shape[0]):
            if self.taus[i] < tau_min:
                tau_min = self.taus[i]
        return tau_min

    cpdef void init_taus_acustic(self):
        cdef size_t i
        for i in range(self.taus.shape[0]):
            self.taus[i] = (self.xs_borders[i+1] - self.xs_borders[i])/(abs(self.us[i])+self.cs[i])

    cpdef void fill_fluxesURP(self):
        self.flux_calculator.fill_fluxesURP(self)

    cpdef void fill_fluxes(self):
        cdef size_t i
        cdef double rayE, rayU, rayR, rayP, M_j, J_j, E_j, ray_W
        for i in range(self.fluxes.shape[1]):
            ray_W = self.Vs_borders[i]
            rayU = self.fluxes[0, i]
            rayR = self.fluxes[1, i]
            rayP = self.fluxes[2, i]
            rayE = self.gasEOS.get_e(rayR, rayP)

            M_j = rayR * (rayU - ray_W)
            J_j = rayP + M_j * rayU
            E_j = (rayE + 0.5*rayU*rayU)*M_j + rayP*rayU
            self.fluxes[0, i] = M_j
            self.fluxes[1, i] = J_j
            self.fluxes[2, i] = E_j

    cpdef void fill_rr(self):
        self.flux_calculator.fill_rr_vals(self)
    
    cpdef void fill_taus(self):
        cdef size_t i
        cdef double dx, v_otn_left, v_otn_right, tau_left, tau_right
        for i in range(self.taus.shape[0]):
            dx = self.xs_borders[i+1] - self.xs_borders[i]
            v_otn_left = self.D_left[i] - self.Vs_borders[i+1]
            if v_otn_left < 1e-12:
                tau_left = 999999
            else:
                tau_left = dx / v_otn_left
            
            v_otn_right = self.Vs_borders[i] - self.D_right[i] 
            if v_otn_right < 1e-12:
                tau_right = 999999
            else:
                tau_right = dx / v_otn_right
            if tau_left < tau_right:
                self.taus[i] = tau_left
            else:
                self.taus[i] = tau_right

    cpdef GasLayer step_simple(self, double tau, double v_left, double v_right):
        self.Vs_borders[0] = v_left
        self.Vs_borders[self.n_cells] = v_right
        self.fill_rr()
        cpdef GasLayer layer1 = self.copy() 
        cdef bint suc = self.grid_strecher.sync_layers(self, layer1, tau, v_left, v_right)
        if not suc:
            print('suc 454')
            return self
        self.init_h()
        self.fill_taus()
        layer1.init_SsdW()
        layer1.taus[:] = self.taus
        if self.get_tau_min() < tau:
            print('suc 455')
            return self
        cdef size_t i, j
        cdef double dx
        self.fill_fluxesURP()
        self.fill_fluxes()
        for i in range(layer1.xs_cells.shape[0]):
            dx = self.xs_borders[i+1] - self.xs_borders[i]
            for j in range(layer1.n_qs):
                layer1.qs[j, i] = (self.qs[j, i]*self.W[i] - tau*(self.S[i+1]*self.fluxes[j, i+1] - self.S[i]*self.fluxes[j, i]) + tau*self.hs[j, i]*dx)/layer1.W[i]
        layer1.init_ropue()
        return layer1

    cpdef GasLayer corrector(self, GasLayer lr_simple05, double tau, double v_left, double v_right):
        lr_simple05.Vs_borders[0] = v_left
        lr_simple05.Vs_borders[self.n_cells] = v_right
        lr_simple05.fill_rr()
        cpdef GasLayer layer2 = lr_simple05.copy()
        cdef double tau2 = self.time+tau-lr_simple05.time
        cdef bint suc = self.grid_strecher.sync_layers(lr_simple05, layer2, tau2, v_left, v_right)
        if not suc:
            print('suc 354')
            return self
        lr_simple05.init_h()
        lr_simple05.fill_taus()
        layer2.init_SsdW()
        layer2.taus[:] = lr_simple05.taus
        if lr_simple05.get_tau_min() < tau:
            print('suc 355')
            return self
        cdef size_t i, j
        cdef double dx
        lr_simple05.fill_fluxesURP()
        lr_simple05.fill_fluxes()
        for i in range(layer2.xs_cells.shape[0]):
            dx = lr_simple05.xs_borders[i+1] - lr_simple05.xs_borders[i]
            for j in range(layer2.n_qs):
                layer2.qs[j, i] = (
                    self.qs[j, i]*self.W[i] - tau*(lr_simple05.S[i+1]*lr_simple05.fluxes[j, i+1] - lr_simple05.S[i]*lr_simple05.fluxes[j, i]) 
                    + tau*lr_simple05.hs[j, i]*dx
                    )/layer2.W[i]
        layer2.init_ropue()
        return layer2

    cpdef GasLayer corrector_bogdanov(self, GasLayer lr_simple, double v_left, double v_right):
        lr_simple.fill_rr()
        cpdef GasLayer layer2 = self.copy()
        cdef double tau = lr_simple.time - self.time
        cdef double vl = 0.5*(v_left + lr_simple.Vs_borders[0])
        cdef double vr = 0.5*(v_right + lr_simple.Vs_borders[lr_simple.n_cells])
        cdef bint suc = self.grid_strecher.sync_layers(self, layer2, tau, vl, vr)
        if not suc:
            print('suc 254')
            return self
        lr_simple.Vs_borders[:] = layer2.Vs_borders
        lr_simple.init_h()
        lr_simple.fill_taus()
        layer2.init_SsdW()
        layer2.taus[:] = lr_simple.taus
        if lr_simple.get_tau_min() < tau:
            print('suc 255')
            return self
        cdef size_t i, j
        cdef double dx, dx2
        lr_simple.fill_fluxesURP()
        lr_simple.fill_fluxes()
        for i in range(layer2.xs_cells.shape[0]):
            dx = lr_simple.xs_borders[i+1] - lr_simple.xs_borders[i]
            dx2 = self.xs_borders[i+1] - self.xs_borders[i]
            for j in range(layer2.n_qs):
                layer2.qs[j, i] = (
                    self.qs[j, i]*self.W[i] 
                    - tau * (0.25*(lr_simple.S[i+1]+self.S[i+1])*(lr_simple.fluxes[j, i+1]+self.fluxes[j, i+1]) 
                           - 0.25*(lr_simple.S[i]+self.S[i])*(lr_simple.fluxes[j, i]+self.fluxes[j, i])) 
                    + tau * 0.5*(lr_simple.hs[j, i]*dx + self.hs[j, i]*dx2)
                    )/layer2.W[i]
        layer2.init_ropue()
        return layer2



    def to_dict(self):
        return {
            'time': self.time,
            'n_qs': self.n_qs,
            'xs_cells': np.array(self.xs_cells),
            'xs_borders': np.array(self.xs_borders),
            'Vs_borders': np.array(self.Vs_borders),
            'ps': np.array(self.ps),
            'ros': np.array(self.ros),
            'us': np.array(self.us),
            'es': np.array(self.es),
            'cs': np.array(self.cs),
            'taus': np.array(self.taus),
            'D_left': np.array(self.D_left),
            'D_right': np.array(self.D_right),
            'U_kr': np.array(self.U_kr),
            'ds': np.array(self.ds),
            'S': np.array(self.S),
            'W': np.array(self.W),

            'qs':np.array(self.qs),
            'hs':np.array(self.hs),
            'fluxes':np.array(self.fluxes),
            'rr_vals':np.array(self.rr_vals),
            'rr_bint':np.array(self.rr_bint),
            'bettas':np.array(self.bettas)
       }
    def from_dict(self, d):
        self.time = np.array(d['time'])
        self.n_cells = np.array(d['xs_cells'].shape[0])
        self.n_qs = np.array(d['n_qs'])
        self.xs_cells = np.array(d['xs_cells'])
        self.xs_borders = np.array(d['xs_borders'])
        self.Vs_borders = np.array(d['Vs_borders'])

        self.S = np.array(d['S'] )

        self.ds = np.array(d['ds'])
        self.W = np.array(d['W'])

        self.ps = np.array(d['ps'])
        self.ros = np.array(d['ros'])
        self.us = np.array(d['us'])
        self.es = np.array(d['es'])
        self.cs = np.array(d['cs'])

        self.taus = np.array(d['taus'])
        self.D_left = np.array(d['D_left'])
        self.D_right = np.array(d['D_right'])
        self.U_kr = np.array(d['U_kr'],)

        self.fluxes = np.array(d['fluxes'])
        self.qs = np.array(d['qs'])
        self.hs = np.array(d['hs'])
        self.rr_vals = np.array(d['rr_vals'])
        self.rr_bint = np.array(d['rr_bint'])
        self.bettas = np.array(d['bettas'])

    def __str__(self):
        dxs = np.asarray(self.xs_borders)[1:] - np.asarray(self.xs_borders)[:-1]
        return f"{self.__class__.__name__}(n_cells={self.n_cells}); {{'p_max':{np.max(self.ps)}, 'tau_min': {self.get_tau_min()}, \
        'u_max': {np.max(self.us)}, 'cs_max': {np.max(self.cs)}, 'dx_min': {np.min(dxs)}, 'x_1': {self.xs_borders[0]}, 'x_2': {self.xs_borders[-1]},\
        'V_1':  {self.Vs_borders[0]}, 'V_2':  {self.Vs_borders[-1]}  }}"

    def __repr__(self):
        return str(self)
    
    




