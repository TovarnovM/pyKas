# distutils: language=c++
# cython: language_level=3, boundscheck=False, nonecheck=False, cdivision=True, initializedcheck=False

from tube cimport Tube, InterpXY
import cython
from libc.math cimport pi, sqrt, abs, copysign, abs
from godunov cimport get_e_13_1, get_p_13_1, get_c_13_8, mega_foo_cython, \
                     get_ray_URP, MegaFooResult, border_wall_URPDD_result, Border_URPDD_Result
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

cpdef inline void roue_to_q_(
    double ro,
    double u,
    double e,
    double[:] q):

    q[0] = ro
    q[1] = ro * u
    q[2] = ro * (e + 0.5 * u * u)

cpdef inline (double, double, double) roue_to_q(
    double ro,
    double u,
    double e):

    return ro, ro * u, ro * (e + 0.5 * u * u)


cpdef inline (double, double, double) q_to_roue(double[:] q):
    cdef double ro = q[0]
    cdef double u = q[1] / q[0]
    cdef double e = q[2] / q[0] - 0.5 * u * u
    return ro, u, e

cpdef inline (double, double, double) q123_to_roue(double q1, double q2, double q3):
    cdef double ro = q1
    cdef double u = q2 / q1
    cdef double e = q3 / q1 - 0.5 * u * u
    return ro, u, e


cpdef inline double roe_to_p(
    double ro,
    double e,
    double gamma,
    double b):

    return (gamma-1)*e*ro/(1-b*ro)



cpdef inline double rop_to_e(
    double ro,
    double p,
    double gamma,
    double b):

    return p*(1/ro-b)/(gamma-1)


cpdef inline double rop_to_csound(
    double ro,
    double p,
    double gamma,
    double b):

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
    double vbi):

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


cdef class GasEOS(object):
    def __init__(self, gamma, kappa=0.0, p_0=0.0, c_0=0.0, kind=1):
        self.gamma = gamma
        self.kappa = kappa
        self.p_0 = p_0
        self.c_0 = c_0
        self.kind = kind

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
    def __init__(self, flux_type=1, left_border_type=1, right_border_type=1):
        self.flux_type = flux_type
        self.left_border_type = left_border_type
        self.right_border_type = right_border_type
        self.epsF=1e-6
        self.n_iter_max = 7

    cpdef void fill_fluxes_taus(self, GasLayer layer):
        if self.flux_type == 1:
            self.fill_fluxes_Ds_Godunov(layer.Vs_borders, layer.ps, layer.ros, layer.us, layer.cs, layer.gasEOS,\
                layer.fluxes[0], layer.fluxes[1], layer.fluxes[2], layer.D_left, layer.D_right)
            self.fill_taus_Godunov(layer.xs_borders, layer.Vs_borders, layer.D_left, layer.D_right, layer.taus)

    cpdef void fill_taus_Godunov(self, double[:] xs_borders, double[:] Vs_borders, double[:] D_left, double[:] D_right, double[:] taus):
        cdef size_t i
        cdef double dx, v_otn_left, v_otn_right, tau_left, tau_right
        for i in range(taus.shape[0]):
            dx = xs_borders[i+1] - xs_borders[i]
            v_otn_left = D_left[i] - Vs_borders[i+1]
            if v_otn_left < 1e-12:
                tau_left = 999
            else:
                tau_left = dx / v_otn_left
            
            v_otn_right = Vs_borders[i] - D_right[i] 
            if v_otn_right < 1e-12:
                tau_right = 999
            else:
                tau_right = dx / v_otn_right
            if tau_left < tau_right:
                taus[i] = tau_left
            else:
                taus[i] = tau_right
            

    cpdef void fill_fluxes_Ds_Godunov(self, double[:] Vs_borders, double[:] ps, double[:] ros, \
            double[:] us, double[:] cs, GasEOS gasEOS, \
            double[:] flux1, double[:] flux2, double[:] flux3, double[:] D_left, double[:] D_right):
        cdef size_t i
        cdef double p_1, ro_1, u_1, c_1, p_2, ro_2, u_2, c_2
        cdef MegaFooResult mfres
        cdef double rayU, rayR, rayP, rayE, ray_W, M_j, J_j, E_j, D_1, D_2
        cdef Border_URPDD_Result border_res

        p_1 = ps[0]
        ro_1 = ros[0] 
        u_1 = us[0] 
        c_1 = cs[0] 
        ray_W = Vs_borders[0]
        
        border_res = border_wall_URPDD_result(True, ray_W, p_1, ro_1, u_1, c_1, gasEOS.p_0, gasEOS.gamma, self.epsF, self.n_iter_max)
        rayU = border_res.rU 
        rayR=border_res.rR 
        rayP=border_res.rP 
        D_1 = border_res.D_1 
        D_2 = border_res.D_2 
        
        rayE = gasEOS.get_e(rayR, rayP)

        M_j = rayR * (rayU - ray_W)
        J_j = rayP + M_j * rayU
        E_j = (rayE + 0.5*rayU*rayU)*M_j + rayP*rayU

        flux1[0] = M_j
        flux2[0] = J_j
        flux3[0] = E_j
        D_left[0] = D_2
        for i in range(1, ps.shape[0]):
            p_1 = ps[i-1]
            ro_1 = ros[i-1] 
            u_1 = us[i-1] 
            c_1 = cs[i-1] 
            p_2 = ps[i] 
            ro_2 = ros[i]
            u_2 = us[i]
            c_2 = cs[i]
            ray_W = Vs_borders[i]
            mfres = mega_foo_cython(p_1, ro_1, u_1, c_1, p_2, ro_2, u_2, c_2, gasEOS.p_0, gasEOS.gamma)
            rayU, rayR, rayP = get_ray_URP(ray_W, mfres.UD_left, mfres.UD_right, 
                                        mfres.D_1, mfres.D_star_1, mfres.U, 
                                        mfres.D_star_2, mfres.D_2, mfres.R_1, mfres.R_2, mfres.P, \
                                        p_1, ro_1, u_1, c_1, \
                                        p_2, ro_2, u_2, c_2, gasEOS.gamma)
            rayE = gasEOS.get_e(rayR, rayP)

            M_j = rayR * (rayU - ray_W)
            J_j = rayP + M_j * rayU
            E_j = (rayE + 0.5*rayU*rayU)*M_j + rayP*rayU

            flux1[i] = M_j
            flux2[i] = J_j
            flux3[i] = E_j
            D_right[i-1] = D_1
            D_left[i] = D_2 

        i=ps.shape[0]-1
        ray_W = Vs_borders[i+1]
        p_1 = ps[i]
        ro_1 = ros[i] 
        u_1 = us[i] 
        c_1 = cs[i] 

        border_res = border_wall_URPDD_result(False, ray_W, p_1, ro_1, u_1, c_1, gasEOS.p_0, gasEOS.gamma, self.epsF, self.n_iter_max)
        
        rayU = border_res.rU 
        rayR=border_res.rR 
        rayP=border_res.rP 
        D_1 = border_res.D_1 
        D_2 = border_res.D_2 

        rayE = gasEOS.get_e(rayR, rayP)
        M_j = rayR * (rayU - ray_W)
        J_j = rayP + M_j * rayU
        E_j = (rayE + 0.5*rayU*rayU)*M_j + rayP*rayU

        flux1[i+1] = M_j
        flux2[i+1] = J_j
        flux3[i+1] = E_j
        D_right[i] = D_1

cpdef inline double min2(double a, double b):
    if a < b:
        return a
    else:
        return b

cpdef inline double max2(double a, double b):
    if a > b:
        return a
    else:
        return b

cdef class GridStrecher(object):
    def __init__(self, strech_type=1, st2_window_part=0.1, st2_adapt_prop=1):
        self.strech_type = strech_type
        self.st2_window_part = st2_window_part
        self.st2_adapt_prop = st2_adapt_prop
        self.bufarr = np.zeros(7)
        self.bufarr_border = np.zeros(8)
        self.interp_smooth = InterpXY([1,2,3], [1,2,3])
        self.interp_adapt = InterpXY([1,2,3,4], [1,2,3,4])

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
            self.strech_layer_adaptive(layer1, layer_prev=layer0, tau=tau)
        cdef size_t i
        for i in range(layer0.Vs_borders.shape[0]):
            layer0.Vs_borders[i] = (layer1.xs_borders[i] - layer0.xs_borders[i])/tau
        return True

    cpdef void strech_layer_adaptive(self, GasLayer layer, GasLayer layer_prev, double tau):
        if self.bufarr.shape[0] != layer.xs_cells.shape[0]:
            self.bufarr = np.zeros(layer.xs_cells.shape[0])
        if self.bufarr_border.shape[0] != layer.xs_borders.shape[0]:
            self.bufarr_border = np.zeros(layer.xs_borders.shape[0])
        self.smooth_arr(layer.xs_cells, layer.es, self.bufarr, self.st2_window_part)
        self.adaptine_borders(layer.xs_borders, self.bufarr, self.bufarr_border)
        cdef size_t i
        cdef double old_x, adapt_x, limit
        for i in range(1, layer.xs_borders.shape[0]-1):
            old_x = layer.xs_borders[i]
            adapt_x = self.bufarr_border[i]
            if adapt_x > old_x:
                limit = (min2(layer.xs_borders[i+1], layer_prev.xs_borders[i+1])-old_x)*self.st2_adapt_prop
                if adapt_x > limit:
                    adapt_x = limit+old_x
            else: 
                limit = (-max2(layer.xs_borders[i-1], layer_prev.xs_borders[i-1])+old_x)*self.st2_adapt_prop
                if adapt_x < limit:
                    adapt_x = limit+old_x
            
            layer.xs_borders[i] = old_x + adapt_x-old_x
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

        self.ps = np.zeros(n_cells, dtype=np.double)
        self.ros = np.zeros(n_cells, dtype=np.double)
        self.us = np.zeros(n_cells, dtype=np.double)
        self.es = np.zeros(n_cells, dtype=np.double)
        self.cs = np.zeros(n_cells, dtype=np.double)

        self.taus = np.zeros(n_cells, dtype=np.double)
        self.D_left = np.zeros(n_cells, dtype=np.double)
        self.D_right = np.zeros(n_cells, dtype=np.double)

        self.fluxes = np.zeros((n_qs, n_cells+1), dtype=np.double)

        self.qs = np.zeros((n_qs, n_cells), dtype=np.double)

        self.hs = np.zeros((n_qs, n_cells), dtype=np.double)
        # self.mega_foo_results = np.zeros((n_cells+1, ), dtype=np.double)

    cpdef GasLayer copy(self):
        cdef GasLayer res = GasLayer(self.n_cells, self.tube, self.gasEOS, self.flux_calculator, self.grid_strecher, self.n_qs)
        res.time = self.time

        res.xs_cells[:] = self.xs_cells
        res.xs_borders[:] = self.xs_borders
        res.Vs_borders[:] = self.Vs_borders

        res.S[:] = self.S 

        res.ds[:] = self.ds
        res.W[:] = self.W

        res.ps[:] = self.ps
        res.ros[:] = self.ros
        res.us[:] = self.us
        res.es[:] = self.es
        res.cs[:] = self.cs

        res.taus[:] = self.taus
        res.D_left[:] = self.D_left
        res.D_right[:] = self.D_right

        res.fluxes[:] = self.fluxes
        res.qs[:] = self.qs
        res.hs[:] = self.hs
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

    cpdef void fill_fluxes_taus(self):
        self.flux_calculator.fill_fluxes_taus(self)

    cpdef GasLayer step_Godunov_simple(self, double v_left, double v_right, double courant, bint init_taus_acustic, double alpha=1):
        if init_taus_acustic:
            self.init_taus_acustic()
        cdef double tau = self.get_tau_min() * courant * alpha  
        self.init_h()
        cpdef GasLayer layer1 = self.copy() 
        cdef bint suc = self.grid_strecher.sync_layers(self, layer1, tau, v_left, v_right)
        if not suc:
            print('suc 554')
            return self
        
        layer1.init_SsdW()
        self.fill_fluxes_taus()
        if self.get_tau_min() < tau/alpha:
            print('tau stuff 560')
            return self.step_Godunov_simple(v_left, v_right, courant, False)
        layer1.taus[:] = self.taus
        cdef size_t i, j
        cdef double dx
        for i in range(layer1.xs_cells.shape[0]):
            dx = self.xs_borders[i+1] - self.xs_borders[i]
            for j in range(layer1.n_qs):
                layer1.qs[j, i] = (self.qs[j, i]*self.W[i] - tau*(self.S[i+1]*self.fluxes[j, i+1] - self.S[i]*self.fluxes[j, i]) + tau*self.hs[j, i]*dx)/layer1.W[i]
            
        layer1.init_ropue()
        return layer1
        
    cpdef GasLayer step_Godunov_corrector2(self, GasLayer layer_simple, double v_left, double v_right):
        cdef double tau = (layer_simple.time - self.time)*2
        cpdef GasLayer layer2 = layer_simple.copy()
        cdef bint suc = self.grid_strecher.sync_layers(layer_simple, layer2, tau/2, v_left, v_right)
        if not suc:
            print('suc 556')
            return self        
        layer2.init_SsdW()
        layer_simple.init_h()
        layer_simple.fill_fluxes_taus()
        cdef size_t i, j
        cdef double dx
        for i in range(layer2.xs_cells.shape[0]):
            dx = layer_simple.xs_borders[i+1] - layer_simple.xs_borders[i]
            for j in range(layer2.n_qs):
                layer2.qs[j, i] = (self.qs[j, i]*self.W[i] \
                    - tau*(layer_simple.S[i+1]*layer_simple.fluxes[j, i+1] \
                         - layer_simple.S[i]*layer_simple.fluxes[j, i]) \
                    + tau*layer_simple.hs[j, i]*dx) / layer2.W[i]
        layer2.init_ropue()
        return layer2

    cpdef GasLayer step_Godunov_corrector(self, GasLayer layer_simple, double v_left, double v_right):
        cdef double tau = layer_simple.time - self.time
        cpdef GasLayer layer2 = self.copy()
        cdef bint suc = self.grid_strecher.sync_layers(self, layer2, tau, v_left, v_right)
        if not suc:
            print('suc 555')
            return self
        layer2.init_SsdW()

        layer_simple.Vs_borders[:] = self.Vs_borders
        layer_simple.fill_fluxes_taus()
        layer_simple.init_h()

        cdef size_t i, j
        cdef double dx, dx_simple
        for i in range(layer2.xs_cells.shape[0]):
            dx = layer_simple.xs_borders[i+1] - layer_simple.xs_borders[i]
            dx_simple = self.xs_borders[i+1] - self.xs_borders[i]
            for j in range(layer2.n_qs):
                layer2.qs[j, i] = ( \
                    self.qs[j, i] * self.W[i] \
                    - tau * ((self.S[i+1] + layer_simple.S[i+1]) * (self.fluxes[j, i+1] + layer_simple.fluxes[j, i+1])*0.25 \
                           - (self.S[i]   + layer_simple.S[i]  ) * (self.fluxes[j, i]   + layer_simple.fluxes[j, i]  )*0.25) \
                    + tau * (self.hs[j, i]*dx +layer_simple.hs[j, i]*dx_simple)*0.5 \
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
            'ds': np.array(self.ds),
            'S': np.array(self.S),
            'W': np.array(self.W),

            'qs':np.array(self.qs),
            'hs':np.array(self.hs),
            'fluxes':np.array(self.fluxes)
       }
    def from_dict(self, d):
        self.time = d['time']
        self.n_cells = d['xs_cells'].shape[0] 
        self.n_qs = d['n_qs']
        self.xs_cells[:] = d['xs_cells']
        self.xs_borders[:] = d['xs_borders']
        self.Vs_borders[:] = d['Vs_borders']

        self.S[:] = d['S'] 

        self.ds[:] = d['ds']
        self.W[:] = d['W']

        self.ps[:] = d['ps']
        self.ros[:] = d['ros']
        self.us[:] = d['us']
        self.es[:] = d['es']
        self.cs[:] = d['cs']

        self.taus[:] = d['taus']
        self.D_left[:] = d['D_left']
        self.D_right[:] = d['D_right']

        self.fluxes[:] = d['fluxes']
        self.qs[:] = d['qs']
        self.hs[:] = d['hs']
    
    

