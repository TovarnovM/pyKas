# distutils: language=c++
# cython: language_level=3, boundscheck=False, nonecheck=False, cdivision=True, initializedcheck=False

from tube cimport Tube
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

    cpdef void fill_fluxes_Ds(self, GasLayer layer):
        if self.flux_type == 1:
            self.fill_fluxes_Godunov(layer.Vs_borders, layer.ps, layer.ros, layer.us, layer.cs, layer.gasEOS,\
                layer.flux1, layer.flux2, layer.flux3, layer.D_left, layer.D_right)

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

cdef class GridStrecher(object):
    def __init__(self, strech_type=1):
        self.strech_type = strech_type
        self.bufarr = np.zeros(7)

    cpdef void init_regular(self, double v1, double v2, double[:] vs):
        cdef int n = vs.shape[0]
        cdef double dv = (v2-v1)/(n-1)
        cdef size_t i
        for i in range(vs.shape[0]):
            vs[i] = v1 + i*dv

    cpdef void fill_euler_vel0_regular(self, double tau, double[:] xs0, double[:] xs1, double[:] vel0_fill):
        cdef size_t i
        for i in range(vel0_fill.shape[0]):
            vel0_fill[i] = (xs1[i] - xs0[i])/tau

    cpdef void fill_xs_cells(self, double[:] xs_borders, double[:] xs_fill):
        cdef size_t i
        for i in range(xs_fill.shape[0]):
            xs_fill[i] = (xs_borders[i+1] + xs_borders[i])/2

    cpdef void smooth_arr(self, double[:] xs, double[:] vs, double[:] vs_smoothed, double window_part=0.1):
        pass

    


cdef class GasLayer(object):
    def __init__(self, n_cells, Tube tube, GasEOS gasEOS, GasFluxCalculator flux_calculator):
        self.tube = tube
        self.gasEOS = gasEOS
        self.time = 0
        self.flux_calculator = flux_calculator
        
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

        self.flux1 = np.zeros(n_cells+1, dtype=np.double)
        self.flux2 = np.zeros(n_cells+1, dtype=np.double)
        self.flux3 = np.zeros(n_cells+1, dtype=np.double)

        self.q1 = np.zeros(n_cells, dtype=np.double)
        self.q2 = np.zeros(n_cells, dtype=np.double)
        self.q3 = np.zeros(n_cells, dtype=np.double)

    cpdef GasLayer copy(self):
        cdef GasLayer res = GasLayer(self.n_cells, self.tube, self.gasEOS, self.flux_calculator)
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

        res.flux1[:] = self.flux1
        res.flux2[:] = self.flux2
        res.flux3[:] = self.flux3

        res.q1[:] = self.q1
        res.q2[:] = self.q2
        res.q3[:] = self.q3
        return res


    cpdef void init_ropue_fromfoo(self, foo_ropu):
        cdef size_t i
        cdef double x
        cdef GasEOS eos = self.gasEOS
        cdef double ro, p, u, e
        for i in range(self.n_cells):
            ro, p, u = foo_ropu(self.xs[i])
            self.ros[i], self.ps[i], self.us[i] = ro, p, u
            self.es[i] = eos.get_e(ro, p)
            self.cs[i] = eos.get_csound(ro, p)
    
    cpdef void init_SsdW(self):
        self.tube.fill_S(self.xs_borders, self.S)
        self.tube.fill_dsdx(self.xs_borders, self.ds)
        self.tube.fill_W(self.xs_borders, self.W)
    
    cpdef void init_q(self):
        cdef size_t i
        cdef size_t n = self.ros.shape[0]
        cdef double q1, q2, q3
        for i in range(n):
            q1, q2, q3 = roue_to_q(self.ros[i], self.us[i], self.es[i])
            self.q1[i] = q1
            self.q2[i] = q2 
            self.q3[i] = q3

    cpdef void init_ropue(self):
        cdef size_t i
        cdef size_t n = self.ros.shape[0]
        cdef double ro, p, u, e
        for i in range(n):
            ro, u, e = q123_to_roue(self.q1[i], self.q2[i], self.q3[i])
            p = self.gasEOS.get_p(ro, e)
            self.ros[i] = ro
            self.ps[i] = p
            self.us[i] = u
            self.es[i] = e
            self.cs[i] = self.gasEOS.get_csound(ro, p)

    # cpdef init_fluxes_AUSM()
    
    

