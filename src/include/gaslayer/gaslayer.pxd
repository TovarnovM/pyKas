# distutils: language=c++
# cython: language_level=3, boundscheck=False, nonecheck=False, cdivision=True, initializedcheck=False

from tube cimport Tube, InterpXY

cpdef double foo()

cpdef inline void roue_to_q_(
    double ro,
    double u,
    double e,
    double[:] q) nogil

cpdef inline double abs(double x) nogil

cpdef inline (double, double, double) q_to_roue(double[:] q) nogil

cpdef inline (double, double, double) q123_to_roue(double q1, double q2, double q3) nogil

cpdef inline (double, double, double) roue_to_q(
    double ro,
    double u,
    double e) nogil

cpdef inline double roe_to_p(
    double ro,
    double e,
    double gamma,
    double b) nogil

cpdef inline double rop_to_e(
    double ro,
    double p,
    double gamma,
    double b) nogil

cpdef inline double rop_to_csound(
    double ro,
    double p,
    double gamma,
    double b) nogil

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
    double vbi) nogil

cpdef inline double min2(double a, double b) nogil

cpdef inline double max2(double a, double b) nogil
    
cdef class GasEOS:
    cdef public double gamma, kappa, p_0, c_0
    cdef public int kind 
    cpdef double get_e(self, double ro, double p)
    cpdef double get_p(self, double ro, double e)
    cpdef double get_csound(self, double ro, double p)  

cdef class Powder(GasEOS):
    cdef public InterpXY psi, dpsi
    cdef public str name
    cdef public double z_k, I_k, alpha_k, ro, f, T_1, nu, R_g
    cpdef double get_e_powder(self, double ro, double p, double z)
    cpdef double get_p_powder(self, double ro, double e, double z)
    cpdef double get_csound_powder(self, double ro, double p, double z)
    cpdef double get_z_k(self)

cdef class GasFluxCalculator:
    cdef public int flux_type, left_border_type, right_border_type, n_iter_max
    cdef public double epsF
    cpdef void fill_fluxesURP_Ds_Godunov(self, double[:] Vs_borders, double[:] ps, double[:] ros, \
            double[:] us, double[:] cs, GasEOS gasEOS, \
            double[:] flux1, double[:] flux2, double[:] flux3, double[:] D_left, double[:] D_right, double[:] U_kr)
    cpdef void fill_fluxesURP(self, GasLayer layer)
    

cdef class GridStrecher:
    cdef public int strech_type
    cdef public double st2_window_part, st2_adapt_prop
    cdef public double[:] bufarr,bufarr_border
    cdef public InterpXY interp_smooth, interp_adapt
    cpdef bint evaluate(self, double tau, GasLayer layer)
    cpdef void init_regular(self, double v1, double v2, double[:] vs)
    cpdef void fill_euler_vel0_regular(self, double tau, double[:] xs0, double[:] xs1, double[:] vel0_fill)
    cpdef void fill_xs_cells(self, double[:] xs_borders, double[:] xs_fill)
    cpdef void smooth_arr(self, double[:] xs, double[:] vs, double[:] vs_smoothed, double window_part=*)
    cpdef void adaptine_borders(self, double[:] xs_borders, double[:] vs, double[:] xs_adapt)
    cpdef void fill_Vs_borders_proportional(self, double v_left, double v_right, double[:] xs_borders, double[:] Vs_borders)
    cpdef void strech_layer_adaptive(self, GasLayer layer, GasLayer layer_prev, double tau)
    cpdef bint sync_layers(self, GasLayer layer0, GasLayer layer1, double tau, double v_left, double v_right)

cdef class GasLayer:
    cdef public double[:] xs_cells, xs_borders, Vs_borders, U_kr,\
        ros, ps, us, es, cs, taus, D_left, D_right,  \
        ds, S, W
    cdef public double[:,:] qs, hs, fluxes
        #pars # 0 - ros, 1 - ps, 2 - us  
    cdef public int n_cells, n_qs
    cdef public Tube tube
    cdef public GasEOS gasEOS
    cdef public double time
    cdef public GasFluxCalculator flux_calculator
    cdef public GridStrecher grid_strecher
    cpdef void copy_params_to(self, GasLayer to_me)
    cpdef GasLayer copy(self)
    cpdef void init_ropue_fromfoo(self, foo_ropu, bint init_q=*, bint init_SsdW=*)
    cpdef void init_SsdW(self)
    cpdef void init_q(self)
    cpdef void init_ropue(self)
    cpdef void init_h(self)
    cpdef double get_tau_min(self)
    cpdef void init_taus_acustic(self)
    cpdef GasLayer step_simple(self, double tau, double v_left, double v_right)
    cpdef GasLayer step_Godunov_simple(self, double v_left, double v_right, double courant, bint init_taus_acustic, double alpha=*)
    cpdef void fill_fluxes(self)
    cpdef void fill_taus(self)
    cpdef GasLayer step_Godunov_corrector(self, GasLayer layer_simple, double v_left, double v_right)
    cpdef GasLayer step_Godunov_corrector2(self, GasLayer layer_simple, double v_left, double v_right)
    cpdef void fill_fluxesURP(self)

cdef class PowderOvLayer(GasLayer):
    cdef public double some
    cdef public double[:] zs
    cpdef void copy_params_to_Ov(self, PowderOvLayer to_me)
    cpdef GasLayer copy(self)
    cpdef void init_ropue_fromfoo(self, foo_ropu, bint init_q=*,  bint init_SsdW=*)
    cpdef void init_q(self)
    cpdef void init_ropue(self)
    cpdef void init_h(self)
    cpdef void fill_fluxes(self)

cdef class ElPistLayer(GasLayer):
    cdef public double[:] tauxx_flux
    cpdef void copy_params_to_elpist(self, ElPistLayer lr)
    cpdef GasLayer copy(self)
    cpdef void init_h(self)

cdef class ElPistEOS(GasEOS):
    cdef public double ro_0, sigma_star, k_0,b_1,b_2,tau_0,mu,tau_s
    cpdef double get_tauu(self, double sigma, double u)
    cpdef double get_kh(self, double h)