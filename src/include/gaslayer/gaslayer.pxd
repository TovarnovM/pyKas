# distutils: language=c++
# cython: language_level=3, boundscheck=False, nonecheck=False, cdivision=True, initializedcheck=False

from tube cimport Tube, InterpXY

cpdef double foo()

cpdef inline void roue_to_q_(
    double ro,
    double u,
    double e,
    double[:] q)

cpdef inline (double, double, double) q_to_roue(double[:] q)

cpdef inline (double, double, double) q123_to_roue(double q1, double q2, double q3)

cpdef inline (double, double, double) roue_to_q(
    double ro,
    double u,
    double e)

cpdef inline double roe_to_p(
    double ro,
    double e,
    double gamma,
    double b)

cpdef inline double rop_to_e(
    double ro,
    double p,
    double gamma,
    double b)

cpdef inline double rop_to_csound(
    double ro,
    double p,
    double gamma,
    double b)

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
    double vbi)

    
cdef class GasEOS:
    cdef public double gamma, kappa, p_0, c_0
    cdef public int kind 
    cpdef double get_e(self, double ro, double p)
    cpdef double get_p(self, double ro, double e)
    cpdef double get_csound(self, double ro, double p)  

cdef class GasFluxCalculator:
    cdef public int flux_type, left_border_type, right_border_type, n_iter_max
    cdef public double epsF
    cpdef void fill_fluxes_Ds_Godunov(self, double[:] Vs_borders, double[:] ps, double[:] ros, \
            double[:] us, double[:] cs, GasEOS gasEOS, \
            double[:] flux1, double[:] flux2, double[:] flux3, double[:] D_left, double[:] D_right)
    cpdef void fill_fluxes_taus(self, GasLayer layer)
    cpdef void fill_taus_Godunov(self, double[:] xs_borders, double[:] Vs_borders, double[:] D_left, double[:] D_right, double[:] taus)

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
    cpdef void strech_layer_adaptive(self, GasLayer layer)
    cpdef bint sync_layers(self, GasLayer layer0, GasLayer layer1, double tau, double v_left, double v_right)

cdef class GasLayer:
    cdef public double[:] xs_cells, xs_borders, Vs_borders, \
        ps, ros, us, es, cs, taus, D_left, D_right, \
        ds, S, W
    cdef public double[:,:] qs, hs, fluxes
    cdef public int n_cells, n_qs
    cdef public Tube tube
    cdef public GasEOS gasEOS
    cdef public double time
    cdef public GasFluxCalculator flux_calculator
    cdef public GridStrecher grid_strecher
    cpdef GasLayer copy(self)
    cpdef void init_ropue_fromfoo(self, foo_ropu, bint init_q=*, bint init_SsdW=*)
    cpdef void init_SsdW(self)
    cpdef void init_q(self)
    cpdef void init_ropue(self)
    cpdef void init_h(self)
    cpdef double get_tau_min(self)
    cpdef void init_taus_acustic(self)
    cpdef GasLayer step_Godunov_simple(self, double v_left, double v_right, double courant, bint init_taus_acustic)
    cpdef void fill_fluxes_taus(self)