# distutils: language=c++
# cython: language_level=3, boundscheck=False, nonecheck=False, cdivision=True, initializedcheck=False

from tube cimport Tube, InterpXY

cpdef double foo()

cpdef  void roue_to_q_(
    double ro,
    double u,
    double e,
    double[:] q) nogil

cpdef  double abs(double x) nogil

cpdef  (double, double, double) q_to_roue(double[:] q) nogil

cpdef  (double, double, double) q123_to_roue(double q1, double q2, double q3) nogil

cpdef  (double, double, double) roue_to_q(
    double ro,
    double u,
    double e) nogil

cpdef  double roe_to_p(
    double ro,
    double e,
    double gamma,
    double b) nogil

cpdef  double rop_to_e(
    double ro,
    double p,
    double gamma,
    double b) nogil

cpdef  double rop_to_csound(
    double ro,
    double p,
    double gamma,
    double b) nogil

cpdef void AUSM_gas(
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
    double vbi,
    double[:] rr_val) nogil

cpdef  double min2(double a, double b) nogil

cpdef  double max2(double a, double b) nogil

cpdef double get_extrapol_2e_v15(double v2, double v3) nogil
cpdef double get_extraopl_3e_v15(double v2, double v3, double v4) nogil
cpdef double get_interp_3i_v15(double v1, double v2, double v3) nogil
cpdef double get_final_limited_value(double v_init, double v1, double v2, double delta_L=*) nogil

cdef class GasEOS:
    cdef public double gamma, kappa, p_0, c_0
    cdef public int kind 
    cpdef double get_e(self, double ro, double p)
    cpdef double get_p(self, double ro, double e)
    cpdef double get_csound(self, double ro, double p)  


cdef class GasFluxCalculator:
    cdef public int x_order, flux_type, left_border_type, right_border_type, n_iter_max, rr_vals_len, rr_bint_len
    cdef public double epsF, delta_L, alpha_1, alpha_2
    cpdef void set_flux_type(self, int new_flux_type)
    cpdef void create_rr_arrs(self, GasLayer layer)
    cpdef void fill_rr_vals(self, GasLayer layer)
    cpdef void fill_rr_vals_borders_fill_fluxesURP(self, GasLayer layer)
    cpdef void fill_rr_vals_Riman(self, GasLayer layer)
    cpdef void fill_fluxesURP_Godunov(self, GasLayer layer)
    cpdef void fill_fluxesURP(self, GasLayer layer)
    cpdef double get_interextra_value(self, double v1, double v2, double v3, double v4)
    cpdef double get_v_right_forRiman(self, double[:] vs, double[:] bettas, int i)
    cpdef double get_v_left_forRiman(self, double[:] vs, double[:] bettas, int i)
    cpdef double get_bogdanov_v(self, double v1, double v2, double v3, double v4, double betta_max)
    cpdef void fill_fluxesURP_AUSM(self, GasLayer layer)
    

cdef class GridStrecher:
    cdef public int strech_type, index_anchor
    cdef public double st2_window_part, st2_adapt_prop, D_mnj
    cdef public double[:] bufarr,bufarr_border
    cdef public InterpXY interp_smooth, interp_adapt
    cpdef bint evaluate(self, double tau, GasLayer layer)
    cpdef void init_regular(self, double v1, double v2, double[:] vs)
    cpdef void fill_euler_vel0_regular(self, double tau, double[:] xs0, double[:] xs1, double[:] vel0_fill)
    cpdef void fill_xs_cells(self, double[:] xs_borders, double[:] xs_fill)
    cpdef void smooth_arr(self, double[:] xs, double[:] vs, double[:] vs_smoothed, double window_part=*)
    cpdef void adaptine_borders(self, double[:] xs_borders, double[:] vs, double[:] xs_adapt)
    cpdef void fill_Vs_borders_proportional(self, double v_left, double v_right, double[:] xs_borders, double[:] Vs_borders)
    cpdef void fill_Vs_borders_anchor(self, double v_left, double v_right, double[:] xs_borders, double[:] Vs_borders)
    cpdef void strech_layer_adaptive(self, GasLayer layer, GasLayer layer_prev, double tau)
    cpdef bint sync_layers(self, GasLayer layer0, GasLayer layer1, double tau, double v_left, double v_right)
    cpdef void sync_xs_borders(self, double[:] xs_borders, double xl, double xr)

cdef class GasLayer:
    cdef public double[:] xs_cells, xs_borders, Vs_borders, U_kr,\
        ros, ps, us, es, cs, taus, D_left, D_right,  \
        ds, S, W, bettas
    cdef public double[:,:] qs, hs, fluxes, rr_vals
        #pars # 0 - ros, 1 - ps, 2 - us  
    cdef public bint[:,:] rr_bint
    cdef public int n_cells, n_qs
    cdef public Tube tube
    cdef public GasEOS gasEOS
    cdef public double time, xl, xr
    cdef public GasFluxCalculator flux_calculator
    cdef public GridStrecher grid_strecher
    cdef public str color_4_plot
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
    cpdef void fill_fluxes(self)
    cpdef void fill_taus(self)
    cpdef void fill_rr(self)
    cpdef void fill_fluxesURP(self)
    cpdef GasLayer corrector(self, GasLayer lr_simple05, double tau, double v_left, double v_right)
    cpdef GasLayer corrector_bogdanov(self, GasLayer lr_simple, double v_left, double v_right)
    cpdef double get_p_left(self)
    cpdef double get_p_right(self)
    cpdef double get_E_potential(self)
    cpdef double get_E_kinetic(self)
    cpdef double get_E_sum(self)
