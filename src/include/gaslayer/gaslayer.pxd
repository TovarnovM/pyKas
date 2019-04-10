# distutils: language=c++
# cython: language_level=3, boundscheck=False, nonecheck=False, cdivision=True, initializedcheck=False

from tube cimport Tube

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
    cdef public int flux_type, left_border_type, right_border_type
    cpdef void fill_fluxes_Godunov(self, double[:] Vs_borders, double[:] ps, double[:] ros, \
            double[:] us, double[:] cs, GasEOS gasEOS, \
            double[:] flux1, double[:] flux2, double[:] flux3)
    cpdef void fill_fluxes(self, GasLayer layer)


cdef class GridStrecher:
    cdef public int strech_type
    cdef public double[:] bufarr
    cpdef void init_regular(self, double v1, double v2, double[:] vs)
    cpdef void fill_euler_vel0_regular(self, double tau, double[:] xs0, double[:] xs1, double[:] vel0_fill)

cdef class GasLayer:
    cdef public double[:] xs_cells, xs_borders, Vs_borders, \
        ps, ros, us, es, cs, \
        flux1, flux2, flux3, \
        q1, q2, q3, \
        ds, S, W

    cdef public Tube tube
    cdef public GasEOS gasEOS
    cdef public double time
    cdef public GasFluxCalculator flux_calculator
    cpdef GasLayer copy(self)
    cpdef void init_ropue_fromfoo(self, foo_ropu)
    cpdef void init_SsdW(self)
    cpdef void init_q(self)
    cpdef void init_ropue(self)
