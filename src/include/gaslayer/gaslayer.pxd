# distutils: language=c++
# cython: language_level=3

cpdef double foo()

cpdef inline void roue_to_q_(
    double ro,
    double u,
    double e,
    double[:] q)

cpdef inline (double, double, double) q_to_roue(double[:] q)

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

cdef class GasLayer:
    cdef public double[:] xs_cells, xs_borders, ps, ros, us, es