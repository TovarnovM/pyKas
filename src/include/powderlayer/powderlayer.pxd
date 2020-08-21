# distutils: language=c++
# cython: language_level=3, boundscheck=False, nonecheck=False, cdivision=True, initializedcheck=False
from tube cimport Tube, InterpXY
from gaslayer cimport GasLayer, GasEOS, GasFluxCalculator, GridStrecher

cpdef void roue_to_q_(double ro, double u, double e, double[:] q) nogil
cpdef (double, double, double) roue_to_q(double ro, double u, double e) nogil
cpdef (double, double, double) q123_to_roue(double q1, double q2, double q3) nogil

cdef class Powder(GasEOS):
    cdef public InterpXY psi, dpsi
    cdef public str name
    cdef public double z_k, I_k, alpha_k, ro, f, T_1, nu, R_g
    cpdef double get_e_powder(self, double ro, double p, double z)
    cpdef double get_p_powder(self, double ro, double e, double z)
    cpdef double get_csound_powder(self, double ro, double p, double z)
    cpdef double get_z_k(self)

cdef class PowderOvLayer(GasLayer):
    cdef public double t_ign
    cdef public double[:] zs
    cpdef void copy_params_to_Ov(self, PowderOvLayer to_me)
    cpdef GasLayer copy(self)
    cpdef void init_ropue_fromfoo(self, foo_ropu, int init_q=*,  int init_SsdW=*)
    cpdef void init_q(self)
    cpdef void init_ropue(self)
    cpdef void init_h(self)
    cpdef void fill_fluxes(self)