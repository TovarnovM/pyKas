# distutils: language=c++
# cython: language_level=3, boundscheck=False, nonecheck=False, cdivision=True, initializedcheck=False
from tube cimport Tube, InterpXY
from gaslayer cimport GasLayer, GasEOS, GasFluxCalculator, GridStrecher


cdef class ElPistLayer(GasLayer):
    cdef public double[:] tauxx_flux
    cpdef void copy_params_to_elpist(self, ElPistLayer lr)
    cpdef GasLayer copy(self)
    cpdef void init_h(self)

cdef class ElPistEOS(GasEOS):
    cdef public double ro_0, sigma_star, k_0,b_1,b_2,tau_0,mu,tau_s
    cdef public bint zeroP, zeroE
    cpdef double get_tauu(self, double sigma, double u)
    cpdef double get_kh(self, double h)
    cpdef double get_e(self, double ro, double p)
    cpdef double get_p(self, double ro, double e)