# distutils: language=c++
# cython: language_level=3, boundscheck=False, nonecheck=False, cdivision=True, initializedcheck=False

from gaslayer cimport GasLayer

cdef class GasPhase(GasLayer):
    cdef public double[:] alphas
    cpdef void copy_params_to_gas_phase(self, GasPhase to_me)
    cpdef GasLayer copy(self)

cdef class SolidPhase(GasLayer):
    pass    
