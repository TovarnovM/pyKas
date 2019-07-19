# distutils: language=c++
# cython: language_level=3, boundscheck=False, nonecheck=False, cdivision=True, initializedcheck=False

from gaslayer cimport GasLayer

cdef class GasPhase(GasLayer):
    pass    
