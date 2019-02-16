# distutils: language=c++
# cython: language_level=3

from tube cimport Tube

cpdef double foo():
    cdef Tube t = Tube([1,2,3],[3,3,3])
    cdef double d = t.get_d(4)
    return d

    
