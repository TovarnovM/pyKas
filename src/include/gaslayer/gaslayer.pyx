# distutils: language=c++
# cython: language_level=3

from tube.tube cimport Tube

cdef void foo():
    cdef Tube t = Tube([1,2,3],[3,3,3])
    cdef double d = t.get_d(4)
    cdef double dd = d*dd
    