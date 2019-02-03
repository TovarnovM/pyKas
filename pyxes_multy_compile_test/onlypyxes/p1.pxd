# distutils: language=c++
# cython: language_level=3

cdef int cfoo(int)
cdef class Bar:
    cdef double[:] a, b
    cdef double sum(self, int a=*)