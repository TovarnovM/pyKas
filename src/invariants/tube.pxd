# distutils: language=c++
# cython: language_level=3
cdef class InterpXY:
    cdef int length, n
    cpdef double get_v(self, double)
    cpdef double[:] get_vs(self, double[:])
    cdef int set_n(self, double)
    cdef void sync_ks_bs(self)

cdef class Tube:
    cpdef double get_d(self, double)
    cpdef double[:] get_ds(self, double[:])
    cpdef double[:] get_xs(self)
    cpdef double[:] get_dsdx(self, double[:])    
    cpdef double[:] get_W(self, double[:])
    cpdef double[:] get_S(self, double[:])
    cpdef double get_x2(self, double , double )
    cpdef double get_W_between(self, double, double)
