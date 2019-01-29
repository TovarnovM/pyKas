

cdef class InterpXY:
    cdef public double[:] xs, ys, ks, bs
    cdef int length, n
    cpdef double get_v(self, double x)
    cpdef double[:] get_vs(self, double[:] xs)
    cdef int set_n(self, double t)
    cdef void sync_ks_bs(self)

cdef class Tube:
    cdef InterpXY d, s, w, w_reverse
    cpdef double get_d(self, double x)
    cpdef double[:] get_ds(self, double[:] xs)
    cpdef double[:] get_xs(self)
    cpdef double[:] get_dsdx(self, double[:] xs)    
    cpdef double[:] get_W(self, double[:] xs)
    cpdef double[:] get_S(self, double[:] xs)
    cpdef double get_x2(self, double x1, double w)
    cpdef double get_W_between(self, double x1, double x2)
