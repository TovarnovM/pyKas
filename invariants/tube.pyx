import numpy as np
from copy import deepcopy
import cython
from libc.math cimport pi

cdef class InterpXY(object):
    cdef double[:] xs, ys, ks, bs
    cdef int length, n
    def __init__(self, xs, ys):
        self.xs = np.array(xs, dtype=np.double)
        self.ys = np.array(ys, dtype=np.double)
        if len(xs) != len(ys):
            raise AttributeError(f'Длины xs ({len(xs)}) и ys ({len(ys)}) не одинаковы')
        cdef double[:] kv = np.zeros(len(xs), dtype=np.double)
        self.ks = kv
        cdef double[:] bv = np.zeros(len(xs), dtype=np.double)
        self.bs = bv
        self.sync_ks_bs()
        self.n = 0
    
    @cython.boundscheck(False)
    @cython.wraparound(False)   
    cpdef double get_v(self, double x):
        cdef double tt = x
        cdef int n = self.set_n(tt)
        if n < 0 or n == (self.length-1) or self.length == 1:
            n = 0 if n<0 else n
            self.n = n
            return self.ys[n]
        self.n = n
        return self.ks[n] * tt + self.bs[n]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef double[:] get_vs(self, double[:] xs):
        cdef double[:] res = np.empty(xs.shape, dtype=np.double)
        cdef int i
        for i in range(xs.shape[0]):
            res[i] = self.get_v(xs[i])
        return res
        
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int set_n(self, double t):
        cdef int n = self.n
        if n < 0:
            if self.xs[0] > t:
                return -1
            n = 0
        cdef int minw, maxw;
        cdef int lengthM1 = self.length - 1;
        if self.xs[n] <= t:
            if n == lengthM1 or self.xs[n + 1] > t:
                return n
            n += 1
            if n == lengthM1 or self.xs[n + 1] > t:
                return n
            if self.xs[self.length-1] <= t:
                return lengthM1
            minw = n
            maxw = lengthM1
        else:
            if n == 0 or self.xs[n - 1] <= t:
                return n-1
            if self.xs[0] > t:
                n = -1
                return n
            minw = 0
            maxw = n
        while minw != maxw:
            n = (minw + maxw) // 2
            if self.xs[n] <= t:
                if self.xs[n + 1] > t:
                    return n
                minw = n
            else:
                maxw = n
        n = minw
        return n
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void sync_ks_bs(self):
        self.length = len(self.xs)
        cdef size_t i
        for i in range(self.length):
            self.ks[i] = (self.ys[i + 1] - self.ys[i]) / (self.xs[i + 1] - self.xs[i])
            self.bs[i] = self.ys[i] - self.ks[i] * self.xs[i]

cdef class Tube(object):
    cdef InterpXY d, s, w, w_reverse
    def __init__(self, xs, ds):
        xs = np.asarray(xs, dtype=np.double)
        dd = np.asarray(ds, dtype=np.double)
        self.s = InterpXY(xs, dd)
        ss = dd ** 2 * pi * 0.25
        self.s = InterpXY(xs, ss)
        def ConeW(d1, d2, h):
            return pi * h * (d1 * d1 + d1 * d2 + d2 * d2) / 12

        xss = [xs[0] - 0.1] + [x for x in xs] + [xs[-1] + 0.1]
        dss = [ds[0]] + [d for d in ds] + [ds[-1]]

        ws = [0]
        for x1, x2, d1, d2 in zip(xss, xss[1:], dss, dss[1:]):
            ws.append(ConeW(d1, d2, x2 - x1) + ws[-1])

        self.w = InterpXY(xss, ws)
        self.w_reverse = InterpXY(ws, xss)
    
    @cython.boundscheck(False)
    @cython.wraparound(False)    
    cpdef double[:] get_dsdx(self, double[:] xs):
        cdef double[:] res = np.empty(xs.shape[0]-1, dtype=np.double)
        cdef int i
        cdef double dx, si ,si1
        si = self.s.get_v(xs[0])
        for i in range(res.shape[0]):
            si1  = self.s.get_v(xs[i+1])
            dx = xs[i+1] - xs[i]
            res[i] = (si1-si)/dx
            si = si1
        return res
    
    @cython.boundscheck(False)
    @cython.wraparound(False)      
    cpdef double[:] get_W(self, double[:] xs):
        cdef double[:] res = np.empty(xs.shape[0]-1, dtype=np.double)
        cdef int i
        cdef double wi ,wi1
        wi = self.w.get_v(xs[0])
        for i in range(res.shape[0]):
            wi1  = self.w.get_v(xs[i+1])
            res[i] = wi1 - wi
            wi = wi1
        return res  
        
    cpdef double[:] get_S(self, double[:] xs):
        return self.s.get_vs(xs) 
    
    cpdef double get_x2(self, double x1, double w):
        cdef double w1 = self.w.get_v(x1)
        return self.w_reverse.get_v(w1 + w)
    
    cpdef double get_W_between(self, double x1, double x2):
        return self.w.get_v(x2) - self.w.get_v(x1)
