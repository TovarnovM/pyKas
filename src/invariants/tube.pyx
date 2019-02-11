# distutils: language=c++
# cython: language_level=3

import numpy as np
from copy import deepcopy
import cython
from libc.math cimport pi

cdef class InterpXY(object):
    """
    __init__(self, xs, ys)
    Класс для БЫСТРОЙ линейной интерполяции точечной одномерной функции
    xs - array like - занчения абсцисс функции
    ys - array like - значения ординат функции
    """

    def __init__(self, xs, ys):
        """
        __init__(self, xs, ys)
        Класс для БЫСТРОЙ линейной интерполяции точечной одномерной функции
        xs - array like - занчения абсцисс функции
        ys - array like - значения ординат функции        
        """
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

    def __call__(self, x):
        """
        Возвращает интерполированное(-ые) значение(-ия) функции в точке(-ах) 'x'
        x - list/ndarray/float/int...
        """
        if isinstance(x, list):
            return np.asarray(self.get_vs(np.array(x, dtype=np.double)))
        if isinstance(x, np.ndarray):
            return np.asarray(self.get_vs(x))
        return self.get_v(x)
    
    def __repr__(self):
        xs = list(np.array(self.xs))
        ys = list(np.array(self.ys))
        return f'InterpXY(xs={xs}, ys={ys})'
    
    def __str__(self):
        return repr(self)
    
    @cython.boundscheck(False)
    @cython.wraparound(False)   
    cpdef double get_v(self, double x):
        """
        double get_v(self, double x):
        Возвращает интерполированное значение функции в точке 'x'
        """
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
        """
        double[:] get_vs(self, double[:] xs)
        Возвращает интерполированные значения функции в точках 'x'
        """
        cdef double[:] res = np.empty(len(xs), dtype=np.double)
        cdef size_t i
        for i in range(len(res)):
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
        cdef int minw, maxw
        cdef int lengthM1 = self.length - 1
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
    """
    Класс 'труба'
    def __init__(self, xs, ds, deltax4ws=1000.0)
        xs - array like - координаты точек по оси
        ds - array like - диаметры трубы в этич точках
        deltax4ws - double - на сколько метров от крайних точек можно так же считать объемы
    """
    def __init__(self, xs, ds, deltax4ws=1000.0):
        """
        xs - array like - координаты точек по оси
        ds - array like - диаметры трубы в этич точках
        deltax4ws - double - на сколько метров от крайних точек можно так же считать объемы
        """
        xs = np.array(xs, dtype=np.double)
        dd = np.array(ds, dtype=np.double)
        self.d = InterpXY(xs, dd)
        ss = dd ** 2 * pi * 0.25
        self.s = InterpXY(xs, ss)
        def ConeW(d1, d2, h):
            return pi * h * (d1 * d1 + d1 * d2 + d2 * d2) / 12

        xss = [xs[0] - 1000.0] + [x for x in xs] + [xs[-1] + 1000.0]
        dss = [ds[0]] + [d for d in ds] + [ds[-1]]

        ws = [0]
        for x1, x2, d1, d2 in zip(xss, xss[1:], dss, dss[1:]):
            ws.append(ConeW(d1, d2, x2 - x1) + ws[-1])

        self.w = InterpXY(xss, ws)
        self.w_reverse = InterpXY(ws, xss)
        
    def __repr__(self):
        xs = list(np.array(self.get_xs()))
        ys = list(np.array(self.get_ds()))
        return f'Tube(xs={xs}, ds={ys})'
    
    def __str__(self):
        return repr(self)
    
    cpdef double get_d(self, double x):
        """
        double get_d(self, double x)
        Возвращает диаметр в точке x [double]
        """
        return self.d.get_v(x)
    
    cpdef double[:] get_ds(self, double[:] xs=None):
        """
        double[:] get_ds(self, double[:] xs=None)
        Возвращает диаметры в точках xs
        Если xs is None -> возвращает диаметры точек (из конструктора), образующих трубу 
        """
        if xs is None:
            return np.array(self.d.ys)
        return self.d.get_vs(xs)
    
    cpdef double[:] get_xs(self):
        """
        double[:] get_xs(self)
        Возвращает координаты 'xs' точек из конструктора
        """
        return np.array(self.d.xs)
    
    @cython.boundscheck(False)
    @cython.wraparound(False)    
    cpdef double[:] get_dsdx(self, double[:] xs):
        """
        double[:] get_dsdx(self, double[:] xs)
        возвращает производные площадей трубы в точках xs, 
        НО это производные (Sx_i+1 - Sx_i)/dx, а не самой трубы. 
        Так же len(result) == len(xs) - 1
        """
        cdef double[:] res = np.empty(len(xs)-1, dtype=np.double)
        cdef size_t i
        cdef double dx, si ,si1
        si = self.s.get_v(xs[0])
        for i in range(len(res)):
            si1  = self.s.get_v(xs[i+1])
            dx = xs[i+1] - xs[i]
            res[i] = (si1-si)/dx
            si = si1
        return res
    
    @cython.boundscheck(False)
    @cython.wraparound(False)      
    cpdef double[:] get_W(self, double[:] xs):
        """
        double[:] get_W(self, double[:] xs)
        метод возвращает объемы трубы, находящиеся между точками
        len(result) == len(xs) - 1
        """
        cdef double[:] res = np.empty(len(xs)-1, dtype=np.double)
        cdef size_t i
        cdef double wi ,wi1
        wi = self.w.get_v(xs[0])
        for i in range(len(res)):
            wi1  = self.w.get_v(xs[i+1])
            res[i] = wi1 - wi
            wi = wi1
        return res  
        
    cpdef double[:] get_S(self, double[:] xs):
        """
        double[:] get_S(self, double[:] xs)
        Возвращает площади сечений трубы в точках xs
        """
        return self.s.get_vs(xs) 
    
    cpdef double get_x2(self, double x1, double w):
        """
        double get_x2(self, double x1, double w)
        возвращает координату x2, такую, что объем трубы между x1 и x2 равен w
        """
        cdef double w1 = self.w.get_v(x1)
        return self.w_reverse.get_v(w1 + w)
    
    cpdef double get_W_between(self, double x1, double x2):
        """
        double get_W_between(self, double x1, double x2)
        Возвращает объем трубы между точками x1 и x2
        """
        return self.w.get_v(x2) - self.w.get_v(x1)
