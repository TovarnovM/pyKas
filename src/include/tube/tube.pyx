# distutils: language=c++
# cython: language_level=3, boundscheck=False, nonecheck=False, cdivision=True, initializedcheck=False
# cimport numpy as cnp
from copy import deepcopy
import cython
import numpy as np

from libc.math cimport pi, sqrt

cimport numpy as np

cdef class InterpXY(object):
    """
    __init__(self, xs, ys)
    Класс для БЫСТРОЙ линейной интерполяции точечной одномерной функции
    xs - array like - занчения абсцисс функции
    ys - array like - значения ординат функции
    """

    def __init__(self, xs, ys, union_tol=1e-9):
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
        self.union_tol = union_tol

    cdef void set_length(self, int length):
        if length == self.length:
            return
        self.length = length
        self.n = 0
        self.xs = np.zeros(length, dtype=np.double)
        self.ys = np.zeros(length, dtype=np.double)
        self.ks = np.zeros(length, dtype=np.double)
        self.bs = np.zeros(length, dtype=np.double)

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

    cpdef void fill_vs(self, double[:] xs, double[:] vs):
        cdef size_t i
        for i in range(xs.shape[0]):
            vs[i] = self.get_v(xs[i])

    cpdef double[:] get_vs(self, double[:] xs):
        """
        double[:] get_vs(self, double[:] xs)
        Возвращает интерполированные значения функции в точках 'x'
        """
        cdef double[:] result = np.empty(len(xs), dtype=np.double)
        self.fill_vs(xs, result)
        return result
        
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
    

    cdef void sync_ks_bs(self):
        self.length = self.xs.shape[0]
        cdef size_t i
        for i in range(self.length-1):
            self.ks[i] = (self.ys[i + 1] - self.ys[i]) / (self.xs[i + 1] - self.xs[i])
            self.bs[i] = self.ys[i] - self.ks[i] * self.xs[i]

    def get_xs_zeros(self):
        """Возвращает список координат абсцисс, где функция обращается в 0
        None, если функция не имеет нулей
        """
        res = []
        for i in range(self.ys.shape[0]-1):
            if self.ys[i] * self.ys[i+1] <= 0:
                res.append(-self.bs[i]/self.ks[i])
        if len(res) == 0:
            return None
        res2 = [res[0]]
        for el in res:
            if abs(el - res2[-1]) > self.union_tol:
                res2.append(el)
        return res2

    def __getitem__(self, key):
        x = self.xs[key]
        y = self.ys[key]
        return (x, y)

    def __setitem__(self, key, value):
        x, y = value
        self.xs[key] = x
        self.ys[key] = y
        self.sync_ks_bs()

    def __iter__(self):
        return zip(self.xs, self.ys)

    def plot(self, fig, ax, **kwargs):
        """Отрисовать интерполятор

            fig, ax = plt.subplots()

        """
        ax.plot(self.xs, self.ys, **kwargs)

    def __reversed__(self):
        return zip(self.xs[::-1], self.ys[::-1])

    cpdef double[:] union_arrs(self, double[:] a1, double[:] a2):
        u = np.sort(np.union1d(a1, a2))
        res = [u[0]]
        for el in u:
            if abs(el - res[-1]) > self.union_tol:
                res.append(el)
        return np.asarray(res)

    def __add__(self, other):
        if isinstance(other, InterpXY):
            xs = self.union_arrs(self.xs, other.xs)
            ys = [self(x) + other(x) for x in xs]
            return InterpXY(xs, ys, self.union_tol)
        elif isinstance(other, list) or isinstance(other, tuple):
            xs = np.concatenate((self.xs, [other[0]]))
            ys = np.concatenate((self.ys, [other[1]]))
            arr1inds = xs.argsort()
            xs = xs[arr1inds]
            ys = ys[arr1inds] 
            return InterpXY(xs, ys, self.union_tol)
        else:
            ys = np.asarray(self.ys) + other
            return InterpXY(self.xs, ys, self.union_tol)

    def __radd__(self, other):
        return self+other
    
    def __abs__(self):
        zeros = self.get_xs_zeros()
        if zeros:
            xs = np.concatenate((self.xs, zeros))
            ys = np.concatenate((self.ys, np.zeros_like(zeros)))
            arr1inds = xs.argsort()
            xs = xs[arr1inds]
            ys = ys[arr1inds] 
            return InterpXY(xs, np.abs(ys), self.union_tol)
        return InterpXY(self.xs, np.abs(self.ys), self.union_tol)

    def __neg__(self):
        return InterpXY(self.xs, -np.asarray(self.ys), self.union_tol)

    def __mul__(self, other):
        if isinstance(other, InterpXY):
            xs = self.union_arrs(self.xs, other.xs)
            ys = [self(x) * other(x) for x in xs]
            return InterpXY(xs, ys, self.union_tol)
        else:
            ys = np.asarray(self.ys) * other
            return InterpXY(self.xs, ys, self.union_tol)
    
    def __rmul__(self, other):
        return self * other

    def __sub__(self, other):
        if isinstance(other, InterpXY):
            xs = self.union_arrs(self.xs, other.xs)
            ys = [self(x) - other(x) for x in xs]
            return InterpXY(xs, ys, self.union_tol)
        else:
            ys = np.asarray(self.ys) - other
            return InterpXY(self.xs, ys, self.union_tol)
    
    def __rsub__(self, other):
        return -(self - other)

    def __truediv__(self, other):
        if isinstance(other, InterpXY):
            xs = self.union_arrs(self.xs, other.xs)
            ys = [self(x) / other(x) for x in xs]
            return InterpXY(xs, ys, self.union_tol)
        else:
            ys = np.asarray(self.ys) / other
            return InterpXY(self.xs, ys, self.union_tol)
    
    def __rtruediv__(self, other):
        i = InterpXY([self.xs[0], self.xs[1]], [other, other], self.union_tol)
        return i/self

    def integrate(self, x1=None, x2=None):
        if x1 is None:
            x1 = self.xs[0]
        if x2 is None:
            x2 = self.xs[-1]
        cdef double res = 0
        cdef double y1 = self.get_v(x1)
        cdef double xx1 = x1
        cdef double y2, xx2
        cdef int i = self.set_n(x1) + 1
        while i < self.xs.shape[0] and self.xs[i] < x2:
            y2 = self.ys[i]
            xx2 = self.xs[i]
            res += 0.5*(y1+y2)*(xx2-xx1)
            i += 1
            xx1 = xx2
            y1 = y2
        y2 = self.get_v(x2)
        xx2 = x2
        res += 0.5*(y1+y2)*(xx2-xx1)
        return res

    def aver(self):
        if self.xs.shape[0] < 2:
            return self.ys[0]
        integr = self.integrate()
        return integr/(self.xs[-1]-self.xs[0])
    
    @classmethod
    def get_diff_report(cls, int1, int2):
        """Получить отчет о различиях двух InterpXY
        
        Arguments:
            int1 {InterpXY} -- относительно этого будет высчитываться относительные величины
            int2 {InterpXY} -- [description]
        
        Returns:
            dict() -- str --> float

                res['diff_maxs_abs']    -- абсолютное отличие максимальных значений
                res['diff_maxs_otn']    -- отностильеное отличие максимальных значений
                res['diff_maxs_x_abs']  -- абсолютное отличие абсцисс максимальных значений
                res['diff_maxs_x_otn']  -- относительное отличие абсцисс максимальных значений
                res['max_diff']         -- максимальное отличие между функциями
                res['max_diff_x']       -- координата максимального отклонения
                res['diff_integrate']   -- интеграл модуля разлийий между функциями
                res['diff_integrate_otn'] -- приведенный интеграл (относительно интеграла int1) модуля различий между функциями
                res['diff_xs_abs']      -- разность между отрезками, на которых определены функции
                res['diff_xs_otn']      -- относительная разность ...
                res['diff_aver_abs']    -- разность медлу среднеинтегральными значениями
                res['diff_aver_otn']    -- относительная разность ...
        """
        res = {}
        res['diff_maxs_abs'] = np.max(int1.ys) - np.max(int2.ys)
        res['diff_maxs_otn'] = res['diff_maxs_abs'] / np.max(int1.ys)
        res['diff_maxs_x_abs'] = int1.xs[np.argmax(int1.ys)] - int2.xs[np.argmax(int2.ys)]
        res['diff_maxs_x_otn'] = (res['diff_maxs_x_abs']-int1.xs[0]) / (int1.xs[-1]-int1.xs[0])
        diff = abs(int1 - int2)
        i_max = np.argmax(diff.ys)
        res['max_diff'] = diff.ys[i_max]
        res['max_diff_x'] = diff.xs[i_max]
        res['diff_integrate'] = diff.integrate()
        res['diff_integrate_otn'] = res['diff_integrate']/int1.integrate()
        res['diff_xs_abs'] = (int1.xs[-1]-int1.xs[0]) - (int2.xs[-1]-int2.xs[0])
        res['diff_xs_otn'] = res['diff_xs_abs'] / (int1.xs[-1]-int1.xs[0])
        res['diff_aver_abs'] = int1.aver() - int2.aver()
        res['diff_aver_otn'] = res['diff_aver_abs'] / int1.aver()
        return res

    def plot_diff(self, int2, fig, ax, **kwargs):
        """Отрисовать различаи между функциями

            fig, ax = plt.subplots()
            kwargs = {'int1_kwargs': {...},
                      'int2_kwargs': {...},
                      'polygon_kwargs': {...}}
        """
        from matplotlib.patches import Polygon
        ix = self.union_arrs(self.xs, int2.xs)
        verts = [(x, self(x)) for x in ix]
        verts+= [(x, int2(x)) for x in ix[::-1]]
        self.plot(fig, ax, **dict({'lw':2}, **kwargs.get('int1_kwargs', {})))
        int2.plot(fig, ax, **dict({'lw':2}, **kwargs.get('int2_kwargs', {})))
        poly = Polygon(verts, **dict({'alpha':0.3}, **kwargs.get('polygon_kwargs', {})))

        ax.add_patch(poly) 

        
        

    


    

cdef class Tube(object):
    """
    Класс 'труба'
    def __init__(self, xs, ds, deltax4ws=1000.0)
        xs - array like - координаты точек по оси
        ds - array like - диаметры трубы в этич точках
        deltax4ws - double - на сколько метров от крайних точек можно так же считать объемы
    """
    
    @classmethod
    def get_standart(cls, tube_dict):
        """
            tube_dict_sample = {
                'tube_points': [[0, 0.023], [1, 0.023], [1.3, 0.023*0.5], [1.5, 0.023*0.5]]
            }        
        """
        if 'tube_points' not in tube_dict:
            raise Exception(f'В функции get_tube(**tube_dict) неправильные аргументы tube_dict={tube_dict}')
        xs, ds = [], []
        for x, d in tube_dict['tube_points']:
            xs.append(x)
            ds.append(d)
        return cls(xs, ds)
    
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
    
    cpdef void fill_ds(self, double[:] xs, double[:] ds):
        self.d.fill_vs(xs, ds)

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

    def get_x_right(self):
        return self.d.xs[-1]

    cpdef void fill_dsdx(self, double[:] xs, double[:] dsdx):
        cdef size_t i
        cdef double dx, si ,si1
        si = self.s.get_v(xs[0])
        for i in range(dsdx.shape[0]):
            si1  = self.s.get_v(xs[i+1])
            dx = xs[i+1] - xs[i]
            dsdx[i] = (si1-si)/dx
            si = si1

    cpdef double[:] get_dsdx(self, double[:] xs):
        """
        double[:] get_dsdx(self, double[:] xs)
        возвращает производные площадей трубы в точках xs, 
        НО это производные (Sx_i+1 - Sx_i)/dx, а не самой трубы. 
        Так же len(result) == len(xs) - 1
        """
        cdef double[:] dsdx = np.empty(xs.shape[0]-1, dtype=np.double)
        self.fill_dsdx(xs, dsdx)
        return dsdx
    
    cpdef void fill_W(self, double[:] xs, double[:] W):
        cdef size_t i
        cdef double s1, s2, x1, x2
        x1 = xs[0]
        s1 = self.s.get_v(x1)
        for i in range(W.shape[0]):
            x2 = xs[i+1]
            s2 = self.s.get_v(x2)
            W[i] = (s1 + s2 + sqrt(s1 * s2)) * (x2-x1)/ 3
            x1 = x2
            s1 = s2
        # cdef size_t i
        # cdef double wi ,wi1
        # wi = self.w.get_v(xs[0])
        # for i in range(W.shape[0]):
        #     wi1  = self.w.get_v(xs[i+1])
        #     W[i] = wi1 - wi
        #     wi = wi1        

    cpdef double[:] get_W(self, double[:] xs):
        """
        double[:] get_W(self, double[:] xs)
        метод возвращает объемы трубы, находящиеся между точками
        len(result) == len(xs) - 1
        """
        cdef double[:] res = np.empty(xs.shape[0]-1, dtype=np.double)
        self.fill_W(xs, res)
        return res  

    cpdef void fill_S(self, double[:] xs, double[:] S):
        self.s.fill_vs(xs, S)
        
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

    cpdef double get_x1(self, double x2, double w):
        """
        double get_x2(self, double x2, double w)
        возвращает координату x1, такую, что объем трубы между x1 и x2 равен w
        """
        cdef double w2 = self.w.get_v(x2)
        return self.w_reverse.get_v(w2 - w)
    
    cpdef double get_W_between(self, double x1, double x2):
        """
        double get_W_between(self, double x1, double x2)
        Возвращает объем трубы между точками x1 и x2
        """
        return self.w.get_v(x2) - self.w.get_v(x1)

    def plot(self, fig, ax, **kwargs):
        """Отрисовать трубу 

            fig, ax = plt.subplots()
            
            y0 = kwargs.get('y0', 0)
            color_muzzle = kwargs.get('color_muzzle', 'black')
            lw_muzzle = kwargs.get('lw_muzzle', 2)
            lw_osei = kwargs.get('lw_osei', 1)
            marker = kwargs.get('marker', 'o')
            markersize = kwargs.get('markersize', 5)
        """
        y0 = kwargs.get('y0', 0)
        xs, ys = np.array(self.get_xs()), np.array(self.get_ds())/2
        color_muzzle = kwargs.get('color_muzzle', 'black')
        lw_muzzle = kwargs.get('lw_muzzle', 2)
        lw_osei = kwargs.get('lw_osei', 1)
        marker = kwargs.get('marker', 'o')
        markersize = kwargs.get('markersize', 5)
        ax.plot(xs, ys+y0, color=color_muzzle, lw=lw_muzzle, marker=marker, markersize=markersize)
        ax.plot(xs[[0,0]], [y0, ys[0]+y0], color='black', lw=lw_osei)
        ax.plot(xs[[-1,-1]], [y0, ys[-1]+y0], color='black', lw=lw_osei)
        ax.plot(xs[[0,-1]], [y0,y0], color='black', lw=lw_osei, ls='-.', marker=marker, markersize=markersize)


    

class BorderSimple(object):
    border_dict_sample = {
    'mass': 0.1, # kg
    'p_f': 100e6, #Pa
    't_init': -1,
    'w_border': 0.005, #m
    }
#     @staticmethod
#     def connect(b_left: BorderSimple, b_right: BorderSimple, f_destroy=0, l_0=None):
#         if l_0 is None:
#             l_0 = b_right.get_x_left() - b_left.get_x_right() 
        
    @classmethod
    def get_standart(cls, border_dict):
        return cls(None, None, **border_dict)
        
    def __init__(self, 
                 lr_left, 
                 lr_right, 
                 mass, 
                 p_f, 
                 t_init=-1, 
                 w_border=0, **kwargs):
        
        self.lr_left= lr_left
        self.lr_right = lr_right
        self.mass = mass
        self.p_f = p_f
        self.t_init = t_init
        self.w_border = w_border
        self.color_4_plot = 'black'
        
    def __str__(self):
        return f'BorderSimple(mass={self.mass}, p_f={self.p_f}, t_init={self.t_init}, w_border={self.w_border})'
    
    def get_p_s_left(self):
        if self.lr_left is not None:
            return self.lr_left.get_p_right(), self.lr_left.S[-1]
        return 0, 0
        
    def get_p_s_right(self):
        if self.lr_right is not None:
            return self.lr_right.get_p_left(), self.lr_right.S[0]
        return 0, 0
    
    def get_right_layer_x1(self):
        return self.lr_left.xs_borders[-1] + self.w_border
    
    def get_x_left(self):
        if self.lr_left is not None:
            return self.lr_left.xs_borders[-1]
        if self.lr_right is not None:
            return self.get_x_right() - self.w_border
        return 0
    
    def get_x_right(self):
        if self.lr_right is not None:
            return self.lr_right.xs_borders[0]
        if self.lr_left is not None:
            return self.get_x_left() + self.w_border
        return 0
    
    def get_v0(self):
        v0 = None
        if self.lr_left is not None:
            v0 = self.lr_left.Vs_borders[-1]
        if self.lr_right is not None:
            if v0 is not None:
                if abs(self.lr_right.xs_borders[0] - self.get_right_layer_x1()) > 1e-5:
                    raise AttributeError(f'Слои не синхронизированы по координатам!\n{self.lr_left}\n{self.lr_right}')
                vr = self.lr_right.Vs_borders[0]
                if abs(vr-v0) > 1e-10:
                    raise AttributeError(f'Слои не синхронизированы по скорости!\n{self.lr_left}\n{self.lr_right}')
            v0 = self.lr_right.Vs_borders[0]
        return v0
    
    def get_v0_f0_a0(self):
        v0 = self.get_v0()
        (pl, sl), (pr, sr) = self.get_p_s_left(), self.get_p_s_right()
        f0 = pl*sl - pr*sr
        if abs(v0) < 1e-13 and abs(pl - pr) < self.p_f:
            a0 = 0
        else:
            a0 = f0/self.mass
        return v0, f0, a0
    
    def update_Vs(self, tau: float):
        v0, f0, a0 = self.get_v0_f0_a0()
        v1 = v0 + tau*a0
        if self.lr_left is not None:
            self.lr_left.Vs_borders[-1] = v1
        if self.lr_right is not None:
            self.lr_right.Vs_borders[0] = v1
            
    def update_xs(self):
        if self.lr_left is not None and self.lr_right is not None:
            self.lr_right.xl = self.lr_left.xs_borders[-1] + self.w_border
            self.lr_left.xr = self.lr_left.xs_borders[-1]
        elif self.lr_left is None:
            self.lr_right.xl = self.lr_right.xs_borders[0]
        elif self.lr_right is None:
            self.lr_left.xr = self.lr_left.xs_borders[-1]
            
    def plot(self, fig, ax, **kwargs):
        """
        fig, ax = plt.subplots()

        n_points = kwargs.get('n_points', 10)
        y0 = kwargs.get('y0', 0)
        color=kwargs.get('color', self.color_4_plot)
        lw = kwargs.get('lw', 2)
        kwargs.get('plot_tube', False)
        """
        from matplotlib.patches import Polygon
        n_points = kwargs.get('n_points', 10)
        y0 = kwargs.get('y0', 0)
        color=kwargs.get('color', self.color_4_plot)
        lw = kwargs.get('lw', 2)
        if self.lr_left is not None:
            tube = self.lr_left.tube
        elif self.lr_right is not None:
            tube = self.lr_right.tube
        else:
            raise AttributeError(f'Нет слоев( ')
        if kwargs.get('plot_tube', False):
            tube.plot(fig, ax, **kwargs)
        if abs(self.w_border) < 1e-9:
            x = self.get_x_left()
            ax.plot([x, x], [y0, tube.get_d(x)/2+y0], color=color, lw=lw)
            return
        ix = np.linspace(self.get_x_left(), self.get_x_right(), n_points)
        iy = np.array(tube.get_ds(ix))/2 + y0
        verts = [(ix[0], y0)] + list(zip(ix, iy)) + [(ix[-1], y0)]
        poly = Polygon(verts, facecolor='1', edgecolor='0.1', hatch=r'//')
        ax.add_patch(poly)
        
    def copy(self):
        return BorderSimple(lr_left=self.lr_left, 
                 lr_right=self.lr_right, 
                 mass=self.mass, 
                 p_f=self.p_f, 
                 t_init=self.t_init, 
                 w_border=self.w_border)
    