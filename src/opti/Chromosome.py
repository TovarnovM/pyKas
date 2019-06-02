import numpy as np
import random as rnd
from math import sqrt
from functional import seq


class DRange(object):
    """
    Класс описывает непрерывный интервал вещественных чисел [a;b], используемый как ген в хромосоме,
    а так же инструменты для кроссовера, мутации, интерполяции значений в интервале
    
    Атрибуты (они же поля):
        tail_part   показывает какая часть от нормального распределения (хвост) находится вне интервала родителей 
                    при кроссовере. НАПРИМЕР значения гена у родителей 1 и 7, tail_part = 0.5, значит параметры 
                    нормального распределения при розыгрыше значения гена потомка будут: МО = 4, СКО = 2
        mutate_part показывает какую долю интервала составляет значение СКО нормального распределения при розыгрыше
                    мутации. НАПРИМЕР значения гена у родителя 1, mutate_part = 0.5, интервал гера [a;b] = [-6;6].
                    значит параметры нормального распределения при розыгрыше значения гена мутанта будут: МО = 1, СКО = 1
        cross_f     функция кроссовера (по-умолчанию с использованием нормального распределения и параметра tail_part)
        mutate_f    функция мутации (по-умолчанию с использованием нормального распределения и параметра mutate_part)
        n_interv    количество единичных отрезков на интервале (необходимо для нормирования при кодировании)
        normalizable флаг, определяющий участвыует ли этот ген в параметрическом вещественном векторе (для детерминированной оптимизации)
    """

    tail_part = 0.15
    mutate_part = 0.5
    normalizable = True

    def __repr__(self):
        return f"DRange({self.a}, {self.b}, name='{self.name}')"

    def get_rnd_value(self):
        """
        функция генерации случайного значения в интервале
        """
        return rnd.uniform(self.a, self.b)

    def __init__(self, a, b, name='drange'):
        """
        a    левая (или правая) граница интервала гена
        b    правая(или левая)  граница интервала гена
        name имя гена
        """
        self.a = min(a, b)
        self.b = max(a, b)
        self.cross_f = self.cross_norm
        self.mutate_f = self.mutate_norm
        self.name = name
        self.n_interv = 1

    def get_vec_val(self, value):
        """
        метод расчета нормализованного значения (для вещественного вектора)
        :param value: значение гена
        :return: нормализованное значение гена
        """
        if value > self.b:
            value = self.b
        elif value < self.a:
            value = self.a
        return (value - self.a) / (self.b - self.a) * self.n_interv

    def get_vec_val_max(self):
        """
        максимальное значение нормализованного значения
        :return: максимальное нормализованное значение
        """
        return self.get_vec_val(self.b)

    def get_vec_val_min(self):
        """
         минимальное значение нормализованного значения
        :return: минимальное нормализованное значение
        """
        return self.get_vec_val(self.a)

    def validate_value(self, value):
        if value > self.b:
            return self.b
        elif value < self.a:
            return self.a
        return value

    def get_val_from_vec(self, value, validate=True):
        v = self.a + (self.b - self.a) * value / self.n_interv
        return self.validate_value(v) if validate else v

    def cross(self, x1, x2):
        """
        Кроссовер двух генов
        
        x1     значение гена первого родителя
        x2     значение гена второго родителя
        
        return значение гена потомка
        """
        return self.cross_f(x1, x2)

    def mutate(self, mu):
        """
        Мутация гена
        
        mu     значение гена до мутации
        
        return значение гена после мутации
        """
        return self.mutate_f(mu)

    def cross_norm(self, x1, x2):
        """
        Кроссовер двух генов (с использованием нормального распределения)
        
        x1     значение гена первого родителя
        x2     значение гена второго родителя
        
        return значение гена потомка
        """
        mu = 0.5 * (x1 + x2)
        xmin = min(x1, x2)
        # xmax = max(x1,x2)
        sko = (mu - xmin) / (3 - 3 * self.tail_part)
        x = rnd.gauss(mu, sko)
        while x < self.a or x > self.b:
            x = rnd.gauss(mu, sko)
        return x

    def mutate_norm(self, mu):
        """
        Мутация гена (с использованием нормального распределения)
        
        mu     значение гена до мутации
        
        return значение гена после мутации
        """
        sko = (self.b - self.a) / 6 * self.mutate_part
        x = rnd.gauss(mu, sko)
        while x < self.a or x > self.b:
            x = rnd.gauss(mu, sko)
        return x

    def interp(self, x1, x2, t):
        """
        Функция интерполяции значения между двух генов 
        
        x1     значение гена первого
        x2     значение гена второго
        t      параметр интерполяции. Наример при t = 0 return x1, при t = 1 return x2, t = 0.5 return (x1+x2)/2
        
        return значение интерполированнного гена
        """
        return (1 - t) * x1 + t * x2


class SRange(object):
    """
    Класс описывает конечный дискретный набор значений, используемый как ген в хромосоме,
    а так же инструменты для кроссовера, мутации, интерполяции этих значений

    normalizable флаг, определяющий участвыует ли этот ген в параметрическом вещественном векторе (для детерминированной оптимизации)
    """

    normalizable = False

    def __init__(self, variants, name='srange'):
        """
        variants  набор значений гена (список)
        name      имя гена
        """
        self.variants = variants
        self.name = name

    def get_vec_val(self, value):
        return self.variants.index(value)

    def get_vec_val_max(self):
        return len(self.variants) - 1

    def get_vec_val_min(self):
        return 0

    def get_val_from_vec(self, value, validate=True):
        return self.variants[int(value)]

    def __repr__(self):
        return f"SRange({repr(self.variants)}, name='{self.name}')"

    def get_rnd_value(self):
        """
        функция генерации случайного значения
        """
        return rnd.choice(self.variants)

    def cross(self, var1, var2):
        """
        Кроссовер двух генов
        
        var1    значение гена первого родителя
        var2    значение гена второго родителя
        
        return значение гена потомка
        """
        return rnd.choice([var1, var2])

    def mutate(self, var1):
        """
        Мутация гена
        
        var1   значение гена до мутации
        
        return значение гена после мутации
        """
        val = self.get_rnd_value()
        while val == var1:
            val = self.get_rnd_value()
        return val

    def interp(self, x1, x2, t):
        """
        Функция интерполяции значения между двух генов 
        
        x1     значение гена первого
        x2     значение гена второго
        t      параметр интерполяции. Наример при t = 0 return x1, при t = 1 return x2
        
        return значение интерполированнного гена
        """
        return x1 if t <= 0.5 else x2


class IRange(object):
    """
    Класс описывает конечный дискретный набор УПОРЯДОЧЕННЫХ значений, которые могут быть интерпритированны как 
    диапазон дискретных значений, используемый как ген в хромосоме,
    а так же инструменты для кроссовера, мутации, интерполяции этих значений 
     Атрибуты (они же поля):
        tail_part   показывает какая часть от нормального распределения (хвост) находится вне интервала родителей 
                    при кроссовере. НАПРИМЕР значения гена у родителей 1 и 7, tail_part = 0.5, значит параметры 
                    нормального распределения при розыгрыше значения гена потомка будут: МО = 4, СКО = 2
        mutate_part показывает какую долю интервала составляет значение СКО равномерного распределения при розыгрыше
                    мутации. НАПРИМЕР значения гена у родителя 1, mutate_part = 3, интервал гена [a;b] = [-6;6].
                    значит параметры равномерного распределения при розыгрыше значения гена мутанта будут: х_min = -2, х_max = 4
    """
    tail_part = 0.5
    mutate_part = 2
    normalizable = True

    def __init__(self, variants, name='irange'):
        """
        variants  набор значений гена (список)
        name      имя гена
        """
        self.variants = variants
        self.name = name
        self.n_interv = len(variants)

    def get_vec_val(self, value):
        return self.variants.index(value) / (self.n_interv - 1)

    def get_val_from_vec(self, value, validate=True):
        ind = int(value * self.n_interv)
        if ind < 0:
            ind = 0
        elif ind >= len(self.variants):
            ind = len(self.variants) - 1
        return self.variants[ind]

    def get_vec_val_max(self):
        return self.get_vec_val(self.variants[-1])

    def get_vec_val_min(self):
        return self.get_vec_val(self.variants[0])

    def __repr__(self):
        return f"IRange({repr(self.variants)}, name='{self.name}')"

    def get_rnd_value(self):
        """
        функция генерации случайного значения
        """
        return rnd.choice(self.variants)

    def cross(self, var1, var2):
        """
        Кроссовер двух генов
        
        var1    значение гена первого родителя
        var2    значение гена второго родителя
        
        return значение гена потомка
        """
        x1 = self.variants.index(var1)
        x2 = self.variants.index(var2)

        mu = 0.5 * (x1 + x2)
        xmin = min(x1, x2)
        # xmax = max(x1,x2)
        sko = (mu - xmin) / (3 - 3 * self.tail_part)
        x = rnd.gauss(mu, sko)
        while x < 0 or x > len(self.variants) - 1:
            x = rnd.gauss(mu, sko)
        return self.variants[round(x)]

    def mutate(self, var1):
        """
        Мутация гена
        
        var1   значение гена до мутации
        
        return значение гена после мутации
        """
        ind = self.variants.index(var1)
        ind_min = ind - self.mutate_part
        ind_max = ind + self.mutate_part
        ind_min = 0 if ind_min < 0 else ind_min
        ind_max = len(self.variants) - 1 if ind_max >= len(self.variants) else ind_max
        indChild = rnd.randrange(ind_min, ind_max + 1)
        while self.variants[indChild] == var1:
            indChild = rnd.randrange(ind_min, ind_max + 1)
        return self.variants[indChild]

    def interp(self, x1, x2, t):
        """
        Функция интерполяции значения между двух генов 
        
        x1     значение гена первого
        x2     значение гена второго
        t      параметр интерполяции. Наример при t = 0 return x1, при t = 1 return x2
        
        return значение интерполированнного гена
        """
        ind1 = self.variants.index(x1)
        ind2 = self.variants.index(x2)
        return round((1 - t) * ind1 + t * ind2)


class ChromoController(object):
    """
    Класс описывает хромосому постоянной структуры
    а так же инструменты для кроссовера, мутации, интерполяции хромосом
    
    Атрибуты (они же поля):
        gr          константа для поиска приемлимых значений
        ranges_dict словарь с описанием генов типа { name1:DRange,name2:SRange,name3:IRange, ... }, гдены 
                    должны поддерживать методы get_rnd_value, cross, mutate, interp
        constraints список фунций для определения пригодности хромосомы, функции вида
                    f(chr) -> float, где chr - хромосома (создается функциями get_chromo, cross, mutate, interp, 
                    find_zero_golden_method). Если возврацаемое значение >= 0, то хромосома приемлимая
        ranges_list список с описанием генов типа (возможно не нужен)
    """

    gr = (sqrt(5) + 1) / 2

    def __init__(self, ranges, constraints=None):
        """
        ranges      список описаний генов (гены должны иметь разные имена!!!)
        constraints список функций для определения пригодности хромосомы, функции вида
                    f(chr) -> float, где chr - хромосома (создается функциями get_chromo, cross, mutate, interp, 
                    find_zero_golden_method). Если возврацаемое значение >= 0, то хромосома приемлимая
        """
        self.ranges_dict = {r.name: r for r in ranges}
        self.constraints = constraints
        self.ranges_list = ranges
        self.normalizable_list = seq(ranges).filter(lambda r: r.normalizable).to_list()
        self.non_normalizable_list = seq(ranges).filter(lambda r: not r.normalizable).to_list()

    def get_vec_param(self, chromo):
        """
        Метод кодирования хромосомы (для normalizable-генов)
        :param chromo: хромосома
        :return: np.array
        """
        res = np.zeros(len(self.normalizable_list))
        for i, ng in enumerate(self.normalizable_list):
            gname = ng.name
            res[i] = ng.get_vec_val(chromo[gname])
        return res

    def get_vec_struct(self, chromo):
        """
        Метод кодирования гена для
        :param chromo:
        :return:
        """
        res = np.zeros(len(self.non_normalizable_list))
        for i, ng in enumerate(self.non_normalizable_list):
            gname = ng.name
            res[i] = ng.get_vec_val(chromo[gname])
        return res

    def get_chromo_from_vecs(self, vec_param, vec_struct, validate=True):
        chromo = self.get_chromo()
        for i, g in enumerate(self.normalizable_list):
            chromo[g.name] = g.get_val_from_vec(vec_param[i], validate)
        for i, g in enumerate(self.non_normalizable_list):
            chromo[g.name] = g.get_val_from_vec(vec_struct[i])
        return chromo

    def __repr__(self):
        return f"ChromoController({repr(self.ranges_list)},constraints = {repr(self.constraints)})"

    def __getitem__(self, key):
        if key in self.ranges_dict:
            return self.ranges_dict[key]
        return self.ranges_list[key]

    def __iter__(self):
        return iter(self.ranges_dict)

    def get_chromo(self):
        """
        метод получения новой хромосомы со случайными, но совокупно приемлимыми значениями генов
        """
        ch = {key: self[key].get_rnd_value() for key in self}
        if not (self.constraints is None):
            while self.f_constr(ch) < 0:
                ch = {key: self[key].get_rnd_value() for key in self}
        return ch

    def cross(self, chromo1, chromo2):
        """
        метод проведения кроссовера двух хромосом
        
        chromo1, chromo2  хромосомы родителей
        
        return            хромосома-потомок
        """
        child = {key: self[key].cross(chromo1[key], chromo2[key]) for key in self}
        if not (self.constraints is None):
            if self.f_constr(child) < 0:
                child = self.find_zero_golden_method(rnd.choice([chromo1, chromo2]), child)
        return child

    def mutate(self, chromo, keys):
        """
        метод проведения мутации хромосомы 
        
        chromo  хромосома-родитель
        keys    имена мутирующих генов
        
        return  хромосома-мутант
        """
        mutant = dict(chromo)
        for key in keys:
            mutant[key] = self[key].mutate(chromo[key])
        if not (self.constraints is None):
            if self.f_constr(mutant) < 0:
                mutant = self.find_zero_golden_method(chromo, mutant)
        return mutant

    def validate_vec_param(self, v0, v_try, s0, eps=1E-10):
        """
        возвращает отвалидированный v_try и ответ, отличается ли отвалидированный v_try от v0
        :param v0: vec_param0
        :param v_try: vec_param
        :param s0: vec_struct
        :param eps: точность
        :return:
        """
        cr0 = self.get_chromo_from_vecs(v0, s0)
        cr_try = self.get_chromo_from_vecs(v_try, s0)
        if self.f_constr(cr_try) < 0:
            cr_try = self.find_zero_golden_method(cr0, cr_try, eps)
        v_try = self.get_vec_param(cr_try)

        return v_try, max(np.abs(v_try - v0)) > eps

    def interp(self, chr1, chr2, t):
        return {key: self[key].interp(chr1[key], chr2[key], t) for key in self}

    def f_constr(self, chromo):
        if not self.constraints:
            return 43
        return min([c(chromo) for c in self.constraints])

    def find_zero_golden_method(self, chr1, chr2, eps=1E-4, n_max=13):
        a, b = 0, 1
        c = b - (b - a) / self.gr
        d = a + (b - a) / self.gr
        fc = abs(self.f_constr(self.interp(chr1, chr2, c)))
        fd = abs(self.f_constr(self.interp(chr1, chr2, d)))
        for i in range(n_max):
            if abs(c - d) < eps:
                break
            if fc < fd:
                b = d
                d = c
                fd = fc
                c = b - (b - a) / self.gr
                fc = abs(self.f_constr(self.interp(chr1, chr2, c)))
            else:
                a = c
                c = d
                fc = fd
                d = a + (b - a) / self.gr
                fd = abs(self.f_constr(self.interp(chr1, chr2, d)))
        return self.interp(chr1, chr2, (b + a) / 2)

    def get_vec_param_min(self):
        lst = [r.get_vec_val_min() for r in self.normalizable_list]
        return np.array(lst)

    def get_vec_param_max(self):
        lst = [r.get_vec_val_max() for r in self.normalizable_list]
        return np.array(lst)

    def get_vec_struct_min(self):
        lst = [r.get_vec_val_min() for r in self.non_normalizable_list]
        return np.array(lst)

    def get_vec_struct_max(self):
        lst = [r.get_vec_val_max() for r in self.non_normalizable_list]
        return np.array(lst)

    def chromo_vel(self, chr0, vel_param):
        """
        делает шаг по скорости, возвращает новую хромосому, удовлетворяющую всем условиям
        :param chr0: хромосома0
        :param vel_param: np.array "скорсоть"
        :return: хромосома = хромосома0 + "скорсоть"
        """
        v0_param = self.get_vec_param(chr0)
        v0_struct = self.get_vec_struct(chr0)

        chr1 = self.get_chromo_from_vecs(v0_param + vel_param, v0_struct)
        if not (self.constraints is None):
            if self.f_constr(chr1) < 0:
                chr1 = self.find_zero_golden_method(chr0, chr1)
        return chr1

    def get_n_4grad(self):
        return len(self.normalizable_list) * 2

    def get_border(self, c0, grad0, inc_f_constr=True):
        v0 = c0['vec_param']
        v_min = self.get_vec_param_min()
        v_max = self.get_vec_param_max()
        to_maxy = grad0 > 0
        to_miny = grad0 < 0
        l = np.zeros_like(v0)
        l[to_maxy] = (v_max - v0)[to_maxy]
        l[to_miny] = (v_min - v0)[to_miny]
        ts = l[grad0 != 0] / grad0[grad0 != 0]
        v1 = v0 + min(ts) * grad0
        if inc_f_constr:
            s0 = c0['vec_struct']
            c1 = self.get_chromo_from_vecs(v1, s0)
            if self.f_constr(c1) < 0:
                c1 = self.find_zero_golden_method(c0, c1)
            return c1['vec_param']
        return v1

    def get_valid_vec_param(self, vec_param):
        res = np.zeros_like(vec_param)
        for i, ng in enumerate(self.normalizable_list):
            res[i] = ng.get_vec_val(ng.get_val_from_vec(vec_param[i]))
        return res

    def get_f_constr_grad(self, wc0, dh, normilize=True):
        if not self.constraints:
            raise AttributeError("В хромосомах нет функциональных ограничений")
        s0 = wc0['vec_struct']
        v0 = wc0['vec_param']

        grad = np.zeros(len(self.normalizable_list))

        def calc_grad_i(i):
            v1 = np.copy(v0)
            v1[i] += dh
            vm1 = np.copy(v0)
            vm1[i] -= dh
            c1 = self.get_chromo_from_vecs(v1, s0)
            cm1 = self.get_chromo_from_vecs(vm1, s0)
            grad_i = (self.f_constr(c1) - self.f_constr(cm1)) / 2 / dh
            return grad_i

        for i in range(len(grad)):
            grad[i] = calc_grad_i(i)

        if normilize:
            grad = grad / np.linalg.norm(grad)

        return grad

    def get_bounds_4scipy(self):
        min_vec = self.get_vec_param_min()
        max_vec = self.get_vec_param_max()
        return [tp for tp in zip(min_vec, max_vec)]

    def get_constr_4scipy(self, chromo):
        if not self.constraints:
            return ()
        s0 = self.get_vec_struct(chromo)

        def func_factory(constr):
            def fun(x):
                chr = self.get_chromo_from_vecs(x, s0, False)
                res = constr(chr)
                return res

            return fun

        return [{'type': 'ineq', 'fun': func_factory(constr)} for constr in self.constraints]


def get_test_chr_contr():
    dr1 = DRange(-1, 3, 'dr1')
    dr2 = DRange(1, 3, 'dr2')
    dr3 = DRange(10, 30, 'dr3')
    ir1 = IRange([1, 2, 3, 4, 5, 6, 7], 'ir1')
    sr1 = SRange(['one', 'two', '3'], 'sr1')
    sr2 = SRange(['1', '2', '33'], 'sr2')

    chr_contr = ChromoController([dr1, dr2, dr3, ir1, sr1, sr2])
    return chr_contr


def main():
    dr1 = DRange(-1, 3, 'dr1')
    dr1.n_interv = 4.0

    v=dr1.get_vec_val(1)
    d=dr1.get_val_from_vec(v)
    print(v)
    print(d)
    # dr2 = DRange(1, 3, 'dr2')
    # dr3 = DRange(10, 30, 'dr3')
    # ir1 = IRange([1, 2, 3, 4, 5, 6, 7], 'ir1')
    # sr1 = SRange(['one', 'two', '3'], 'sr1')
    # sr2 = SRange(['1', '2', '33'], 'sr2')
    #
    # chr_contr = ChromoController([dr1, dr2, dr3, ir1, sr1, sr2], constraints=[lambda ch: ch['dr1']])
    # cr1 = chr_contr.get_chromo()
    # cr2 = chr_contr.get_chromo()
    # cr3 = chr_contr.cross(cr1, cr2)
    # print(cr1, cr2, cr3)
    #
    # v_p = chr_contr.get_vec_param(cr1)
    # v_s = chr_contr.get_vec_struct(cr1)
    # print(v_p, v_s)
    # cr11 = chr_contr.get_chromo_from_vecs(v_p, v_s)
    # v_p11 = chr_contr.get_vec_param(cr11)
    # v_s11 = chr_contr.get_vec_struct(cr11)
    # print(v_p == v_p11, v_s == v_s11)


if __name__ == '__main__':
    main()
