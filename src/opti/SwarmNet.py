from multiprocessing import RLock
import numpy as np
import random as rnd
from math import ceil
from heapq import *

def hash_arr(arr):
    return arr.tostring()

class SurfPoint(object):
    def __init__(self, dist, ind, vec_param):
        self.vec_param = vec_param
        self.dist = dist
        self.ind = ind

    def __lt__(self, other):
        return self.dist < other.dist

class SwarmNet(object):
    """
    Класс для детерминированного параллельного поиска локального максимума
    """
    def __init__(self, op, dx=1E-3):
        """

        :param op: начальная точка поиска
        :param dx: шаг по vec_param
        """
        self.op = op
        wchr0, self.fit_max = op.get_best(fit_too=True)
        self.v0 = wchr0['vec_param']
        self.s0 = wchr0['vec_struct']
        self.ind_best = np.zeros_like(self.v0, dtype=np.int)
        self.dx = dx
        self.locker = RLock()
        self.all_dict = {}
        self.surf_p_heap = []

        key = hash_arr(self.ind_best)
        di = {'ind': self.ind_best }
        wchr_id = wchr0['id']
        wchr0['dop_info'] = di
        self.all_dict[key] = {
            'fitness': self.fit_max,
            'id': wchr_id,
            'ind': self.ind_best,
            'vec_param': self.v0 }

        neibs_inds = self.get_neibs_indexes(self.ind_best)
        v_center = self.all_dict[hash_arr(self.ind_best)]['vec_param']
        for neib_ind in neibs_inds:
            if self.is_in_dict(neib_ind):
                continue
            neib_v, isgood = self.get_good_neib_v(v_center, neib_ind)
            if isgood:
                self.add_to_surf_heap(neib_ind, neib_v)
        if not self.fit_max:
            self.add_to_surf_heap(self.ind_best, self.v0)

    def dist_func(self, ind):
        """
        Метод получения расстояния между решением под индексом ind и до самого лучшего решения
        Чем ближе точка находится к лучшему решению, тем вероятнее ее посчитают
        :param ind: индеск точки
        :return: дистанция
        """
        return np.linalg.norm(ind - self.ind_best)

    def add_to_surf_heap(self, ind, vec_param):
        """
        метод добавления новой потенциальной точки в "поверхность" выбора
        :param ind: индекс новой точки
        :param vec_param:
        :return: none
        """
        dist = self.dist_func(ind)
        heappush(self.surf_p_heap, SurfPoint(dist, ind, vec_param))

    def get_v(self, indexes):
        """
        получить vec_param по индексу
        :param indexes: int or nd-int array of delta xs
        :return:
        """
        v = np.copy(self.v0)
        dv = np.array(indexes, dtype=np.int) * self.dx
        v = v + dv
        return v

    def get_good_neib_v(self, v_center, neib_ind):
        """
        получить отвалидированный vec_param
        :param ind:
        :param neib_ind:
        :return: (neib_v, true/false); true=good neib_v
        """
        neib_v = self.get_v(neib_ind)
        return self.op.chr_contr.validate_vec_param(v_center, neib_v, self.s0)

    def get_neibs_indexes(self, ind):
        """
        Получить список соседей (неотвалидированный)
        :param ind: индеск центра
        :return: [neib_ind, ..]
        """
        neibs = []
        for i in range(len(ind)):
            n1 = np.array(ind, dtype=np.int)
            n1[i] += 1
            neibs.append(n1)
            n2 = np.array(ind, dtype=np.int)
            n2[i] -= 1
            neibs.append(n2)
        return neibs

    def is_in_dict(self, ind):
        """

        :param ind: ndarray, int
        :return: true/false
        """
        key = hash_arr(ind)
        return key in self.all_dict

    def add_to_all_dict(self, ind, vec_param):
        """
        ДОбавляет новую точку в мега-словарь
        :param ind:
        :param vec_param:
        :return:
        """
        key = hash_arr(ind)
        if key in self.all_dict:
            return
        chromo = self.op.chr_contr.get_chromo_from_vecs(vec_param, self.s0)
        di = {'ind': ind }
        wchr_id = self.op.add_new_chromo(chromo, dop_info=di)
        self.all_dict[key] = {
            'fitness': None,
            'id': wchr_id,
            'ind': ind,
            'vec_param': self.op[wchr_id]['vec_param'] }

    def get_wchr_4calc_thread_save(self):
        """
        главынй метод получения точки для расчета (потокобезопасный)
        :return: wchromo
        """
        res = None
        try:
            self.locker.acquire()
            res = self.get_wchr_4calc()
        finally:
            self.locker.release()
        return res


    def get_wchr_4calc(self):
        """
        главынй метод получения точки для расчета (потокоНЕбезопасный)
        :return: wchromo
        """
        surf_p = heappop(self.surf_p_heap)
        ind = surf_p.ind
        v_center = surf_p.vec_param
        self.add_to_all_dict(ind, v_center)

        wchr_id = self.all_dict[hash_arr(ind)]['id']
        wchr = self.op[wchr_id]

        neibs_inds = self.get_neibs_indexes(ind)
        rnd.shuffle(neibs_inds)
        for neib_ind in neibs_inds:
            if self.is_in_dict(neib_ind):
                continue
            neib_v, isgood = self.get_good_neib_v(v_center, neib_ind)
            if isgood:
                self.add_to_surf_heap(neib_ind, neib_v)

        return wchr

    def take_calced_wchr_thread_save(self, wchr, fit):
        """
        Потокобезопастный метод принятия посчитанного варианта
        :param wchr:
        :param fit:
        :return:
        """
        res = None
        try:
            self.locker.acquire()
            res = self.take_calced_wchr(wchr, fit)
        finally:
            self.locker.release()
        return res

    def take_calced_wchr(self, wchr, fit):
        """
        ПотокоНЕбезопастный метод принятия посчитанного варианта
        :param wchr:
        :param fit:
        :return:
        """
        if wchr['name'] != self.op.name:
            return
        wchr_id = wchr['id']
        self.op[wchr_id]['fitness'] = fit
        ind = np.array(wchr['dop_info']['ind'], dtype=np.int)
        value = self.all_dict[hash_arr(ind)]
        value['fitness'] = fit

        res = (not self.fit_max) or fit > self.fit_max

        if res:
            self.fit_max = fit
            self.ind_best = value['ind']
            for sp in self.surf_p_heap:
                dist = self.dist_func(sp.ind)
                sp.dist = dist
            heapify(self.surf_p_heap)

            n = len(self.ind_best)*2
            hl = len(self.surf_p_heap)
            if n < hl:
                h_new = []
                for _ in range(n):
                    heappush(h_new, heappop(self.surf_p_heap))
                self.surf_p_heap = h_new

        return res

    def on_improve(self):
        pass

    def get_fit(self, ind):
        """
        мтод получения фитнеса точки по индексу
        :param ind:
        :return:
        """
        key = hash_arr(ind)
        if key in self.all_dict:
            return self.all_dict[key]['fitness']
        else:
            return None

    def is_extremum(self, ind=None):
        """
        Уже всё?
        :param ind:
        :return:
        """
        if not ind:
            ind = self.ind_best
        fit = self.get_fit(ind)
        v_center = None
        if not fit:
            return False
        inds_neibs = self.get_neibs_indexes(ind)
        for ind_neib in inds_neibs:
            fit_neib = self.get_fit(ind_neib)
            if not fit_neib:
                if not v_center:
                    v_center = self.get_v(ind)
                neib_v, isgood = self.get_good_neib_v(v_center, ind_neib)
                if isgood:
                    return False
            elif fit_neib > fit:
                return False
        return True

    def get_best_cromo(self):
        """
        получить лучшую chromo
        :return:
        """
        key = hash_arr(self.ind_best)
        wchr_id = self.all_dict[key]['id']
        wchr = self.op[wchr_id]
        cromo = wchr['chromo']
        return cromo




