import random as rnd
import numpy as np
from functional import seq
from math import sqrt
from .OptiPerson import OptiPerson


def calc_diff_euclide(op1, op2, sel_sol_func):
    """
    функция вычисления Эвклидового параметрического расстояния между двумя решениями OptiPerson

    :param op1: OptiPerson
    :param op2: OptiPerson
    :param sel_sol_func: функция выбора решения
                         sel_sol_func(OptiPerson) -> OptiPerson.history[*]
    :return: double
    """

    def get_vecs(op):
        sol = sel_sol_func(op)
        return sol['vec_param'], sol['vec_struct']

    vp1, vs1 = get_vecs(op1)
    vp2, vs2 = get_vecs(op2)
    diff_par = np.square(vp1 - vp2)
    diff_struct = vs1 != vs2
    return sqrt(np.sum(diff_par) + np.sum(diff_struct))


def get_diff_matr(op_lst, sel_sol_func, diff_fuc):
    n = len(op_lst)
    matr = np.zeros((n, n))
    for i in range(n):
        op1 = op_lst[i]
        for j in range(i + 1, n):
            op2 = op_lst[j]
            matr[i, j] = diff_fuc(op1, op2, sel_sol_func)
            matr[j, i] = matr[i, j]
    return matr


def find_best_neibs(op_lst, neib_count):
    """
    Функция находит для каждого OptiPerson из списка op_lst лучшего соседа (из ближайших neib_count соседей)
    так же находит лучшую OptiPerson
    :param op_lst: список из OptiPerson
    :param neib_count: количество ближайших соседей, среди которых будут искаться лучшие
    :return: кортеж best_neibs_i - np.array из индексов лучших соседей
                    best_i - индекс лучшего решения вообще
    """
    best_fitnesses = np.array(seq(op_lst).map(lambda op: op.best_fitness).to_list())

    diff_matr = get_diff_matr(op_lst, lambda op: op.get_last(), calc_diff_euclide)

    def get_best_neib_i(my_index):
        row = diff_matr[my_index, :]
        neib_indexes = np.array(np.argsort(row)[1:neib_count + 1])
        neibs_fitns = best_fitnesses[neib_indexes]
        return np.argmax(neibs_fitns)

    best_neibs_i = [get_best_neib_i(i) for i in range(len(op_lst))]
    return best_neibs_i, np.argmax(best_fitnesses)


class StepUpSwarm(object):
    def __init__(self, chr_contr):
        self.chr_contr = chr_contr
        self.neib_count = 5
        self.b_I = 0.7298
        self.b_C = 1.49618
        self.b_S = self.b_C

    def step_up(self, pop, seed=None):
        """
        pop - словарь "имя - OptiPerson"
        """
        if seed:
            rnd.seed(seed)
            np.random.seed(seed)

        op_lst = seq(pop.values()).map(lambda op: op.copy()).to_list()
        best_neibs_i, best_i = find_best_neibs(op_lst, self.neib_count)
        vels0 = [op.vel for op in op_lst]
        x_i_lst = [op.get_last()['vec_param'] for op in op_lst]
        x_mybest_lst = [op.get_best()['vec_param'] for op in op_lst]
        x_best_lst = [x_mybest_lst[i] for i in best_neibs_i]
        vec_length = len(vels0[0])

        def calc_vel(v0, x_i, x_mybest, x_best):
            r1 = np.random.random(vec_length)
            r2 = np.random.random(vec_length)
            return v0 * self.b_I + self.b_C * r1 * (x_mybest - x_i) + r2 * self.b_S * (x_best - x_i)

        vel_lst = seq(zip(vels0, x_i_lst, x_mybest_lst, x_best_lst)).map(lambda tp: calc_vel(*tp)).to_list()
        for op, vel in zip(op_lst, vel_lst):
            last = op.get_last()
            chromo_new = self.chr_contr.chromo_vel(last['chromo'], vel)
            op.add_new_chromo(chromo_new)
            op.clear_swarm_hist()

        return OptiPerson.lst_to_dict(op_lst)