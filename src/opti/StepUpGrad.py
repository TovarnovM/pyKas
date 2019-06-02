from .StepUpSwarm import get_diff_matr
import numpy as np
from math import sqrt
from functional import seq
from math import ceil


def calc_diff_struct(op1, op2, sel_sol_func):
    """
    функция вычисления параметрического расстояния между двумя решениями OptiPerson,
    Если у них есть структурные параметры - то отличаи будут считаться только по ним

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
    if len(vs1) > 0:
        diff_struct = vs1 != vs2
        return sqrt(np.sum(diff_struct))
    diff_par = np.square(vp1 - vp2)
    return np.sum(diff_par)


def reduce_pop(pop_lst, n_uniquest, elite_top_perc=50):
    lst_elite = sorted(pop_lst, key=lambda op: op.get_best()['fitness'], reverse=True)
    n_part = int(len(lst_elite) * elite_top_perc / 100)
    if n_uniquest <= n_part:
        lst_elite = lst_elite[:n_part]
    mat = get_diff_matr(lst_elite, lambda op: op.get_last(), calc_diff_struct)
    n = len(lst_elite)
    inds = np.arange(n)
    loosers = np.repeat(True, n)
    elite = [0]
    for _ in range(n_uniquest - 1):
        loosers[elite] = False
        res = np.zeros_like(mat[0, loosers])
        for r in (mat[i, loosers] for i in elite):
            res += r
        next_elite = inds[loosers][np.argmax(res)]
        elite.append(next_elite)
    return [lst_elite[ei].copy() for ei in elite]


def get_grad(op, id_x0=None):
    wchr_gr = seq(op.history.values()) \
        .filter(lambda wchr: wchr['dop_info']) \
        .filter(lambda wchr: 'title' in wchr['dop_info']) \
        .filter(lambda wchr: wchr['dop_info']['title'] == '4grad') \
        .to_list()

    if len(wchr_gr) == 0:
        raise AttributeError("Нет таких окрестностей")

    if id_x0:
        wchr_gr = seq(wchr_gr) \
            .filter(lambda wchr: wchr['dop_info']['id_x0'] == id_x0) \
            .to_list()
    else:
        ids = seq(wchr_gr).map(lambda wchr: wchr['dop_info']['id_x0']).distinct().to_list()
        if len(ids) != 1:
            raise AttributeError("Слишком много разных окрестностей")

        id_x0 = ids[0]

    f0 = op[id_x0]['fitness']
    extremum = seq(wchr_gr).count(lambda wchr: wchr['fitness'] > f0) == 0
    grad = np.zeros_like(wchr_gr[0]['vec_param'])

    tups = []
    while len(wchr_gr) != 0:
        c1 = wchr_gr.pop()
        v1 = c1['vec_param']
        for i, c2 in enumerate(wchr_gr):
            v2 = c2['vec_param']
            dv = v2 - v1
            inds = np.argwhere(np.abs(dv) > 1E-10)
            if len(inds) != 1:
                continue
            c2 = wchr_gr.pop(i)
            tups.append((c1, c2, inds[0]))
            break
        else:
            raise AttributeError("Хромосомы для градиента так себе")

    def fill_one_number(cc1, cc2, ind):
        f1 = cc1['fitness']
        f2 = cc2['fitness']
        x1 = cc1['vec_param'][ind]
        x2 = cc2['vec_param'][ind]
        grad[ind] = (f1 - f2) / (x1 - x2)

    for tp in tups:
        fill_one_number(*tp)

    return grad, extremum, id_x0


def write_grad(op, grad, extremum, id_x0, del_chr_4grad=True):
    dop_info = {'title': 'grad',
                'id_x0': id_x0,
                'value': grad,
                'extremum': extremum}
    op[id_x0]['dop_info'] = dop_info
    if del_chr_4grad:
        del_chromos_4grad(op, id_x0)


def del_chromos_4grad(op, id_x0=None):
    chr_4grad = seq(op.history.values()) \
        .filter(lambda wc: wc['dop_info']) \
        .filter(lambda wchr: 'title' in wchr['dop_info']) \
        .filter(lambda wchr: wchr['dop_info']['title'] == '4grad') \
        .to_list()
    if id_x0:
        chr_4grad = seq(chr_4grad) \
            .filter(lambda wchr: wchr['dop_info']['id_x0'] == id_x0) \
            .to_list()
    for c in chr_4grad:
        op.remove(c['id'])


def del_chromos_4line(op, id_x0=None):
    chr_4grad = seq(op.history.values()) \
        .filter(lambda wc: wc['dop_info']) \
        .filter(lambda wchr: 'title' in wchr['dop_info']) \
        .filter(lambda wchr: wchr['dop_info']['title'] == '4line') \
        .to_list()
    if id_x0:
        chr_4grad = seq(chr_4grad) \
            .filter(lambda wchr: wchr['dop_info']['id_x0'] == id_x0) \
            .to_list()
    for c in chr_4grad:
        op.remove(c['id'])


def get_and_write_grad(op, id_x0=None, del_chr_4grad=True):
    grad, extremum, id_x0 = get_grad(op, id_x0)
    write_grad(op, grad, extremum, id_x0, del_chr_4grad)
    return grad, extremum, id_x0


def get_line_vecs(op, n, length, h_min=1E-5, id_x0=None, grad=None, f_constr_eps=1E-5, inc_f_constr=True):
    if not id_x0:
        c0_lst = seq(op.history.values()) \
            .filter(lambda wc: wc['dop_info']) \
            .filter(lambda wc: 'title' in wc['dop_info']) \
            .filter(lambda wc: wc['dop_info']['title'] == 'grad') \
            .sorted(lambda wc: wc['id']) \
            .to_list()
        id_x0 = c0_lst[-1]['id']
    c0 = op[id_x0]
    if not grad:
        grad = c0['dop_info']['value']
    chr_contr = op.chr_contr
    v0 = c0['vec_param']
    s0 = c0['vec_struct']
    g0 = grad / np.linalg.norm(grad)
    v1 = chr_contr.get_border(c0, g0, inc_f_constr)
    to_border_l = np.linalg.norm(v1 - v0)
    if to_border_l < h_min:
        # Если мы уперлись в границу
        if inc_f_constr and chr_contr.f_constr(c0['chromo']) < f_constr_eps:
            # Если мы уперлис в границу f_constr
            n0_constr = chr_contr.get_f_constr_grad(c0, h_min)
            p_back = v0 + n0_constr
            p_front = v0 + g0
            c_back = chr_contr.get_chromo_from_vecs(p_back, s0, False)
            c_front = chr_contr.get_chromo_from_vecs(p_front, s0, False)
            c1 = chr_contr.find_zero_golden_method(c_back, c_front)
            p1 = c1['vec_param']
            gr_new = p1 - v0
            gr_new_length = np.linalg.norm(gr_new)
            if gr_new_length < 1E-7:
                return None
            g0 = gr_new / gr_new_length
            points = get_line_vecs(op, n, length, h_min, id_x0, g0, f_constr_eps, inc_f_constr=False)
            if points:
                def valid_p(p_i):
                    c_i = chr_contr.get_chromo_from_vecs(p_i, s0, False)
                    if chr_contr.f_constr(c0['chromo']) > f_constr_eps:
                        return p_i
                    c_i_val = chr_contr.find_zero_golden_method(c_back, c_i)
                    return c_i_val['vec_param']

                return [valid_p(p_i) for p_i in points]
            else:
                return None

        else:
            # Если мы уперлись в границу глобальную
            p_front = v0 + g0
            p1 = chr_contr.get_valid_vec_param(p_front)
            gr_new = p1 - v0
            gr_new_length = np.linalg.norm(gr_new)
            if gr_new_length < 1E-7:
                return None
            g0 = gr_new / gr_new_length
            v1 = chr_contr.get_border(c0, g0, True)
            to_border_l = np.linalg.norm(v1 - v0)
            if to_border_l < h_min:
                return None

    length = length if length < to_border_l else to_border_l
    v1 = v0 + g0 * length
    h = length / n
    xs = np.linspace(0, length, num=n, endpoint=False, retstep=False) \
        if h > h_min \
        else np.arange(0, length, h_min)
    return [v1 - x * g0 for x in xs]


def linspace_vec(v1, v2, n):
    ts = np.linspace(0, 1, n + 1, endpoint=False)[1:]
    return [(1 - t) * v1 + t * v2 for t in ts]

def get_good_linspace_vec(v1, v2, n, h_min):
    l = np.linalg.norm(v1 - v2)
    h = l / (n + 1)
    if h < h_min:
        n = ceil(l / h_min)
        h = l / (n + 1)
    return linspace_vec(v1, v2, n), h


def get_line_max(op, h_min, n):
    line_lst = seq(op.history.values()) \
        .filter(lambda wc: wc['dop_info']) \
        .filter(lambda wc: 'title' in wc['dop_info']) \
        .filter(lambda wc: wc['dop_info']['title'] == '4line') \
        .sorted(lambda wc: wc['id'])
    line_lst = list(enumerate(line_lst))
    max_tp = seq(line_lst).max_by(lambda tp: tp[1]['fitness'])
    x0 = max_tp[1]['vec_param']
    if max_tp[0] == 0:
        xr = line_lst[1][1]['vec_param']
        xl = x0
    elif max_tp[0] == (len(line_lst) - 1):
        xl = line_lst[-2][1]['vec_param']
        xr = x0
    else:
        xr = line_lst[max_tp[0] + 1][1]['vec_param']
        xl = line_lst[max_tp[0] - 1][1]['vec_param']
        n = int(n / 2)



class StepUpGrad(object):
    def __init__(self):
        self.h = 1E-5

    def prep_4grad(self, op):
        best = op.get_best()
        best_id = best['id']
        x0 = best['vec_param']
        s0 = best['vec_struct']

        def get_diff_vec(ind, dh):
            res = np.copy(x0)
            res[ind] += dh
            return res

        ph = [get_diff_vec(i, self.h) for i in range(len(x0))]
        mh = [get_diff_vec(i, -self.h) for i in range(len(x0))]

        cs = ph + mh
        ss = [s0] * len(cs)
        dop_infos = [{'title': '4grad', 'id_x0': best_id}] * len(cs)
        return op.add_some(cs, ss, dop_infos)

    def step_up(self, pop):
        pass
