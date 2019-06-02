import numpy as np
from scipy.optimize import minimize
from rx.subjects import Subject, ReplaySubject, AsyncSubject
from functional import seq

class OptiSLSQP(object):
    def __init__(self, op, sub_in, sub_compl):
        self.op = op
        self.s0 = self.op.get_last()['vec_struct']
        self.name = op.name
        self.task_subject = Subject()
        self.compl_task_id_subject = Subject()
        self.compl_tasks_ids = []
        s1 = self.task_subject.subscribe(sub_in)
        self.h = 1E-7
        s2 = sub_compl.subscribe(on_next=self.gimme_res)
        self.all_subscriptions = [s1, s2]

    def dispose_conn(self):
        for s in self.all_subscriptions:
            s.dispose()

    def run(self, maxiter=1E3):
        wchr0 = self.op.get_last()
        x0 = wchr0['vec_param']
        constr = self.op.chr_contr.get_constr_4scipy(wchr0['chromo'])
        bounds = self.op.chr_contr.get_bounds_4scipy(wchr0['chromo'])
        res = minimize(self.objective_func, x0, method='SLSQP', bounds=bounds, constraints=constr, jac=self.jac_func,
                       options={'maxiter': maxiter})
        return res

    def objective_func(self, x):
        chromo = self.op.chr_contr.get_chromo_from_vecs(x, self.s0, False)
        id = self.op.add_new_chromo(chromo, dop_info={'title': '4func'})
        wchr = self.op[id]
        self.task_subject.on_next(wchr)
        self.wait_compl_ids([id])
        return -self.op[id]['fitness']

    def jac_func(self, x):
        x0 = x
        s0 = self.s0

        def get_diff_vec(ind, dh):
            res = np.copy(x0)
            res[ind] += dh
            return res, ind, int(np.sign(dh))

        ph = [get_diff_vec(i, self.h) for i in range(len(x0))]
        mh = [get_diff_vec(i, -self.h) for i in range(len(x0))]

        tps = ph + mh
        cs = seq(tps).map(lambda tp: tp[0]).to_list()
        ss = [s0] * len(tps)
        dop_infos = seq(tps)\
            .map(lambda tp: {'title': '4grad',
                             'x0': x0,
                             'ind': tp[1],
                             'sign': tp[2]})\
            .to_list()
        ids = self.op.add_some(cs, ss, dop_infos)
        for id in ids:
            self.task_subject.on_next(self.op[id])
        self.wait_compl_ids(ids)
        ids_p = ids[:len(ph)]
        ids_m = ids[len(ph):]
        res = np.zeros_like(x0)
        for i in range(len(ph)):
            wchr_p = self.op[ids_p[i]]
            wchr_m = self.op[ids_m[i]]
            f_p = -wchr_p['fitness']
            f_m = -wchr_m['fitness']
            x_p = wchr_p['vec_param'][i]
            x_m = wchr_m['vec_param'][i]
            res[i] = (f_p - f_m)/(x_p - x_m)
        return res

    def gimme_res(self, result):
        c, fit = result
        cid = c['id']
        name = c['name']
        if name != self.name:
            return
        self.compl_tasks_ids.append(cid)
        self.op.history[cid]['fitness'] = fit
        self.compl_task_id_subject.on_next(cid)

    def wait_compl_ids(self, ids):
        ids = set(ids)
        tick_tack = ReplaySubject()

        def onn(id):
            ids.discard(id)
            if len(ids) == 0:
                tick_tack.on_completed()

        s3 = self.compl_task_id_subject.subscribe(on_next=onn)
        for i in self.compl_tasks_ids:
            ids.discard(i)
        if len(ids) == 0:
            tick_tack.on_next(0)

        tick_tack.to_blocking().first_or_default(0)
        s3.dispose()
