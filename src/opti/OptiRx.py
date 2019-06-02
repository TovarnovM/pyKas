from functional import seq
from rx.concurrency import ThreadPoolScheduler
from rx.subjects import Subject, ReplaySubject
from rx import Observable

from opti.Generation import Generation
from opti.OptiSLSQP import OptiSLSQP
from opti.StepUpGenetic import StepUpGenetic
from opti.StepUpGrad import reduce_pop
from .Pauser import Pauser

import threading

default_options = {
    'pop_count': 100,
    'elite_count': 5,
    'prob_cross': 0.7,
    'prob_mut': 0.2,
    'prob_mut_gene': 0.5,
    'seed0': None,
    'n_gener_switch': 10,
    'n_gener_done': 20,
    'n_op_direct': 5,
    'n_gener_max': 300
}


class OptiRx(object):
    def __init__(self, chr_contr, executor, gen0=None, opts=default_options):
        self.gen0 = gen0
        self.chr_contr = chr_contr
        self.pauser = Pauser()
        self.executor = executor
        self.executor.kpp_fun = self.pauser.kpp
        self.sub_chr_best = ReplaySubject()
        self.sub_gen_solved = Subject()
        self.sub_gen_unsolved = Subject()
        self.subj_done = Subject()

        self.su_genetic = StepUpGenetic(chr_contr)

        self.opts = dict(default_options)
        self.init_opts(opts)
        self.stop_flag = True
        self.gen_flag = True
        self.init_logic()

    def init_opts(self, opts):
        for k in opts:
            self.opts[k] = opts[k]
        self.su_genetic.pop_count = self.opts['pop_count']
        self.su_genetic.prob_mut = self.opts['prob_mut']
        self.su_genetic.prob_mut_gene = self.opts['prob_mut_gene']
        self.su_genetic.prob_cross = self.opts['prob_cross']
        self.su_genetic.elite_count = self.opts['elite_count']

    def init_logic(self):
        # (fit, op)...
        def best_fun(gener):
            fit, op = gener.get_best(True)
            return gener.num_g, fit, op

        best_stream = self.sub_gen_solved.map(best_fun).publish()
        improve_best = best_stream.distinct(lambda tp: tp[1])
        improve_best.subscribe(self.sub_chr_best)

        def wind_switch(beep_inds):
            beep_set = set(beep_inds)

            def inner(s):
                inds = Observable.interval(0)
                return s.zip(inds, lambda x, ind: ind) \
                    .filter(lambda ind: ind in beep_set)

            return inner

        switch = best_stream\
            .window(improve_best.skip(1))\
            .flat_map(wind_switch([self.opts['n_gener_switch'], self.opts['n_gener_done']-2]))

        done = best_stream\
            .window(improve_best.skip(1))\
            .flat_map(wind_switch([self.opts['n_gener_done']]))

        def stop(ind):
            self.stop_flag = True

        def switch_f(ind):
            self.gen_flag = False

        switch.subscribe(switch_f)
        done.subscribe(stop)
        best_stream.skip(self.opts['n_gener_max']-1).subscribe(stop)

        best_stream.connect()

    def paused(self):
        return self.pauser.paused

    def pause(self):
        self.pauser.pause()

    def unpause(self):
        self.pauser.play()

    def run_sync(self):
        if self.gen0:
            gen0 = self.gen0
        else:
            gen0 = Generation(self.chr_contr, 0)
            gen0.get_init_pop(self.opts['pop_count'], self.opts['seed0'])
        self.stop_flag = False
        self.gen_flag = False
        while not self.stop_flag:
            self.sub_gen_unsolved.on_next(gen0)
            if self.gen_flag:
                gen0 = self.stepup_gen(gen0)
            else:
                gen0 = self.stepup_slsqp(gen0)
                self.gen_flag = True
        self.subj_done.on_completed()

    def run(self):
        cs1 = threading.Thread(name='calc_thread', target=self.run_sync)
        cs1.start()

    def stepup_gen(self, gen, seed=None):
        self.calc_gen(gen)
        self.sub_gen_solved.on_next(gen)
        gen_dict = self.su_genetic.step_up(gen.pop, seed)
        gen2 = Generation(self.chr_contr, gen.num_g + 1, gen_dict)
        return gen2

    def calc_gen(self, gen):
        fitnessless = gen.get_fitlessness()
        if len(fitnessless) == 0:
            return
        tick_tack = ReplaySubject()
        fset = seq(fitnessless) \
            .map(lambda wchr: (wchr['name'], wchr['id'])) \
            .to_set()
        results = []

        def remove_fset(result):
            wchr, fit = result
            tp = wchr['name'], wchr['id']
            if tp in fset:
                fset.discard(tp)
                results.append(result)
            if len(fset) == 0:
                tick_tack.on_next(0)

        s1 = self.executor.sub_compl.subscribe(on_next=remove_fset)

        for fn in fitnessless:
            self.executor.sub_in.on_next(fn)
        if len(fset) == 0:
            tick_tack.on_next(0)
        tick_tack.to_blocking().first()

        tick_tack.dispose()
        s1.dispose()

        gen.init_fitnesses(results)
        if len(gen.get_fitlessness()) > 0:
            raise AssertionError("Посчитались не все гены в хромосоме")

    def stepup_slsqp(self, gen0):
        self.calc_gen(gen0)
        lst = reduce_pop(gen0.pop_list, self.opts['n_op_direct'])
        opti_lst = [OptiSLSQP(op, self.executor.sub_in, self.executor.sub_compl) for op in lst]
        tp = ThreadPoolScheduler()
        calc_stream = Observable \
            .from_(opti_lst) \
            .flat_map(lambda op: Observable.just(op).observe_on(tp).map(lambda op: op.run()))

        calc_stream.to_blocking().last_or_default(0)
        for op in lst:
            op.remove_exept_best()

        for o_sls in opti_lst:
            o_sls.dispose_conn()
        gen1 = Generation(self.chr_contr, gen0.num_g + 1, lst)
        self.sub_gen_solved.on_next(gen1)
        gen1.fill_pop(self.opts['pop_count'])
        return gen1
