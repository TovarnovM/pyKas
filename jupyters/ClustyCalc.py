from rx.concurrency import CurrentThreadScheduler, ThreadPoolScheduler
from rx.subjects import Subject
from rx import Observable


class ClustyCalc(object):
    def __init__(self, fun, n_workers):
        """

        :param pauser:
        :param fun: wchr -> lbv.apply_async()
        """
        self.fun = fun
        self.curr_tr_sched = CurrentThreadScheduler()
        self.thr_pool_sched = ThreadPoolScheduler()

        self.sub_new_jobs = Subject()
        self.sub_done_jobs = Subject()
        self.sub_freeworkers = Subject()
        self.n_workers = n_workers

        self.sub_new_jobs \
            .map(self.on_next_njob_map) \
            .flat_map(lambda ar: Observable.just(ar).observe_on(self.thr_pool_sched).map(lambda ar2: ar2.result())) \
            .observe_on(self.curr_tr_sched) \
            .subscribe(on_next=self.job_done)

    def on_next_njob_map(self, wchr):
        res = self.fun(wchr)
        self.n_workers -= 1
        self.sub_freeworkers.on_next(self.n_workers)
        return res

    def job_done(self, res):
        self.n_workers += 1
        self.sub_freeworkers.on_next(self.n_workers)
        self.sub_done_jobs.on_next(res)
