from rx import Observable
from rx.concurrency import ThreadPoolScheduler
from rx.subjects import Subject
from .TestFunctions import get_funcs_diap


class ExecLocal(object):
    def __init__(self, fun, kpp_fun=None):
        self.fun = fun
        self.kpp_fun = kpp_fun
        self.sub_in = Subject()
        self.sub_compl = Subject()
        self.sub_in.subscribe(on_next=self.comp)

    def comp(self, wchr):
        if self.kpp_fun:
            self.kpp_fun()
        res = self.fun(wchr)
        self.sub_compl.on_next(res)

    @classmethod
    def get_standart(cls, fun_name):
        d = get_funcs_diap()
        fun = d[fun_name][0]

        def wfun(wchr):
            return wchr, fun(wchr['chromo'])

        return cls(wfun)

    @classmethod
    def get_standatr_names(cls):
        d = get_funcs_diap()
        return list(d.keys())


class ExecCluster(object):
    def __init__(self, fun, kpp_fun=None):
        """

        :param fun: должна возвращать AsyncResult объект
                    что-то типа:
                    def fun(wchr):
                        return lbv.apply_async(task_fun , wchr)
        """
        self.fun = fun
        self.kpp_fun = kpp_fun
        self.sub_in = Subject()
        self.sub_compl = Subject()

        tp = ThreadPoolScheduler()

        def f2(ar):
            if self.kpp_fun:
                self.kpp_fun()
            return ar.result()

        self.sub_in \
            .map(lambda wchr: self.fun(wchr)) \
            .flat_map(lambda ar: Observable.just(ar).observe_on(tp).map(f2)) \
            .subscribe(self.sub_compl)
