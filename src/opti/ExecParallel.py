from rx import Observable
from rx.concurrency import ThreadPoolScheduler, CurrentThreadScheduler
from enum import Enum

ExState = Enum('ExState', 'paused playing error')


class ExecParallel(object):
    """
    класс-обертка для вычислительного кластера
    Поля класса:
        client     - ссылка на ipyparallel.Client()
        lbv        - ссылка на client.load_balanced_view()
        fit_fun    - ссылка на фитнесс-функцию. Функция принимает 1 аргумент (хромосому)
        async_res  - ссылка на async result от lbv.map_async(...)
        stream     - Observable поток кортежей результатов вычисления (патаметр, результат вычисления)
    """

    def __init__(self, client, fit_fun):
        """
        Конструктор -__- класса ExecParallel
        
        Аргумены:
            client     - ссылка на ipyparallel.Client()
            fit_fun    - ссылка на фитнесс-функцию. Функция принимает 1 аргумент (хромосому), возвращать она должна
                         кортеж (хромосома, фитнесс-функция)
        """
        try:
            self.client = client
            self.lbv = client.load_balanced_view()
            self.fit_fun = fit_fun
            self._state = ExState.paused
            self.scheduler = ThreadPoolScheduler()
            self.async_res = None
            self.stream = None
        except:
            self._state = ExState.error

    def state(self):
        """
        Состояние, в котором находится объект:
            ExState.paused 
            ExState.playing 
            ExState.error   
        """
        return self._state

    def execute(self, params):
        """
        основная функция. Функция высшего порядка map для функции fit_fun. Для вычисления используется кластер client/lbv.
        Функция возвращает Observable-поток, на который можно подписаться) Поток состоит из кортежа
        (параметр, fit_fun(параметр))
        
        Аргументы:
            params - список/iterable объект, содержащий элементы-аргументы для fit_fun 
        """

        self.async_res = self.lbv.map_async(self.fit_fun, params, ordered=False)

        def callb():
            self._state = ExState.paused

        self.stream = Observable.from_(iter(self.async_res), scheduler=self.scheduler) \
            .publish().auto_connect(2)
        self.stream.subscribe(on_completed=callb)
        self._state = ExState.playing
        return self.stream

    def wait(self, timeout=-1):
        """
        Wait until the result is available or until `timeout` miliseconds pass.
        
        Возвращает 'beep beep', если сработал таймер (таймер истек раньше, чем закончились расчеты),
        или 0, если расчеты закончились раньше
        """
        if self.stream is None:
            return 0

        if timeout < 0:
            lasty = self.stream.to_blocking().last_or_default(0)
            return 0

        timer_message = 'beep beep'
        timer = Observable.timer(timeout).map(lambda x: timer_message)
        lasty = self.stream.last_or_default(0) \
            .amb(timer) \
            .observe_on(CurrentThreadScheduler()) \
            .to_blocking() \
            .first()
        return timer_message if lasty == timer_message else 0

    def abort(self):
        """
        Отменяет расчет,
        врзвращает количество отмененных задач
        """
        try:
            if self.async_res is None:
                return 0
            if self.stream is None:
                return 0
            # if self.async_res.
            results = []
            self.stream \
                .observe_on(CurrentThreadScheduler()) \
                .filter(lambda x: not isinstance(x, tuple)) \
                .subscribe(on_next=lambda x: results.append(x))

            self.async_res.abort()

            self.wait()
            return len(results)
        except:
            return 0


if __name__ == '__main__':
    from ipyparallel import Client

    cl = Client()
    lbv = cl.load_balanced_view()

    with cl[:].sync_imports():
        import time
        import random


    def param_res_wrapp(f):
        """
        Функция для обертки обычной фитнесс-функции,
        возвращает кортеж (патаметр, f(параметр))
        """

        def inner(par):
            return par, f(par)

        return inner

    def tst_fun(par):
        tm = random.random() * 3 + 1
        time.sleep(tm)
        return f"result of {par}, time {tm}"

    print('запущено ', len(cl[:]), ' ядер')
    ex = ExecParallel(cl, param_res_wrapp(tst_fun))
    n = 33
    obs = ex.execute(range(n))
    results = []

    def onnext(x):
        print(x)
        results.append(x)

    obs.filter(lambda x: isinstance(x, tuple)) \
        .subscribe(on_next=onnext,
                   on_completed=lambda: print('done'),
                   on_error=lambda e: print('error ', e))

    print('запущено ', n, ' задач')
    t = 2000
    res = ex.wait(t)
    print('отменено спустя', t, ' мс, код = ', res)
    aborted = ex.abort()
    print('отменено (не посчитано) ', aborted, ' шт.')
    print('посчитано               ', len(results), ' шт.')
    print('всего                   ', aborted + len(results), ' шт.')




