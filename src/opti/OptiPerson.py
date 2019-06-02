from . import NameGenerator
from functional import seq
import json
import numpy as np
import copy


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


name_generator = NameGenerator.NameGenerator.get_standart()
fit_key = 'fitness'


class OptiPerson(object):
    """
    Класс для хранения "истории" передвижения особи в параметрическом пространстве
    
    Поля класса: 
        name       - имя особи
        chr_contr  - ссылка на ChromoController, в котором хранится структура хромосомы особи
        idcount    - счетчик id для записей в истории
        history    - словарь, хранящий "истории" передвижения особи в параметрическом пространстве, т.е. хромосомы 
                     словарь имеет следующие ключи-значения: 
                         [fit_key](см. выше) - хранит посчитанную фитнесс-функцию особи (чем больше, тем лучше)
                         ['id']      - хранит id записи 
                         ['dop_info'] - хранит дополнительную информацию о записи (например номер поколения/эпохи)
                         ['name']    - хранит имя особи (например для поиска особи при расчете фитнесс-функции)
                         ['chromo']  - хранит собственно хромосому
                         ['vec_param'] - хранит numpy вектор с нормализованным представлением хромосомы (для normalizable генов)
                         ['vec_struct'] - хранит numpy вектор с нормализованным представлением хромосомы (для normalizable ==False генов)
    """

    def __init__(self, chr_contr, chromo0=None, fitness=None, name=None):
        """
        Конструктор класса для хранения "истории" передвижения особи в параметрическом пространстве

        Параметры: 
            name       - имя особи
            chr_contr  - ссылка на ChromoController, в котором хранится структура хромосомы особи
        """
        if not name:
            name = name_generator.new()
        self.name = name
        self.chr_contr = chr_contr
        self.idcount = 0
        if not chromo0:
            chromo0 = chr_contr.get_chromo()
        self.history = {}
        self.add_new_chromo(chromo0, fitness)

    @classmethod
    def from_vecs(cls, chr_contr, vec_param, vec_struct, fitness=None, name=None):
        c = chr_contr.get_chromo_from_vecs(vec_param, vec_struct)
        return cls(chr_contr, c, fitness, name)

    def copy(self):
        res = OptiPerson(self.chr_contr, name=self.name)
        res.history = copy.deepcopy(self.history)
        res.idcount = self.idcount
        return res

    def change_name(self, name=None):
        if not name:
            name = name_generator.new()
        self.name = name
        for h in self.history.values():
            h['name'] = name
        return name

    def get_wchomo(self, chromo, fitness=None, dop_info=None):
        vec_param = self.chr_contr.get_vec_param(chromo)
        vec_struct = self.chr_contr.get_vec_struct(chromo)
        wr_chr = {fit_key: fitness,
                  'name': self.name,
                  'chromo': chromo,
                  'id': self.idcount,
                  'dop_info': dop_info,
                  'vec_param': vec_param,
                  'vec_struct': vec_struct
                  }
        return wr_chr

    def add_new_chromo(self, chromo, fitness=None, dop_info=None):
        """
        Метод добавления новой хромосомы ("обернутой функцией wrap_chromo")

        Параметры: 
            chromo  - новая хромосома
            dop_info - хранит дополнительную информацию о записи (например номер поколения/эпохи)
        """
        wr_chr = self.get_wchomo(chromo, fitness, dop_info)
        self.idcount += 1
        self.history[wr_chr['id']] = wr_chr
        return wr_chr['id']

    def get_fitlessness(self):
        """
        Метод получения списка записей истории, у которых еще не посчитана фитнесс-функция
        """
        return seq(self.history.values()).where(lambda h: not h[fit_key]).to_list()

    def __getitem__(self, key):
        return self.history[key]

    def __iter__(self):
        return iter(self.history)

    def get_best(self, fit_too=False):
        """
        Метод получения лучшего решения за всю историю)
        """
        try:
            res = (seq(self.history.values())
                   .where(lambda h: h[fit_key])
                   .max_by(lambda h: h[fit_key]))
            return (res, res[fit_key]) if fit_too else res
        except ValueError:
            res = seq(self.history.values()).last()
            return (res, None) if fit_too else res

    def get_best_chromo(self):
        bsol = self.get_best()
        return bsol['chromo']

    @property
    def best_fitness(self):
        bf = self.get_best()
        return bf[fit_key]

    def remove(self, key):
        """
        Метод удаления записи
        """
        if isinstance(key, list):
            for k in key:
                self.history.pop(k, None)
        else:
            self.history.pop(key, None)

    def remove_exept(self, keys):
        nhis = {k: self[k] for k in keys}
        self.history = nhis

    def remove_exept_best(self):
        best = self.get_best()['id']
        self.remove_exept([best])

    @property
    def get_last_with_fit(self):
        """
        Метод получения последней посчитанной записи (считается по id, чем id больше, тем запись новее)
        """
        id = seq(self.history).where(lambda id: self.history[id][fit_key]).order_by(lambda id: id).last()
        return self[id]

    @property
    def vel(self):
        if len(self.history) < 2:
            last = self.get_last()
            return np.random.uniform(-1, 1, len(last['vec_param']))
        id_last, id_pre_last = seq(self.history.keys()).sorted(reverse=True).take(2).to_list()
        return self[id_last]['vec_param'] - self[id_pre_last]['vec_param']

    def get_last(self):
        id = seq(self.history).max_by(lambda k: k)
        return self[id]

    def to_str(self):
        return json.dumps(self.history, cls=NumpyEncoder)

    def __str__(self):
        return self.name + ' ' + self.to_str()

    def __repr__(self):
        return f"OptiPerson(name={self.name}, best={self.best_fitness})"

    def step_vel(self, vel):
        last = self.get_last()
        chr_new = self.chr_contr.chromo_vel(last, vel)
        self.add_new_chromo(chr_new)

    def clear_swarm_hist(self, except_ids=[]):
        if len(self.history) < 2:
            return
        id_last, id_pre_last = seq(self.history.keys()).sorted(reverse=True).take(2).to_list()
        id_best = self.get_best()['id']
        surv = except_ids + [id_last, id_pre_last, id_best]
        self.history = {id: self.history[id] for id in set(surv)}

    @classmethod
    def from_string(cls, string, chr_contr):
        hist_str_id = json.loads(string)
        hist = {int(k): hist_str_id[k] for k in hist_str_id}
        for r in hist.values():
            r['vec_param'] = np.array(r['vec_param'])
            r['vec_struct'] = np.array(r['vec_struct'])
        name = next(iter(hist.values()))['name']
        idmax = seq(hist).max()
        res = cls(chr_contr, name=name)
        res.history = hist
        res.idcount = idmax + 1
        return res

    @classmethod
    def lst_to_dict(cls, lst):
        names = set()
        for op in lst:
            if not (op.name in names):
                names.update([op.name])
                continue
            name = name_generator.new()
            while name in names:
                name = name_generator.new()
            op.change_name(name)
            names.update([op.name])
        return {op.name: op for op in lst}

    @classmethod
    def get_fitnessless_chromo_list(cls, pop):
        """
        Возвращает список хромосом для расчета
        :param pop:  словарь "имя - OptiPerson"
        :return: [ {fit_key: None, 'name': name1, 'chromo': chromo1, 'id': ..., 'dop_info': dop_info},
                   {fit_key: None, 'name': name1, 'chromo': chromo1, 'id': ..., 'dop_info': dop_info}, ...]
        """
        return seq(pop.values()).flat_map(lambda op: op.get_fitlessness()).to_list()

    @classmethod
    def init_fitnesses(cls, pop, results):
        """
        Функция инициализации (фиксации) фитнесс-функций поколения pop уже посчитанными результатами results
        :param pop: словарь "имя - OptiPerson"
        :param results: список кортежей ({fit_key: None, 'name': name1, 'chromo': chromo1, 'id': ..., 'dop_info': dop_info}, calc_fitness)
        :return: None
        """
        for c, fit in results:
            name = c['name']
            cid = c['id']
            pop[name].history[cid][fit_key] = fit



    def add_some(self, vecs_param, vecs_struct, dop_infos):
        return [self.add_new_chromo(self.chr_contr.get_chromo_from_vecs(vp, vs),
                                    dop_info=di) for vp, vs, di in zip(vecs_param, vecs_struct, dop_infos)]

    def prep_chromos_4line(self, id_x0, grad, n, length, h_min=1E-5):
        vecs, step = self.get_line_vecs(id_x0, grad, n, length, h_min)
        s0 = self[id_x0]['vec_struct']
        return [self.add_new_chromo(self.chr_contr.get_chromo_from_vecs(v, s0),
                                    dop_info={'title': '4line', 'id_x0': id_x0}) for v in vecs]

    def __cmp__(self, other):
        try:
            f_me, f_other = self.best_fitness
            if f_me < f_other:
                return -1
            if f_me > f_other:
                return 1
            return 0
        except Exception as e:
            return 0

    def __lt__(self, other):
        return self.best_fitness < other.best_fitness

#
# def main():
#     chr_contr = get_test_chr_contr()
#
#     op = OptiPerson(chr_contr)
#     for i in range(10):
#         op.add_new_chromo(chr_contr.get_chromo())
#     op[0]['fitness'] = 11
#     op[1]['fitness'] = 10
#     op[3]['fitness'] = 9
#
#     print('Все записи в истории')
#     [print('id = ', h, '   ', op[h]) for h in op]
#     print('лушчее решение')
#     print(op.get_best())
#     print('последнее решение')
#     print(op.get_last_with_fit)
#
#     string = op.to_str()
#     op2 = OptiPerson.from_string(string, chr_contr)
#     [print('id = ', h, '   ', op2[h]) for h in op2]
#

# if __name__ == '__main__':
#     main()
