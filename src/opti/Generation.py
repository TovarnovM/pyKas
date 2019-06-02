from functional import seq
from . import OptiPerson
from . import NameGenerator
import numpy as np
import random as rnd


def op_lst_to_dict(lst):
    names = set()
    for op in lst:
        if not (op.name in names):
            names.update([op.name])
            continue
        name = NameGenerator.name_generator.new()
        while name in names:
            name = NameGenerator.name_generator.new()
        op.change_name(name)
        names.update([op.name])
    return {op.name: op for op in lst}


class Generation(object):
    def __init__(self, chr_contr, num_g, init_pop=None):
        self.chr_contr = chr_contr
        if init_pop:
            if isinstance(init_pop, list):
                self.pop = op_lst_to_dict(init_pop)
            if isinstance(init_pop, dict):
                self.pop = init_pop
        else:
            self.pop = {}
        self.num_g = num_g

    def __len__(self):
        return len(self.pop)

    def __str__(self):
        return f"Generation {self.num_g}"

    @property
    def pop_list(self):
        return list(self.pop.values())

    def get_best(self, fit_too=False):
        fit, op = seq(self.pop.values()).map(lambda op: (op.best_fitness, op)).max_by(lambda tp: tp[0])
        return (fit, op) if fit_too else op

    def save_db(self, opti_saver):
        opti_saver.delete_pops(self.num_g)
        save_list = [op.to_str() for op in self.pop.values()]
        opti_saver.write_pops(self.num_g, save_list)

    def load_db(self, opti_saver):
        lst = opti_saver.load_generation(self.num_g)
        self.pop = op_lst_to_dict([OptiPerson.from_string(s) for s in lst])

    def load_slst(self, slst):
        self.pop = op_lst_to_dict([OptiPerson.from_string(s) for s in slst])

    def get_init_pop(self, n, seed=None, elite_lst=[]):
        if seed:
            np.random.seed(seed)
            rnd.seed(seed)
        self.pop = op_lst_to_dict([OptiPerson.OptiPerson(self.chr_contr) for _ in range(n - len(elite_lst))]+elite_lst)

    def get_fitlessness(self):
        """
        Метод возвращает список всех хромосом (решений) в поколении, у которых нету фитнесс-функций
        :return:
        """
        return seq(self.pop.values()).flat_map(lambda op: op.get_fitlessness()).to_list()

    def init_fitnesses(self, results):
        """
        :param results: список кортежей ({fit_key: None, 'name': name1, 'chromo': chromo1, 'id': ..., 'dop_info': dop_info}, calc_fitness)
        :return:
        """
        for c, fit in results:
            name = c['name']
            cid = c['id']
            self.pop[name].history[cid]['fitness'] = fit

    def fill_pop(self, n, seed=None):
        if n > len(self.pop):
            self.get_init_pop(n, seed, self.pop_list)
        else:
            self.get_init_pop(n, seed, self.pop_list[:n])




