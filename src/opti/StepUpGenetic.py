import random as rnd
import numpy as np
from functional import seq
from . import OptiPerson


class StepUpGenetic(object):
    def __init__(self, chr_contr):
        self.chr_contr = chr_contr
        self.prob_cross = 0.8
        self.prob_mut = 0.3
        self.prob_mut_gene = 0.3
        self.pop_count = 10
        self.elite_count = 2
        self.select_parents = self.select_parents_roulette
        self.select_mutants = self.select_mutants_uniform
        self.select_survivers = self.select_survivers_roulette

    def make_children(self, parents):
        """
        parents - [(OptiPerson1,OptiPerson1), ...]
        returns [OptiPerson,OptiPerson, ...]
        """

        def create_one(tup):
            op1, op2 = tup
            chrom_child = self.chr_contr.cross(op1.get_best_chromo(), op2.get_best_chromo())
            return OptiPerson.OptiPerson(self.chr_contr, chrom_child)

        return [create_one(t) for t in parents]

    def make_mutants(self, candidates):
        """
        candidates - [OptiPerson, ...]
        return - [OptiPerson, ...]
        """
        gene_keys = list(self.chr_contr.ranges_dict.keys())

        def mutate(c):
            gene_mut_count = 0
            while gene_mut_count == 0:
                gene_mut_count = np.random.binomial(len(gene_keys), self.prob_mut_gene)
            gene_mutate = np.random.choice(gene_keys, gene_mut_count, replace=False)
            mutant_chromo = self.chr_contr.mutate(c.get_best_chromo(), gene_mutate)
            return OptiPerson.OptiPerson(self.chr_contr, mutant_chromo)

        return [mutate(c) for c in candidates]

    def select_parents_roulette(self, pop):
        """
        pop - словарь "имя - OptiPerson"
        returns [(OptiPerson1,OptiPerson2), ...]
        """
        slst = seq(pop.values()).map(lambda op: op.get_best()).sorted(lambda wr_chr: wr_chr[OptiPerson.fit_key]).to_list()
        fitnss = seq(slst).map(lambda c: c[OptiPerson.fit_key]).to_list()
        delt_fit = fitnss[-1] - fitnss[0]
        p = seq(fitnss).map(lambda f: (f - fitnss[0] + 0.1 * delt_fit) / 1.1 / delt_fit).to_list()
        summ = seq(p).sum()
        p = seq(p).map(lambda f: f / summ).to_list()
        n_pairs = np.random.binomial(len(slst), self.prob_cross)

        def get_by_inds(nar):
            name1 = slst[nar[0]]['name']
            name2 = slst[nar[1]]['name']
            if name1 in pop.keys() and name2 in pop.keys():
                return pop[name1], pop[name2]
            return None

        res = seq.range(n_pairs) \
            .map(lambda i: np.random.choice(len(slst), 2, replace=False, p=p)) \
            .map(get_by_inds) \
            .filter(lambda par: par) \
            .to_list()
        return res

    def select_mutants_uniform(self, pop):
        """
        pop - словарь "имя - OptiPerson"
        return [OptiPerson, ...]
        """
        lst = list(pop.values())
        mut_count = np.random.binomial(len(lst), self.prob_mut)
        return list(np.random.choice(lst, mut_count, replace=False))

    def select_survivers_roulette(self, pop, children, mutants):
        """
        pop - словарь "имя - OptiPerson"
        children - [OptiPerson, ...]
        mutants - [OptiPerson, ...]
        returns [OptiPerson,OptiPerson, ...]
        """
        slst = seq(pop.values()).sorted(lambda op: op.get_best()[OptiPerson.fit_key], reverse=True).to_list()
        elita = [e.copy() for e in slst[:self.elite_count]]
        loosers = [l.copy() for l in slst[self.elite_count:]]
        loosers_p = [0.5] * len(loosers)
        chindren_p = [1] * len(children)
        mutants_p = [1.5] * len(mutants)
        noobies_all = loosers + children + mutants
        noobies_p = loosers_p + chindren_p + mutants_p
        sum_p = sum(noobies_p)
        noobies_p = seq(noobies_p).map(lambda p: p / sum_p).to_list()
        if self.pop_count - self.elite_count >= len(noobies_all):
            return elita+ loosers+ children + mutants
        noobies = list(np.random.choice(noobies_all, self.pop_count - self.elite_count, replace=False, p=noobies_p))
        return elita + noobies

    def step_up(self, pop, seed=None):
        """
        pop - словарь "имя - OptiPerson"
        """
        if seed:
            rnd.seed(seed)
            np.random.seed(seed)
        parents = self.select_parents(pop)
        children = self.make_children(parents)
        mut_cand = self.select_mutants(pop)
        mutants = self.make_mutants(mut_cand)
        new_gener = self.select_survivers(pop, children, mutants)
        return OptiPerson.OptiPerson.lst_to_dict(new_gener)
