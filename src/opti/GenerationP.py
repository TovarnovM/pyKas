from .OptiPerson import OptiPerson
import random as rnd
import numpy as np
import bisect
from rx.subjects import Subject



class GenerationP(object):
    def __init__(self, chr_contr, n_max):
        self.chr_contr = chr_contr
        self.n_max = n_max
        self.bests = []
        self.prob_cross = 0.7
        self.prob_mut = 0.3
        self.prob_mut_gene = 0.3
        self.fit_threshold = 100
        self.candi_prob = 0.4
        self.sub_wchr4calc = Subject()
        self.sub_reqest_4newwchr = Subject()
        self.sub_take_calced = Subject()
        self.sub_update_pop = Subject()
        self.sub_update_best_op = Subject()

        self.sub_reqest_4newwchr\
            .map(lambda t: self.get_new_wcr())\
            .subscribe(self.sub_wchr4calc)

        self.sub_take_calced.subscribe(on_next=self.on_next4sub_take_calced)


    def on_next4sub_take_calced(self, result):
        wchr, fit = result
        op = OptiPerson(self.chr_contr, chromo0=wchr['chromo'], fitness=fit,  name=wchr['name'])
        if self.take_calc_op(op, fit):
            self.sub_update_pop.on_next((fit, op))
            if fit == self.bests[-1][0]:
                self.sub_update_best_op.on_next((fit, op))



    def get_new_wcr(self):
        if len(self.bests) < self.n_max:
            op = OptiPerson(self.chr_contr)
            return op.get_last()

        if rnd.random() < self.prob_cross:
            op1, op2 = self.select_2par()
            child = self.make_child(op1, op2)
            return child.get_last()

        if rnd.random() < self.prob_mut:
            fit, op = rnd.choice(self.bests)
            mut = self.make_mutant(op)
            return mut.get_last()

        op = OptiPerson(self.chr_contr)
        return op.get_last()



    def take_calc_op(self, op, fit):
        """

        :param op:
        :param fit:
        :return: True - if update bests
        """
        if fit < self.fit_threshold:
            return False
        if len(self.bests) == 0:
            self.bests.append((fit, op))
            return True
        if len(self.bests) >= self.n_max and self.bests[0][0] > fit:
            return False

        bisect.insort(self.bests, (fit, op))
        if len(self.bests) > self.n_max:
            self.bests = self.bests[1:]
        return True

    def select_2par(self):
        """

        :return: (OptiPerson ,OptiPerson)
        """
        def get_winer():
            candi = [rnd.choice(self.bests)]
            while rnd.random() < self.candi_prob:
                candi.append(rnd.choice(self.bests))
            return max(candi)
        p1, p2 = get_winer(), get_winer()
        return (p1[1], p2[1]) if p1[0] != p2[0] and p1[1].name != p2[1].name else self.select_2par()


    def make_child(self, op1, op2):
        chrom_child = self.chr_contr.cross(op1.get_best_chromo(), op2.get_best_chromo())
        return OptiPerson(self.chr_contr, chrom_child)

    def make_mutant(self, op):
        """
        candidates - [OptiPerson, ...]
        return - [OptiPerson, ...]
        """
        gene_keys = list(self.chr_contr.ranges_dict.keys())
        gene_mut_count = 0
        while gene_mut_count == 0:
            gene_mut_count = np.random.binomial(len(gene_keys), self.prob_mut_gene)
        gene_mutate = np.random.choice(gene_keys, gene_mut_count, replace=False)
        mutant_chromo = self.chr_contr.mutate(op.get_best_chromo(), gene_mutate)
        return OptiPerson(self.chr_contr, mutant_chromo)

    # def set_done_op(self, ):

