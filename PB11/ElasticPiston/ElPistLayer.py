import numpy as np
import copy
import math
from Pneumatics.PnLayer import PnLayer


class ElPistLayer(PnLayer):
    def __init__(self, solver_grid, func_in, get_flux_left, get_flux_right, tube):
        super().__init__(solver_grid, func_in, get_flux_left, get_flux_right, tube)
        self.k0 = solver_grid['consts']['k0']
        self.c0 = solver_grid['consts']['c0']
        self.ro0 = solver_grid['consts']['ro0']
        self.sigmas = solver_grid['consts']['sigmas']
        self.b1 = solver_grid['consts']['b1']
        self.b2 = solver_grid['consts']['b2']
        self.mu = solver_grid['consts']['mu']
        self.taus = solver_grid['consts']['taus']

        self.e = self.get_energ(self.p, self.ro)
        self.q = self.init_arr_q()
        self.h = self.get_arr_h()

    def get_arr_h(self):
        """
        Получение вектора h
        :return: список нумпи массивов h1, h2, h3
        """
        u = self.q[1] / self.q[0]
        u_right = np.append(0.5 * (u + np.roll(u, -1))[:-1], self.V[-1])
        u_left = np.roll(u_right, 1)
        u_left[0] = self.V[0]

        ha = 2 * (u_right - u_left) / self.dx - u * self.tube.get_stuff(self.x) / self.tube.get_S(self.x_c)
        hh = (0.5 / math.sqrt(3)) * np.abs(ha)

        tauxx = np.zeros_like(ha, dtype=np.float64)

        cond1 = np.abs(hh) > 1e-11
        tauxx[cond1] = (2 / 3) * (self.mu + 0.5 * self.taus / hh[cond1]) * ha[cond1]

        qw = np.zeros_like(tauxx, dtype=np.float64)

        sigmaxx = -self.p + tauxx
        sigmannw = sigmaxx - 1.5 * tauxx

        sigmantw = np.zeros_like(sigmaxx)
        cond2 = self.p < self.sigmas

        sigmantw[cond2] = -self.p[cond2] * self.k0 * (1 + self.b1 * u[cond2]) * np.exp(-self.b2 * u[cond2]) * \
                          np.sign(self.q[1][cond2])

        sigmantw[cond2 == False] = -self.sigmas * self.k0 * (1 + self.b1 * u[cond2 == False]) * \
                                   np.exp(-self.b2 * u[cond2 == False]) * np.sign(self.q[1][cond2 == False])

        h1 = np.zeros_like(sigmantw)
        h2 = -self.tube.get_stuff(self.x) * sigmannw + math.pi * self.tube.d(self.x_c) * sigmantw
        h3 = math.pi * self.tube.d(self.x_c) * (sigmantw * u - qw)

        return [h1, h2, h3]

    def get_energ(self, p, ro):
        """
        Получение энергии
        :param p: давление
        :param ro: плотность
        :return: энергия
        """
        return (p - self.c0 * self.c0 * (ro - self.ro0)) / (self.const.g[9] * ro)

    def get_pressure(self, q):
        """
        Получение давления
        :param q:
        :return: давление
        """
        p = self.const.g[9] * q[0] * (q[2] / q[0] - 0.5 * np.square(q[1] / q[0])) \
                  + self.c0 * self.c0 * (q[0] - self.ro0)
        # p = self.const.g[9] * q[0] * self.get_energ(self.p, q[0]) + self.c0 * self.c0 * (q[0] - self.ro0)
        cond = p < 0
        # cond = q[0] <= self.ro0
        p[cond] = 101325
        return p

    def get_Csound(self, ro, p):
        """
        Получение скорости звука
        :param ro: плотность
        :param p: давление
        :return: скорость звука
        """
        return np.sqrt((self.const.g[0] * p + self.ro0 * self.c0 * self.c0) / ro)

    def clone(self, l):
        """
        Копирование слоя
        :param l: слой
        :return: скопированный слой
        """
        l1 = ElPistLayer(self.solver_grid, self.func_in, self.flux_left, self.flux_right, self.tube)
        l1.q = [np.copy(ar) for ar in l.q]
        l1.h = [np.copy(ar) for ar in l.h]
        l1.x = np.copy(l.x)
        l1.V = np.copy(l.V)
        l1.W = l.W
        l1.S = l.S
        l1.ds = l.ds
        l1.p = np.copy(l.p)
        l1.time = copy.copy(l.time)
        return l1

    def get_a(self, j, p):
        """
        Получение ускорения поршня
        :return: ускорение поршня
        """
        u = self.q[1][j] / self.q[0][j]
        if j == -1:
            ha = 4 * (self.V[-1] - u) / self.dx[-1]
            hh = (0.5 / math.sqrt(3)) * abs(ha)
            if hh != 0:
                tauxx = (2 / 3) * (self.mu + 0.5 * self.taus / hh) * ha
            else:
                tauxx = 0
            a = (self.p[-1] - p) * self.S[-1] / self.m2
            return a
        elif j == 0:
            ha = 4 * (u - self.V[0]) / self.dx[0]
            hh = (0.5 / math.sqrt(3)) * abs(ha)
            if hh != 0:
                tauxx = (2 / 3) * (self.mu + 0.5 * self.taus / hh) * ha
            else:
                tauxx = 0
            a = (p - self.p[0]) * self.S[0] / self.m1
            return a
