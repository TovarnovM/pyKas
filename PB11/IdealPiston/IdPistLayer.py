import numpy as np
import copy
from Pneumatics.PnLayer import PnLayer


class IdPistLayer(PnLayer):
    def __init__(self, param_piston, gas, name, geom, get_x_v_left, get_x_v_right, get_flux_left, get_flux_right,
                 func_in, m):
        super().__init__(gas, name, geom, get_x_v_left, get_x_v_right, get_flux_left, get_flux_right, func_in, m)

        self.param_piston = param_piston
        self.B = param_piston['B']
        self.Ck = param_piston['Ck']
        self.ro0 = param_piston['ro0']

    def get_energ(self, p, ro):
        """
        Получение экнергии
        :param p: давление
        :param ro: плотнность
        :return: энергия
        """
        return (p / self.const.g[9]) * (1 / ro)

    def get_pressure(self, q):
        """
        Получение давления
        :param q:
        :return: давление
        """
        Ro = q[0] / self.ro0
        p = np.zeros(self.n, dtype=np.float64)
        cond = q[0] >= self.ro0
        p[cond] = self.B * Ro[cond] * (Ro[cond] - 1) / ((self.Ck - Ro[cond]) ** 2)
        return p

    def get_Csound(self, ro, p):
        """
        Получение скорости звука
        :param ro: плотность
        :param p: давление
        :return: скорость звука
        """
        return np.sqrt(self.B * self.ro0 * (self.ro0 - self.Ck * (2 * ro - self.ro0)) / ((ro - self.Ck * self.ro0) ** 3))

    def clone(self, l):
        """
        Копирование слоя
        :param l: слой
        :return: скопированный слой
        """
        l1 = IdPistLayer(self.param_piston, self.gas, self.name, self.geom, self.left_border, self.right_border, self.flux_left, self.flux_right, self.func_in, self.m)
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
