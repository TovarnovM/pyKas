import numpy as np


class Powder(object):
    """
    Описание пороха
    """
    def __init__(self, param_powder):
        self.I_k = param_powder['I_k']
        self.T_1 = param_powder['T_1']
        self.z_k = param_powder['Z_k']
        self.alpha_k = param_powder['alpha_k']
        self.etta = param_powder['etta']
        self.f = param_powder['f']
        self.k_1 = param_powder['k_1']
        self.k_2 = param_powder['k_2']
        self.k_f = param_powder['k_f']
        self.k_l = param_powder['k_l']
        self.lambda_1 = param_powder['lambda_1']
        self.lambda_2 = param_powder['lambda_2']
        self.name = param_powder['name']
        self.ro = param_powder['ro']

        self.f = self.f * 1e6
        self.alpha_k = self.alpha_k / 1e3
        self.ro = self.ro * 1e3
        self.I_k = self.I_k * 1e6

    def psi(self, z):
        """
        Функция газоприхода (относительный объем сгоревшего свода)
        :param z: относительная толщина сгоревщего свода
        :return: функция газоприхода
        """
        cond1 = z < 0
        cond2 = np.logical_and(z >= 0, z < 1)
        cond3 = np.logical_and(1 <= z, z < self.z_k)
        cond4 = z >= self.z_k

        psi = np.zeros_like(z)

        psi[cond1] = 0
        psi[cond2] = self.k_1 * z[cond2] * (1 + self.lambda_1 * z[cond2])
        psi[cond3] = self.k_1 * (1 + self.lambda_1) + self.k_2 * (z[cond3] - 1) * (1 + self.lambda_2 * (z[cond3] - 1))
        psi[cond4] = 1

        return psi

    def dpsi_dz(self, z):
        """
        Производна psi по z
        :param z: относительная толщина сгоревщего свода
        :return: производна psi по z
        """
        cond1 = z < 0
        cond2 = np.logical_and(z >= 0, z < 1)
        cond3 = np.logical_and(1 <= z, z < self.z_k)
        cond4 = z >= self.z_k

        dpsi = np.zeros_like(z)

        dpsi[cond1] = 0
        dpsi[cond2] = self.k_1 * (1 + 2 * self.lambda_1 * z[cond2])
        dpsi[cond3] = self.k_2 * (1 + 2 * self.lambda_2 * (z[cond3] - 1))
        dpsi[cond4] = 0

        return dpsi
