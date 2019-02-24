import numpy as np
import copy
from OneVelocity.OvAUSMplus import get_f
from Invariants.Constants import Constants
from Invariants.Tube import Tube
from Invariants.Powder import Powder


class OvLayer(object):
    def __init__(self, solver_grid, func_in, get_flux_left, get_flux_right, tube):
        # Время расчета
        self.time = 0
        self.solver_grid = solver_grid
        self.func_in = func_in
        # Количество ячеек
        self.n = solver_grid['n_cells']
        # Массив координат узлов, размерностью n+1
        self.x = np.linspace(solver_grid['xl'], solver_grid['xr'], self.n + 1, dtype=np.float64)
        self.x_right = np.roll(self.x, -1)
        # Шаги по пространству
        self.dx = (self.x_right - self.x)[:-1]
        # Массив координат центров ячеек, размерностью n
        self.x_c = ((self.x + self.x_right) / 2)[:-1]
        # Массив скоростей узлов, размерностью n+1
        self.V = np.linspace(solver_grid['vl'], solver_grid['vr'], self.n + 1, dtype=np.float64)
        # Константы с показателем адиабаты и коволюм
        self.const = Constants(solver_grid['consts']['gamma'])
        # Параметры трубы
        self.tube = tube
        # Параметры пороха
        self.powd = Powder(solver_grid['consts']['param_powder'])

        self.nu = solver_grid['consts']['nu']

        self.ds = self.tube.get_stuff(self.x)  # Нумпи массив dS/dx, размерностью n
        self.S = self.tube.get_S(self.x)  # Нумпи массив площадей в координатах узлов, размерностью n+1
        self.W = self.tube.get_W(self.x)  # Нумпи массив объемов, размерностью n

        self.ro = np.zeros(self.n)
        self.u = np.zeros(self.n)
        self.p = np.zeros(self.n)
        self.z = np.zeros(self.n)

        for i in range(self.n):
            self.ro[i], self.u[i], self.p[i], self.z[i] = func_in(self.x_c[i], solver_grid, 0)

        self.y = self.powd.psi(self.z)
        self.e = self.get_energ(self.p, self.ro, self.y)

        self.q = self.init_arr_q()  # Список нумпи массивов q1, q2, q3, размерностями n

        self.h = self.get_arr_h()  # Список нумпи массивов h1, h2, h3, размерностями n

        self.flux_left = get_flux_left  # Получение потока через левую границу
        self.flux_right = get_flux_right  # Получение потока через правую границу

    def init_arr_q(self):
        """
        Инициализация вектора q
        :return: список q1, q2, q3, q4
        """
        q1 = self.ro
        q2 = self.ro * self.u
        q3 = self.ro * (self.e + 0.5 * np.square(self.u))
        q4 = self.ro * self.z
        return [q1, q2, q3, q4]

    def get_arr_h(self):
        """
        Получение вектора h
        :param l: слой
        :return: список h1, h2, h3, h4
        """
        h1 = np.zeros(self.n, dtype=np.float64)
        h2 = self.p * self.ds
        h3 = np.zeros(self.n, dtype=np.float64)
        h4 = self.q[0] * self.tube.get_S(self.x_c) * (self.p ** self.nu) / self.powd.I_k
        return [h1, h2, h3, h4]

    def get_energ(self, p, ro, y):
        return (p / self.const.g[9]) * (1 / ro - (1 - y) / self.powd.ro - self.powd.alpha_k * y) + \
               (1 - y) * self.powd.f / self.const.g[9]

    def get_pressure(self, q):
        """
        Получение давления
        :param q: слой
        :return: давление
        """
        self.y = self.powd.psi(q[3] / q[0])
        return (self.const.g[9] * (q[2] / q[0] - 0.5 * np.square(q[1] / q[0])) -
                (1 - self.y) * self.powd.f) / (1 / q[0] - (1 - self.y) / self.powd.ro - self.powd.alpha_k * self.y)

    def get_Csound(self, ro, p, y):
        """
        Получение скорости звука
        :param q: список q
        :param p: давление
        :return: скорость звука
        """
        # psi = self.powd.psi(z)
        return np.sqrt(p / (self.const.g[8] * (1 / ro - (1 - y) / self.powd.ro - self.powd.alpha_k * y))) / ro

    def get_param(self, q):
        """
        Пересчет параметров газа из вектора q
        :param q: список q1, q2, q3, q4
        :return: плотность, скорость, внутреннюю энергию, относительную толщину сгоревшего свода
        """
        ro = q[0]
        u = q[1] / q[0]
        e = q[2] / q[0] - 0.5 * (u ** 2)
        z = q[3] / q[0]
        return ro, u, e, z

    def time_step(self):
        """
        Получение максимального шага по времени
        :return: шаг по времени
        """
        self.Cs = self.get_Csound(self.q[0], self.p, self.q[3] / self.q[0])
        Vmax = max(self.Cs + np.abs(self.q[1]) / self.q[0])
        Vmax = max(Vmax, max(self.V))
        dx = min([self.x[i] - self.x[i - 1] for i in range(1, len(self.x))])
        tau = dx / Vmax
        return tau

    def clone(self, l):
        """
        Копирование слоя
        :param l: слой
        :return: скопированный слой
        """
        l1 = OvLayer(self.solver_grid, self.func_in, self.flux_left, self.flux_right, self.tube)
        l1.q = [np.copy(ar) for ar in l.q]
        l1.h = [np.copy(ar) for ar in l.h]
        l1.x = np.copy(l.x)
        l1.x_right = np.copy(l.x_right)
        l1.x_c = np.copy(l.x_c)
        l1.dx = np.copy(l.dx)
        l1.V = np.copy(l.V)
        l1.W = l.W
        l1.S = l.S
        l1.ds = l.ds
        l1.p = np.copy(l.p)
        l1.time = l.time
        return l1

    def get_dQs(self):
        """
        Функция пересчета правой части диф уравнения
        :return:
        """
        self.ds = self.tube.get_stuff(self.x)
        self.p = self.get_pressure(self.q)
        f_left, f_right = get_f(self)
        S_left = self.tube.get_S(self.x[:-1])
        S_right = self.tube.get_S(self.x_right[:-1])
        self.h = self.get_arr_h()
        df = [self.h[i] * self.dx - (f_right[i] * S_right - f_left[i] * S_left) for i in range(len(f_right))]
        return df

    def euler_step_new(self, tau, x_left, V_left, x_right, V_right, isignite):
        """
        Шаг вперед по времени
        :param p_right: давление справа от правой границы
        :param p_left: давление слева от левой границы
        :param l: слой
        :param tau: шаг по времени
        :return: слой на следующем шаге
        """
        l1 = self.clone(self)
        l1.stretch_me_new(tau, x_left, V_left, x_right, V_right)
        l1.W = l1.tube.get_W(l1.x)
        df = self.get_dQs()
        l1.q = [(self.q[num] * self.W + tau * df[num]) / l1.W for num in range(len(self.q) - 1)]
        if isignite:
            l1.q.append((self.q[3] * self.W + tau * df[3]) / l1.W)
        else:
            l1.q.append(self.q[3])
        l1.p = l1.get_pressure(l1.q)
        return l1

    def move_to(self, tau, x_left, V_left, x_right):
        l1 = self.clone(self)
        l1.stretch_me_new(tau, x_left, V_left, x_right, V_left)
        l1.u[:] = V_left
        l1.init_arr_q()
        l1.p = l1.get_pressure(l1.q)
        return l1


    def stretch_me_new(self, tau, x_left, V_left, x_right, V_right):
        """
        Прибавление общего счетчика времени, пересчет скоростей и координат узлов
        :param p_right: давление справа от правой границы
        :param p_left: давление слева от левой границы
        :param tau: шаг по времени
        :return:
        """
        self.time += tau
        self.x = np.linspace(x_left, x_right, self.n + 1)
        self.x_right = np.roll(self.x, -1)
        self.dx = (self.x_right - self.x)[:-1]
        self.x_c = ((self.x + self.x_right) / 2)[:-1]
        self.V = np.linspace(V_left, V_right, self.n + 1)
