# distutils: language=c++
# cython: language_level=3, boundscheck=False, nonecheck=False, cdivision=True, initializedcheck=False


from tube cimport Tube, InterpXY
from gaslayer cimport GasLayer, GasEOS, GasFluxCalculator, GridStrecher

import cython
from libc.math cimport pi, sqrt, copysign, exp, pow, fabs
import numpy as np
cimport numpy as np
import json

cpdef inline void roue_to_q_(
    double ro,
    double u,
    double e,
    double[:] q) nogil:

    q[0] = ro
    q[1] = ro * u
    q[2] = ro * (e + 0.5 * u * u)

cpdef inline (double, double, double) roue_to_q(
    double ro,
    double u,
    double e) nogil:

    return ro, ro * u, ro * (e + 0.5 * u * u)


cpdef inline (double, double, double) q123_to_roue(double q1, double q2, double q3) nogil:
    cdef double ro = q1
    cdef double u = q2 / q1
    cdef double e = q3 / q1 - 0.5 * u * u
    return ro, u, e

class Singleton(type):
    """Метакласс для объектов с паттеном синглтон
    """

    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]

class PowderBD(metaclass=Singleton):
    def __init__(self):
        import os
        with open(os.path.dirname(os.path.abspath(__file__))+'\\gpowders.json') as f:
            self.all_powders_dict = json.load(f)


cdef class Powder(GasEOS):
    @classmethod
    def from_bd(cls, dict_or_name, n_points=300):
        """Создать образец Powder из словаря БД
        
        Arguments:
            d {[type]} -- словарь типа 
                        {   'name': '4\\7',
                            'f': 1.027,
                            'etta': 0.228,
                            'alpha_k': 1.008,
                            'T_1': 3006.0,
                            'ro': 1.6,
                            'I_k': 0.32,
                            'Z_k': 1.488,
                            'k_1': 0.811,
                            'lambda_1': 0.081,
                            'k_2': 0.505,
                            'lambda_2': -1.024,
                            'k_f': 0.0003,
                            'k_l': 0.0016   }


                        ИЛИ это просто имя пороха типа r'4\7'
        """
        if isinstance(dict_or_name, str):
            dict_or_name = PowderBD().all_powders_dict[dict_or_name]
        Z_k = dict_or_name['Z_k']
        k1 = dict_or_name['k_1']
        l1 = dict_or_name['lambda_1']
        k2 = dict_or_name['k_2']
        l2 = dict_or_name['lambda_2']
        zs = np.linspace(0, Z_k, n_points)
        dpsidz = np.empty_like(zs)
        for i, z in enumerate(zs):
            dpsidz[i] = k1*(1+2*l1*z) if z<=1 else k2*(1+2*l2*(z-1)) if z<Z_k else 0
        return cls(name=dict_or_name['name'], I_k=dict_or_name['I_k'], alpha_k=dict_or_name['alpha_k'], 
                   ro=dict_or_name['ro'], f=dict_or_name['f'], etta=dict_or_name['etta'], 
                   T_1=dict_or_name['T_1'], zs=zs, dpsi_dz=dpsidz)

    def __init__(self, name, I_k, alpha_k, ro, f, etta, T_1, zs, dpsi_dz):
        """Конструктор
        
        Arguments:
            I_k  -- импульс конца горения МПа*с
            alpha_k  -- коволюм (как в БД, т.е. настоящий коволюм =/ 1000)
            ro  -- плотность пороха г/см^3
            f  -- сила пороха, МДж/кг
            etta  -- k -1
            T_1  -- темп. горения ?
            zs {list} -- точки z для интерполяции dpsi
            dpsi_dz {list} -- точки dpsi_dz для интерполяции dpsi_dz
        """
        self.name = name
        super().__init__(gamma=etta+1, kappa=alpha_k/1000, p_0=0, c_0=0, kind=33)
        self.I_k = I_k * 1e6
        self.ro = ro * 1000
        self.f = f * 1e6
        self.T_1 = T_1
        self.R_g = self.f / self.T_1
        self.nu = 1
        
        psis = [0]
        dpsis = [dpsi_dz[0]]
        for i in range(1, len(dpsi_dz)):
            dz = zs[i] - zs[i-1]
            dpsi = (dpsi_dz[i] + dpsi_dz[i-1])*0.5
            psis.append(psis[-1] + dz*dpsi)
            dpsis.append(dpsi_dz[i])
        mnj = 1/psis[-1]
        for i in range(len(zs)):
            psis[i] = psis[i]*mnj
            dpsis[i] = dpsis[i]*mnj
        self.dpsi = InterpXY(zs, dpsis)    
        self.psi = InterpXY(zs, psis)
        self.z_k = self.get_z_k()
        
    
    cpdef double get_e_powder(self, double ro, double p, double z):
        cdef double psi = self.psi.get_v(z)
        cdef double e = (p/(self.gamma-1)) * (1/ro - ((1-psi)/self.ro + self.kappa*psi)) +(1-psi)*self.f/(self.gamma-1)
        return e if e > 0 else 0
    
    cpdef double get_p_powder(self, double ro, double e, double z):
        cdef double psi = self.psi.get_v(z)
        cdef double chsl = e * (self.gamma-1) - (1-psi)*self.f
        cdef double znam = 1/ro - ((1-psi)/self.ro + self.kappa*psi)
        if fabs(znam) < 1e-13:
            return 0
        cdef double p = chsl/znam
        return p if p > 0 else 0
    
    cpdef double get_csound_powder(self, double ro, double p, double z): 
        cdef double psi = self.psi.get_v(z)
        cdef double chsl = self.gamma * p
        cdef double znam = 1/ro - ((1-psi)/self.ro + self.kappa*psi)
        if znam < 1e-13:
            return sqrt(self.gamma*p/ro)
        return sqrt(chsl/znam)/ro

    cpdef double get_z_k(self):
        return self.dpsi.xs[-1]


cdef class PowderOvLayer(GasLayer):
    @classmethod
    def get_standart(cls, tube: Tube, x_1: float, powder_layer_dict, calc_settings):
        """Получить и проиницилаизировать стандартный, равномерный полиэтиленовый слой
        
        Arguments:
            tube {Tube} -- труба)
            x_1 {float} -- левая граница полиэтиленовой области
            powder_layer_dict {dict} -- powder_layer_dict_sample = {
                                            'type': 'powder',
                                            'powder': {
                                                'name': '4\\7',
                                                'f': 1.027,
                                                'etta': 0.228,
                                                'alpha_k': 1.008,
                                                'T_1': 3006.0,
                                                'ro': 1.6,
                                                'I_k': 0.32,
                                                'Z_k': 1.488,
                                                'k_1': 0.811,
                                                'lambda_1': 0.081,
                                                'k_2': 0.505,
                                                'lambda_2': -1.024,
                                                'k_f': 0.0003,
                                                'k_l': 0.0016
                                            },
                                            'omega': 35,  # грамм
                                            'delta': 700, # г/cm^3
                                            'p_0': 5e6, # начальное давление
                                            't_ign': 0.00001, # начало горения
                                            'u_0': 0,     #начальная скорость
                                        }
            calc_settings {dict} -- calc_settings_sample = {
                                        'cell_dx': 0.0025,
                                        'n_cells_min': 13,
                                        'n_cells_max': 300,
                                        'GasFluxCalculator_kwargs': {},
                                        'GridStrecher_kwargs': {}
                                    }
        """
        def get_n_cells(x_1, x_2, calc_settings):
            n_cells = round(abs(x_2-x_1)/calc_settings['cell_dx'])
            n_cells = min(calc_settings['n_cells_max'], n_cells)
            n_cells = max(calc_settings['n_cells_min'], n_cells)
            return n_cells
        powder = Powder.from_bd(powder_layer_dict['powder'])
        flux_calculator = GasFluxCalculator(**calc_settings['GasFluxCalculator_kwargs'])
        grid_strecher = GridStrecher(**calc_settings['GridStrecher_kwargs'])
        omega = powder_layer_dict['omega']*1e-3
        delta = powder_layer_dict['delta']
        w = omega/delta
        x_2 = tube.get_x2(x_1, w)
        n_cells = get_n_cells(x_1, x_2, calc_settings)
        pow_layer = PowderOvLayer(n_cells=n_cells, tube=tube, powder=powder, 
                                flux_calculator=flux_calculator, grid_strecher=grid_strecher,
                                t_ign=powder_layer_dict['t_ign'])
        pow_layer.xs_borders = np.linspace(x_1, x_2, n_cells+1)
        u_0 = powder_layer_dict['u_0']
        def foo_ropu(*args):
            return delta, powder_layer_dict['p_0'], u_0, 0
        pow_layer.init_ropue_fromfoo(foo_ropu)
        return pow_layer

    def __init__(self, n_cells, Tube tube, Powder powder, GasFluxCalculator flux_calculator, GridStrecher grid_strecher, double t_ign=0):
        super().__init__(n_cells, tube, powder, flux_calculator, grid_strecher, 4)
        self.zs = np.zeros(n_cells, dtype=np.double)
        self.t_ign = t_ign
        self.color_4_plot = '#6b400e'

    cpdef void copy_params_to_Ov(self, PowderOvLayer to_me):
        self.copy_params_to(to_me)
        to_me.zs[:] = self.zs
        to_me.t_ign = self.t_ign

    cpdef GasLayer copy(self):
        cdef PowderOvLayer res = PowderOvLayer(self.n_cells, self.tube, <Powder>self.gasEOS, self.flux_calculator, self.grid_strecher)
        self.copy_params_to_Ov(res)
        return res
    
    cpdef void init_ropue_fromfoo(self, foo_ropu, bint init_q=True,  bint init_SsdW=True):
        self.grid_strecher.fill_xs_cells(self.xs_borders, self.xs_cells)
        cdef size_t i
        cdef double x
        cdef Powder eos = <Powder>self.gasEOS
        cdef double ro, p, u, e, z
        for i in range(self.n_cells):
            ro, p, u, z = foo_ropu(self.xs_cells[i], self)
            self.ros[i], self.ps[i], self.us[i], self.zs[i] = ro, p, u, z
            self.es[i] = eos.get_e_powder(ro, p, z)
            self.cs[i] = eos.get_csound_powder(ro, p, z)
        if init_q:
            self.init_q()
        if init_SsdW:
            self.init_SsdW()
        self.init_taus_acustic()

    cpdef void init_q(self):
        cdef size_t i
        cdef size_t n = self.ros.shape[0]
        for i in range(n):
            q1, q2, q3 = roue_to_q(self.ros[i], self.us[i], self.es[i])
            self.qs[0, i] = q1
            self.qs[1, i] = q2 
            self.qs[2, i] = q3
            self.qs[3, i] = self.ros[i]*self.zs[i]

    cpdef void init_ropue(self):
        cdef size_t i
        cdef size_t n = self.ros.shape[0]
        cdef double ro, p, u, e, z
        cdef Powder eos = <Powder>self.gasEOS
        cdef double z_k = eos.z_k
        for i in range(n):
            ro, u, e = q123_to_roue(self.qs[0, i], self.qs[1, i], self.qs[2, i])
            z = self.qs[3, i] / ro
            if z > z_k:
                z = z_k
            p = eos.get_p_powder(ro, e, z)
            self.ros[i] = ro
            self.ps[i] = p
            self.us[i] = u
            self.es[i] = e
            self.zs[i] = z
            self.cs[i] = eos.get_csound_powder(ro, p, z)

    cpdef void init_h(self):
        cdef size_t i
        cdef Powder powder = <Powder>self.gasEOS
        for i in range(self.hs.shape[1]):
            self.hs[0, i] = 0
            self.hs[1, i] = self.ps[i] * self.ds[i] 
            self.hs[2, i] = 0
            self.hs[3, i] = 0
        if self.time >= self.t_ign:
            for i in range(self.hs.shape[1]):
                self.hs[3, i] = self.ros[i] * 0.5*(self.S[i]+self.S[i+1]) * pow(self.ps[i], powder.nu) / powder.I_k

    cpdef void fill_fluxes(self):
        cdef size_t i
        cdef Powder powder = <Powder>self.gasEOS
        cdef double rayE, rayU, rayR, rayP, rayZ, M_j, J_j, E_j, ray_W
        cdef double[:] rr_val
        cdef double Mrf, pf, cs, r1, r2, u1, u2, H1, H2, vbi
        cdef double z1, z2
        if self.flux_calculator.flux_type == 1:
            for i in range(self.fluxes.shape[1]):
                ray_W = self.Vs_borders[i]
                rayU = self.fluxes[0, i]
                rayR = self.fluxes[1, i]
                rayP = self.fluxes[2, i]

                if i == 0:
                    rayZ = self.zs[0]
                elif i == self.fluxes.shape[1]-1:
                    rayZ = self.zs[self.zs.shape[0]-1]
                elif ray_W > self.U_kr[i]:
                    rayZ = self.zs[i]
                else:
                    rayZ = self.zs[i-1]

                rayE = powder.get_e_powder(rayR, rayP, rayZ)

                M_j = rayR * (rayU - ray_W)
                J_j = rayP + M_j * rayU
                E_j = (rayE + 0.5*rayU*rayU)*M_j + rayP*rayU
                self.fluxes[0, i] = M_j
                self.fluxes[1, i] = J_j
                self.fluxes[2, i] = E_j
                self.fluxes[3, i] = M_j * rayZ
        if self.flux_calculator.flux_type == 2:
            z1 = self.zs[0]
            for i in range(self.fluxes.shape[1]):
                rr_val = self.rr_vals[i]
                Mrf = rr_val[0]
                pf = rr_val[1]
                cs = rr_val[2]
                r1 = rr_val[3]
                r2 = rr_val[4]
                u1 = rr_val[5] 
                u2 = rr_val[6]
                H1 = rr_val[7]
                H2 = rr_val[8]
                vbi = rr_val[9]
                if i != self.fluxes.shape[1]-1:
                    z2 = self.zs[i]

                self.fluxes[0, i] = 0.5*(cs*Mrf*(r1+r2)-cs*abs(Mrf)*(r2-r1))
                self.fluxes[1, i] = 0.5*(cs*Mrf*(r1*u1+r2*u2)-cs*abs(Mrf)*(r2*u2-r1*u1)) + pf
                self.fluxes[2, i] = 0.5*(cs*Mrf*(r1*H1+r2*H2)-cs*abs(Mrf)*(r2*H2-r1*H1)) + pf*vbi
                self.fluxes[3, i] = 0.5*(cs*Mrf*(r1*z1+r2*z2)-cs*abs(Mrf)*(r2*z2-r1*z1))

                z1 = z2

    def to_dict(self):
        res = super().to_dict()
        res['zs'] = np.array(self.zs)
        return res

    def from_dict(self, d):
        super().from_dict()
        self.zs = np.array(d['zs'])

    def __str__(self):
        return super().__str__() + f"\n        {{ 'powder': r'{self.gasEOS.name}', 'z_max': {np.max(self.zs)}, 't_ign'={self.t_ign} }}"