# distutils: language=c++
# cython: language_level=3, boundscheck=False, nonecheck=False, cdivision=True, initializedcheck=False


from tube cimport Tube, InterpXY
from gaslayer cimport GasLayer, GasEOS, GasFluxCalculator, GridStrecher
from godunov cimport get_p_0

import cython
from libc.math cimport pi, sqrt, copysign, exp
import numpy as np
cimport numpy as np


cdef class ElPistLayer(GasLayer):
    @classmethod
    def get_standart(cls, tube: Tube, x_1: float, pist_layer_dict, calc_settings):
        """Получить и проиницилаизировать стандартный, равномерный полиэтиленовый слой
        
        Arguments:
            tube {Tube} -- труба)
            x_1 {float} -- левая граница полиэтиленовой области
            pist_layer_dict {dict} -- pist_layer_dict_sample = {
                                                'type': 'pist',
                                                'omega': 0.03,  # кг
                                                'p_0': 1e5, # начальное давление в полиэтилене
                                                'u_0': 0,     #начальная скорость
                                                'ElPistEOS_kwargs': dict(
                                                    k=1.63098, 
                                                    c_0=2308, 
                                                    ro_0=919.03, 
                                                    sigma_star=25.2, 
                                                    k_0=0.054, 
                                                    b_1=0.027, 
                                                    b_2=0.00675, 
                                                    tau_0=1.36, 
                                                    mu=0.001, 
                                                    tau_s=1,                 
                                                    zeroP=False, 
                                                    zeroE=True)
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
        flux_calculator = GasFluxCalculator(**calc_settings['GasFluxCalculator_kwargs'])
        grid_strecher = GridStrecher(**calc_settings['GridStrecher_kwargs'])
        epistEOS = ElPistEOS(**pist_layer_dict['ElPistEOS_kwargs'])
        ro_0 = pist_layer_dict['ElPistEOS_kwargs']['ro_0']
        
        W0 = pist_layer_dict['omega']/ro_0
        x_2 = tube.get_x2(x_1, W0)
        n_cells = get_n_cells(x_1, x_2, calc_settings)
        pist_layer = cls(n_cells=n_cells, tube=tube, epistEOS=epistEOS, flux_calculator=flux_calculator, grid_strecher=grid_strecher)
        pist_layer.xs_borders = np.linspace(x_1, x_2, n_cells+1)
        u_0 = pist_layer_dict['u_0']
        def foo_ropu(*args):
            return ro_0, pist_layer_dict['p_0'], u_0
        pist_layer.init_ropue_fromfoo(foo_ropu)
        return pist_layer

    def __init__(self, n_cells, Tube tube, ElPistEOS epistEOS, GasFluxCalculator flux_calculator, GridStrecher grid_strecher):
        super().__init__(n_cells, tube, epistEOS, flux_calculator, grid_strecher, 3)
        self.tauxx_flux = np.zeros(n_cells+1, dtype=np.double)
        self.color_4_plot = '#0022cf'

    cpdef double get_p_left(self):
        return self.ps[0] - self.tauxx_flux[0]

    cpdef double get_p_right(self):
        return self.ps[-1]- self.tauxx_flux[-1]

    cpdef void copy_params_to_elpist(self, ElPistLayer lr):
        self.copy_params_to(lr)
        lr.tauxx_flux[:] = self.tauxx_flux

    cpdef GasLayer copy(self):
        cdef ElPistLayer res = ElPistLayer(self.n_cells, self.tube, self.gasEOS, self.flux_calculator, self.grid_strecher)
        self.copy_params_to_elpist(res)
        return res

    cpdef void init_h(self):
        cdef double S, R, sigma_nt_w, sigma_nn_w, tauxx, dudx, skobka, h
        cdef size_t i
        cdef ElPistEOS eos = <ElPistEOS>self.gasEOS
        for i in range(self.hs.shape[1]):
            S = 0.5*(self.S[i]+self.S[i+1])
            R = sqrt(S/pi)
            sigma_nt_w = copysign(eos.get_tauu(self.ps[i], self.us[i]), -self.us[i])

            if i == self.hs.shape[1]-1:
                dudx = (self.us[i] - self.us[i-1])/(self.xs_cells[i] - self.xs_cells[i-1])
            else:
                dudx = (self.us[i+1] - self.us[i])/(self.xs_cells[i+1] - self.xs_cells[i])
            skobka = 2*dudx - self.us[i]*self.ds[i]/S
            h = 1/2/sqrt(3)*abs(skobka)
            if abs(h) < 1e-10:
                tauxx = 0.0
            else:
                tauxx = 2/3*eos.get_kh(h)*skobka
            sigma_nn_w = - self.ps[i] - 0.5 * tauxx 
            self.hs[0, i] = 0
            self.hs[1, i] = 2*pi*R*sigma_nt_w - self.ds[i]*sigma_nn_w
            self.hs[2, i] = 2*pi*R*sigma_nt_w*self.us[i]

            self.tauxx_flux[i] = tauxx
        self.tauxx_flux[self.tauxx_flux.shape[0]-1] = tauxx

    cpdef void fill_fluxes(self):
        cdef size_t i
        cdef double rayE, rayU, rayR, rayP, M_j, J_j, E_j, ray_W, tauxx
        for i in range(self.fluxes.shape[1]):
            ray_W = self.Vs_borders[i]
            rayU = self.fluxes[0, i]
            rayR = self.fluxes[1, i]
            rayP = self.fluxes[2, i]
            rayE = self.gasEOS.get_e(rayR, rayP)
            tauxx = self.tauxx_flux[i]

            M_j = rayR * (rayU - ray_W)
            J_j = rayP + M_j * rayU + tauxx
            E_j = (rayE + 0.5*rayU*rayU)*M_j + (rayP + tauxx)*rayU
            self.fluxes[0, i] = M_j
            self.fluxes[1, i] = J_j
            self.fluxes[2, i] = E_j


cdef class ElPistEOS(GasEOS):
    def __init__(self, k=1.63098, c_0=2308, ro_0=919.03, sigma_star=25.2, \
                k_0=0.054, b_1=0.027, b_2=0.00675, tau_0=1.36, mu=0.001, tau_s=1, \
                zeroP=True, zeroE=True):
        
        super().__init__(gamma=k, kappa=0, p_0=get_p_0(ro_0, c_0, k), c_0=c_0, kind=2)
        self.ro_0 = ro_0
        self.sigma_star = sigma_star * 1e6
        self.k_0 = k_0
        self.b_1 = b_1
        self.b_2 = b_2
        self.tau_0 = tau_0 * 1e6
        self.mu = mu
        self.tau_s = tau_s * 1e6
        self.zeroP = zeroP
        self.zeroE = zeroE
    
    def __repr__(self):
        return f'''ElPistEOS(k={self.gamma}, c_0={self.c_0}, ro_0={self.ro_0}, 
        sigma_star={self.sigma_star*1e-6}, k_0={self.k_0}, b_1={self.b_1}, 
        b_2={self.b_2}, tau_0={self.tau_0*1e-6}, mu={self.mu}, tau_s={self.tau_s*1e-6},
        zeroP={bool(self.zeroP)}, zeroE={bool(self.zeroE)})'''

    def __str__(self):
        return repr(self)

    cpdef double get_e(self, double ro, double p):
        cdef double e = (p-self.c_0*self.c_0*(ro-self.ro_0))/((self.gamma-1)*ro)
        if self.zeroE:
            return e if e > 0 else 0
        return e

    cpdef double get_p(self, double ro, double e):
        cdef double p = e*ro*(self.gamma-1) + self.c_0*self.c_0*(ro-self.ro_0)
        if self.zeroP:
            return p if p > 0 else 0
        return p

    cpdef double get_tauu(self, double sigma, double u):
        if sigma < self.sigma_star:
            return self.k_0*(1+self.b_1*u)*exp(-self.b_2*u)*sigma
        else:
            return self.tau_0*(1+self.b_1*u)*exp(-self.b_2*u)
    
    cpdef double get_kh(self, double h):
        return self.mu + self.tau_s *0.5 / h