from gaslayer import GasEOS, GasLayer


class MixEOS(GasEOS):
    def __init__(self, v_H2, v_O2, v_H20, v_He) -> None:
        super().__init__(1.4, kappa=0.0, p_0=0.0, c_0=0.0, kind=1)
        self.v_H2 = v_H2
        self.v_O2 = v_O2
        self.v_H20 = v_H20
        self.v_He = v_He
        self.set_param_before()

    def set_param_before(self):
        pass

class ReactantGasLayer(GasLayer):
    @classmethod
    def get_standart(tube, x_left, lr_dict, calc_settings):
        pass

    def __init__(self, n_cells, tube, flux_calculator, grid_strecher, gasEOS):
        super().__init__(n_cells, tube, gasEOS, flux_calculator, grid_strecher, 3)


