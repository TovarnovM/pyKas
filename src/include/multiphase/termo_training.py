import numpy as np

class Termo1d:
    def __init__(self, ys, Ts, delta_b, c_b, lambda_b, time):
        self.ys = ys
        self.Ts = Ts
        self.delta_b = delta_b
        self.c_b = c_b
        self.lambda_b = lambda_b
        self.time = time

    def copy(self, create_new_arrs=True):
        ys = np.array(self.ys) if create_new_arrs else self.ys
        Ts = np.array(self.Ts) if create_new_arrs else self.Ts
        return Termo1d(ys, Ts, self.delta_b, self.c_b, self.lambda_b, self.time)

    def step_up(self, tau, q0, q1):
        

    