import numpy as np

class Termo1d:
    @classmethod
    def get_standart(cls, **kw):
        r_0 = kw.get('r_0', 0.023/2)
        r_1 = kw.get('r_1', 0.03/2)
        n = kw.get('n', 33)
        T_0 = kw.get('T_0', 293)
        s = r_1 - r_0
        q = kw.get('q', 0.85)
        b = s*(1-q)/(1-q**(n-1))
        rs = np.zeros(n)
        rs[-1] = r_1
        rs[-2] = r_1 - b
        for i in range(rs.shape[0]-3,0,-1):
            b *= q
            rs[i] = rs[i+1] - b
        rs[0] = r_0
        Ts = np.zeros(n)
        Ts[:] = T_0
        return cls(rs, Ts,  
            delta_b=kw.get('delta_b', 7800), 
            c_b=kw.get('c_b', 480),
            lambda_b=kw.get('lambda_b', 27),
            time=kw.get('time', 0.0))

    def __init__(self, rs, Ts, delta_b, c_b, lambda_b, time):
        self.rs = rs
        self.Ts = Ts
        self.delta_b = delta_b
        self.c_b = c_b
        self.lambda_b = lambda_b
        self.time = time

    def copy(self, create_new_arrs=True):
        rs = np.array(self.rs) if create_new_arrs else self.rs
        Ts = np.array(self.Ts) if create_new_arrs else self.Ts
        return Termo1d(rs, Ts, self.delta_b, self.c_b, self.lambda_b, self.time)

    def step_up(self, tau, q_0, T_up):
        res = self.copy()
        alphas = res.Ts
        bettas = res.rs
        a = self.lambda_b/(self.c_b*self.delta_b)
        B = self.rs[0] - self.rs[1]
        C = -B
        F = -q_0/self.lambda_b
        alphas[0] = -C/B
        bettas[0] = F/B
        for i in range(1 ,self.rs.shape[0]-1):
            T_k_i = self.Ts[i]
            r_k = self.rs[i]
            r_km1=self.rs[i-1]
            r_kp1=self.rs[i+1]
            m_1 = a*(r_kp1+r_k)/(r_kp1-r_km1)/(r_kp1-r_k)
            m_2 = a*(r_k+r_km1)/(r_kp1-r_km1)/(r_k-r_km1)
            A = m_2
            B = -m_1 - m_2 - 1/tau
            C = m_1
            F = -T_k_i/tau
            alphas[i] = -C/(A*alphas[i-1] + B)
            bettas[i] = (F - A*bettas[i-1])/(A*alphas[i-1] + B)
        res.Ts[-1] = T_up
        for i in range(self.Ts.shape[0]-2,-1,-1):
            alpha = alphas[i]
            betta = bettas[i]
            res.Ts[i] = res.Ts[i+1]*alpha + betta
        res.rs[:] = self.rs
        res.time += tau
        return res

if __name__ == "__main__":
    l1 = Termo1d.get_standart()
    for i in range(100):
        l1 = l1.step_up(0.0001, 0.2, 273)
    i=0

