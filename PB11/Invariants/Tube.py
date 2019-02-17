"""
Created on Wed Jan 17 10:53:30 2018

Создание геометрии ствола
"""

import numpy as np
from scipy import interpolate
from math import pi


class Tube(object):
    def __init__(self, xs, ds):
        self.d = interpolate.interp1d(xs, ds, bounds_error=False, fill_value=(ds[0], ds[-1]))
        dd = np.array(ds, dtype=np.float64)
        ss = dd ** 2 * pi * 0.25
        self.s = interpolate.interp1d(xs, ss, bounds_error=False, fill_value=(ss[0], ss[-1]))

        # Something new
        def ConeW(d1, d2, h):
            return pi * h * (d1 * d1 + d1 * d2 + d2 * d2) / 12

        xss = [xs[0] - 0.1] + [x for x in xs] + [xs[-1] + 0.1]
        dss = [ds[0]] + [d for d in ds] + [ds[-1]]

        ws = [0]
        for x1, x2, d1, d2 in zip(xss, xss[1:], dss, dss[1:]):
            ws.append(ConeW(d1, d2, x2 - x1) + ws[-1])

        self.w = interpolate.interp1d(xss, ws, bounds_error=False, fill_value="extrapolate")
        self.w_reverse = interpolate.interp1d(ws, xss, bounds_error=False, fill_value="extrapolate")

    def get_stuff(self, xs):
        """
        return (ds)
        """
        x = np.array(xs, dtype=np.float64)

        x_right = np.roll(x, -1)

        dx = x_right - x

        s_l = self.s(x)
        s_r = self.s(x_right)

        ds = (s_r - s_l) / dx
        return ds[:-1]

    def get_W(self, xs):
        """
        ????????? ?????? ?????
        :param xs:
        :return:
        """
        x = np.array(xs, dtype=np.float64)

        x_right = np.roll(x, -1)

        dx = x_right - x

        s_l = self.s(x)
        s_r = self.s(x_right)

        W = (s_r + s_l + np.sqrt(s_l * s_r)) * dx / 3
        return W[:-1]

    def get_S(self, xs):
        """
        ????????? ???????? ? ?????
        :param xs:
        :return:
        """
        x = np.array(xs, dtype=np.float64)
        s = self.s(x)
        return s

    def get_W_between(self, x1, x2):
        return self.w(x2) - self.w(x1)

    def get_x2(self, x1, w):
        w1 = self.w(x1)
        return self.w_reverse(w1 + w)
