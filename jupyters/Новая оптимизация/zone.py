
from math import *
import numpy.random as rnd
# import matplotlib.pyplot as plt
import numpy as np

class Zone:
    @classmethod
    def get_x2min(cls, x1, r1, r2, alpha_max=30):
        k = tan(-pi * self.alpha_max / 180)
        b = self.r1 - k * self.x1 
        return (self.r2 - b)/k
    
    def __init__(self, x1, x2, r1, r2, alpha_max=30):
        self.x1 = min(x1, x2)
        self.x2 = max(x1, x2)
        self.r1 = max(r1, r2)
        self.r2 = min(r1, r2)
        self.alpha_max = alpha_max
        
    def plot(self, ax):
        if self.is_valid():
            xs = np.array([self.x1, 0, self.x2, 0.0, self.x1])
            ys = np.array([self.r1, self.r2, self.r2, self.r1, self.r1], dtype=float)
            
            k = tan(-pi * self.alpha_max / 180)
            b = self.r1 - k * self.x1 
            xs[1] = (self.r2 - b)/k
            print(k, b)
            b = self.r2 - k * self.x2 
            xs[3] = (self.r1 - b)/k
            print(k, b)
            ax.plot(xs, ys)
        
    def plot_diag(self, ax):
        if self.is_valid():
            xs = np.array([self.x1, self.x2])
            ys = np.array([self.r1, self.r2], dtype=float)
        
            ax.plot(xs, ys)
    
    def is_valid(self):
        alpha = atan((self.r1 - self.r2)/(self.x2 - self.x1)) * 180 / pi
        return alpha <= self.alpha_max
    
    def get_rnd_x_r(self):
        x = rnd.uniform(self.x1, self.x2)
        rmin, rmax = self.get_rmin_rmax(x)
        
        r = rnd.uniform(self.r2, self.r1)
        while r > rmax or r < rmin:
            x = rnd.uniform(self.x1, self.x2)
            rmin, rmax = self.get_rmin_rmax(x)
            r = rnd.uniform(self.r2, self.r1)
        return x, r
        
    def get_rmin_rmax(self, x):
        k = -tan(pi * self.alpha_max / 180)
        b = self.r1 - k * self.x1 
        r = k * x + b
        rmin = max(r, self.r2)
        
        b = self.r2 - k * self.x2
        r = k * x + b
        rmax = min(r, self.r1)
        
        return rmin, rmax
    
    def is_in(self, x, r):
        if self.x1 <= x <= self.x2:
            rmin, rmax = self.get_rmin_rmax(x)
            return rmin <= r <= rmax
        return False
    
    def get_split_zones(self, x, r):
        if self.is_in(x, r):
            return Zone(self.x1, x, self.r1, r, self.alpha_max), Zone(x, self.x2, r, self.r2, self.alpha_max)
        
    @property
    def area(self):
        return float((self.x2 - self.x1) * (self.r1 - self.r2))

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    z = Zone(4,10, 1,2, alpha_max=30)
    fig, ax = plt.subplots()
    z.plot(ax)
    in_xs, in_ys, out_xs, out_ys = [],[],[],[]
    for i in range(1000):
        x, r = z.get_rnd_x_r()
        in_xs.append(x)
        in_ys.append(r)
    ax.scatter(in_xs, in_ys)
    plt.show()
    
    z = Zone(4,10, 1,2, alpha_max=15)
    fig, ax = plt.subplots()
    z.plot(ax)
    in_xs, in_ys, out_xs, out_ys = [],[],[],[]
    for i in range(1000):
        x = rnd.uniform(4, 10)
        r = rnd.uniform(1, 2)
        if z.is_in(x, r):
            in_xs.append(x)
            in_ys.append(r)
        else:
            out_xs.append(x)
            out_ys.append(r)

    ax.scatter(in_xs, in_ys)
    ax.scatter(out_xs, out_ys)
    plt.show()
    
    z = Zone(4.0 ,10, 1,2, alpha_max=45)
    fig, ax = plt.subplots()
    zs = [z]
    for i in range(4):
        areas = np.array([z.area for z in zs])
        areas /= np.sum(areas)
        i = rnd.choice(range(len(zs)), p=areas)
        z = zs.pop(i)
        x, r = z.get_rnd_x_r()
        z1, z2 = z.get_split_zones(x, r)
        zs.append(z1)
        zs.append(z2)
    for z in zs:
        z.plot_diag(ax)

    plt.plot()
