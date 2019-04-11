import os
import sys
from math import *
wd = os.path.abspath(__file__) # os.path.dirname(os.path.dirname(os.getcwd()))
wd = os.path.dirname(os.path.dirname(wd))
sys.path.append(wd+"\\src\\")
import numpy as np
import pytest
from pytest import approx

from tube import Tube
from gaslayer import GasLayer, GasEOS, GridStrecher, GasFluxCalculator

def test_adapt():
    gs = GridStrecher()
    n = 100#np.random.randint(10,100)
    for i in range(1000):
        n = np.random.randint(10,1000)
        xx1 = np.random.uniform(-1000,2000)
        l = np.random.uniform(0.001, 2000)
        xs = np.linspace(xx1, xx1+l, n)
        xs_borders = np.linspace(xx1, xx1+l, n+1)

        ys = 100000*np.sin(xs)  
        xs_adapted = np.zeros_like(xs_borders)

        gs.adaptine_borders(xs_borders, ys, xs_adapted)
        assert xs_adapted[0] == approx(xs_borders[0])
        assert xs_adapted[-1] == approx(xs_borders[-1])
        for x1, x2 in zip(xs_adapted[:-1], xs_adapted[1:]):
            assert x2 > x1


def _plot():
    import matplotlib.pyplot as plt

    n = 1000
    xs = np.linspace(-10, 10, n)
    xs += np.random.uniform(-0.01,0.01,n)

    ys = np.sin(xs/2) + 0.3*np.random.normal(size=n)
    plt.plot(xs, ys)    
    
    ys_smoothed = np.zeros_like(ys)

    gs = GridStrecher()
    gs.smooth_arr(xs, ys, ys_smoothed,0.01)
    plt.plot(xs, ys_smoothed)
    
    gs.smooth_arr(xs, ys, ys_smoothed,0.1)
    plt.plot(xs, ys_smoothed)

    gs.smooth_arr(xs, ys, ys_smoothed,0.2)
    plt.plot(xs, ys_smoothed)

    plt.show()

def _plot2():
    import matplotlib.pyplot as plt

    n = 1000
    xs = np.linspace(-10, 10, n)
    xs_borders = np.linspace(-10, 10, n+1)

    ys = 100000*np.sin(xs)
    plt.plot(xs, ys)    
    
    xs_adapted = np.zeros_like(xs_borders)

    gs = GridStrecher()
    gs.adaptine_borders(xs_borders, ys, xs_adapted)
    plt.scatter(xs_adapted, np.zeros_like(xs_adapted))

    plt.show()

if __name__ == "__main__":
    _plot()
    _plot2()