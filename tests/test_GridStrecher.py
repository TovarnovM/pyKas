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

def _plot():
    import matplotlib.pyplot as plt

    n = 100
    xs = np.linspace(-10, 10, n)
    xs += np.random.uniform(-0.1,0.1,n)

    ys = np.sin(xs/2) + 0.3*np.random.normal(size=n)
    plt.plot(xs, ys)    
    
    ys_smoothed = np.zeros_like(ys)

    gs = GridStrecher()
    gs.smooth_arr(xs, ys, ys_smoothed,0.1)
    plt.plot(xs, ys_smoothed)
    
    gs.smooth_arr(xs, ys, ys_smoothed,0.15)
    plt.plot(xs, ys_smoothed)

    gs.smooth_arr(xs, ys, ys_smoothed,0.2)
    plt.plot(xs, ys_smoothed)

    plt.show()

if __name__ == "__main__":
    _plot()