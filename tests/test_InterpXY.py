import sys
import os
wd = os.path.abspath(__file__) # os.path.dirname(os.path.dirname(os.getcwd()))
wd = os.path.dirname(os.path.dirname(wd))
sys.path.append(wd+"\\src\\")
from tube import Tube, InterpXY
from math import pi, sqrt
from pytest import approx
import numpy as np
import matplotlib.pyplot as plt

def _main():
    xs = [1,2,3,4,5]
    ys = [4,5,6,1,2]
    interp = InterpXY(xs, ys)
    fig, ax = plt.subplots()
    interp.plot(fig, ax, marker='o')
    plt.show()

def _main2():
    a1 = np.array([1,2,3,4.1,-1])
    a2 = np.array([4,4.100000001])
    print(np.sort(np.union1d(a1,a2)))

if __name__ == "__main__":
    _main2()