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
import pprint as pp

def test_integrate():
    xs = [1,2,3,4,5]
    ys = [1,2,0,-2,-1]
    interp = InterpXY(xs, ys)
    assert 0 == approx(interp.integrate(0,6))

def test_integrate2():
    xs = [1,2,3,4,5]
    ys = [1,2,0,-2,-1]
    interp = InterpXY(xs, ys)
    assert 10 == approx(interp.integrate(-10,0))


def test_integrate3():
    xs = [1,2,3,4,5]
    ys = [1,2,0,-2,-1]
    interp = InterpXY(xs, ys)
    assert -10 == approx(interp.integrate(10,20))

def test_integrate4():
    xs = [1,2,3,4,5]
    ys = [1,2,0,-2,-1]
    interp = InterpXY(xs, ys)
    assert 0 == approx(interp.integrate(1.5,4.5))

def test_get_xs_zeros():
    xs = [1,2,3,4,5]
    ys = [1,2,0,-2,-1]
    interp = InterpXY(xs, ys)
    zeros = interp.get_xs_zeros()
    assert [3] == approx(zeros)

def test_get_xs_zeros2():
    xs = [1,2,3,4,5]
    ys = [-2,2,0,-2,-1]
    interp = InterpXY(xs, ys)
    zeros = interp.get_xs_zeros()
    assert [1.5, 3] == approx(zeros)

def test_get_xs_zeros3():
    xs = [1,2,3,4,5]
    ys = [2,2,0,2,1]
    interp = InterpXY(xs, ys)
    zeros = interp.get_xs_zeros()
    assert [3] == approx(zeros)

def test_get_xs_zeros4():
    xs = [1,2,3,4,5]
    ys = [2,2,1,2,1]
    interp = InterpXY(xs, ys)
    zeros = interp.get_xs_zeros()
    assert zeros is None

def _main():
    xs = [-2,2,3,4,9]
    ys = [1,5,3,1,2]
    interp1 = InterpXY(xs, ys)

    xs = [-1,2,3,5.1,7]
    ys = [-1,2,6,-3,0]
    interp2 = InterpXY(xs, ys) + (-4,3)


    fig, ax = plt.subplots()
    interp1.plot_diff(interp2, fig, ax)
    # interp1.plot(fig, ax, marker='o')
    # interp2.plot(fig, ax, marker='x')
    # abs(interp1-interp2).plot(fig, ax, marker='*')
    plt.grid(True)
    pp.pprint(InterpXY.get_diff_report(interp1, interp2))
    plt.show()

def _main2():
    a1 = np.array([1,2,3,4.1,-1])
    a2 = np.array([4,4.100000001])
    print(a1[a1.argsort()])

if __name__ == "__main__":
    _main()