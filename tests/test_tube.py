import sys
import os
sys.path.append(os.getcwd()+"\\src\\")
from tube import Tube, InterpXY
from math import pi, sqrt
from pytest import approx
import numpy as np
import random

def test_tube_exists():
    t = Tube([1,2,3], [1,0,0])
    assert t is not None

def test_InterpXY():
    ixy = InterpXY([1,2,3], [4,5,6])
    assert ixy.get_vs(np.array([-10,1,1.5,2,2.25,2.5,3,33])) == approx([4,4,4.5,5,5.25,5.5,6,6.])

def test_tube_d():
    ds = np.array([1,1,0.5])
    xs = [1,2,4]
    t = Tube(xs, ds)
    assert t.get_ds() == approx(ds)
    assert t.get_xs() == approx(xs)
    assert t.get_d(3) == approx(0.75)

def test_tube_s():
    d = sqrt(4 / pi) # s = 1
    ds = np.array([d,d,0.5 * d])
    xs = [-1,2,4]
    t = Tube(xs, ds)
    assert t.get_S(np.array([-10,-1,2], dtype=np.double)) == approx([1,1,1])
    assert t.get_S(np.array([2,4,55], dtype=np.double)) == approx([1,pi*(0.5*d)**2/4, pi*(0.5*d)**2/4])
    assert t.get_S(np.array([3], dtype=np.double)) == approx([pi*(0.5*d)**2/4*0.5+0.5])

def test_tube_ws():
    xs = [1,2,4]
    ds = [1,1,1]
    s = pi/4
    t = Tube(xs, ds)
    for i in range(33):
        x1 = -10 + random.random() * 30
        x2 = -10 + random.random() * 30
        x1, x2 = min(x1, x2), max(x1, x2)
        assert t.get_W_between(x1, x2) == approx(s*(x2-x1))

def test_tube_dsdx_1():
    xs = [1,2,4]
    ds = [1,1,1]
    s = pi/4
    t = Tube(xs, ds)
    for i in range(33):
        x1 = -10 + random.random() * 30
        x2 = -10 + random.random() * 30
        x1, x2 = min(x1, x2), max(x1, x2)
        assert t.get_dsdx(np.array([x1, x2])) == approx([0.0])

def test_tube_dsdx_2():
    xs = [0,2,4,6]
    ds = [1,2,4,2]
    t = Tube(xs, ds)
    for i in range(33):
        x1 = 0 + random.random() * 2
        x2 = 0 + random.random() * 2
        x1, x2 = min(x1, x2), max(x1, x2)
        assert t.get_dsdx(np.array([x1, x2])) == approx([(pi*2**2/4-pi*1**2/4)/2])
    for i in range(33):
        x1 = 2 + random.random() * 2
        x2 = 2 + random.random() * 2
        x1, x2 = min(x1, x2), max(x1, x2)
        assert t.get_dsdx(np.array([x1, x2])) == approx([(pi*4**2/4-pi*2**2/4)/2])
    for i in range(33):
        x1 = 4 + random.random() * 2
        x2 = 4 + random.random() * 2
        x1, x2 = min(x1, x2), max(x1, x2)
        assert t.get_dsdx(np.array([x1, x2])) == approx([(pi*2**2/4-pi*4**2/4)/2])

if __name__ == "__main__":
    xs = [1,2,4]
    ixy = InterpXY([1,2,3], [4,5,6])
    ans = np.array(ixy.get_vs(np.array([-1,1,1.5,2,2.25,2.5,3,33])))
    print(ans)
    print([4,4,4.5,5,5.25,5.5,6,6.])
    # ds = [1,1,1]
    # s = pi/4
    # t = Tube(xs, ds)
    # for i in range(33):
    #     x1 = -10 + random.random() * 30
    #     x2 = -10 + random.random() * 30
    #     x1, x2 = min(x1, x2), max(x1, x2)
    #     print( np.asarray(t.get_dsdx(np.array([x1, x2])) ))