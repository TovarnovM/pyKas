import sys
import os
wd = os.path.abspath(__file__) # os.path.dirname(os.path.dirname(os.getcwd()))
wd = os.path.dirname(os.path.dirname(wd))
sys.path.append(wd+"\\src\\")
from gaslayer import foo, q_to_roue, roue_to_q_, roe_to_p, rop_to_e, AUSM_gas_, rop_to_csound
import numpy as np
from pytest import approx

def test_cimport():
    a = foo()
    assert a == 3

def test_q_to_roue():
    q = np.array([1,2,3.0])
    ro, u, e = q_to_roue(q)
    assert ro == approx(q[0])
    assert u == approx(q[1]/q[0])
    assert e == approx(q[2]/q[0] - 0.5 * u * u)

def test_roue_to_q_():
    q = np.array([0,0,0.0])
    ro, u, e = 1,2,3
    roue_to_q_(ro, u, e, q)
    assert ro == approx(q[0])
    assert u == approx(q[1]/q[0])
    assert e == approx(q[2]/q[0] - 0.5 * u * u)

def test_roe_to_p():
    ro, e, gamma, b = 1, 2, 1.2, 0.001
    p = roe_to_p(ro, e, gamma, b)

    e2 = rop_to_e(ro, p , gamma, b)
    assert e == approx(e2)

def test_AUSM_gas_():
    ro, e, gamma, b = 2, 300, 1.2, 0.0001
    p = roe_to_p(ro, e, gamma, b)
    c = rop_to_csound(ro, p, gamma, b)
    u = 0.0

    f1, f2, f3 = AUSM_gas_(
        p1=p, 
        ro1=ro, 
        u1=u, 
        e1=e, 
        c1=c,
        p2=p, 
        ro2=ro, 
        u2=u, 
        e2=e, 
        c2=c,
        vbi=u)
    assert f1 == approx(0)
    assert f2 == approx(0)
    assert f3 == approx(0)


if __name__ == "__main__":
    ro, e, gamma, b = 1, 2, 1.2, 0.001
    p = roe_to_p(ro, e, gamma, b)
    test_AUSM_gas_()