import sys
import os
wd = os.path.abspath(__file__) # os.path.dirname(os.path.dirname(os.getcwd()))
wd = os.path.dirname(os.path.dirname(wd))
sys.path.append(wd+"\\src\\")
from gaslayer import foo, q_to_roue, roue_to_q_, roe_to_p, rop_to_e, AUSM_gas_, rop_to_csound
import numpy as np
from pytest import approx
from math import *

def test_cimport():
    a = foo()
    assert a == 3

def test_q_to_roue():
    for _ in range(1000):
        q = np.random.uniform(0.1, 10000, size=3)
        ro, u, e = q_to_roue(q)
        assert ro == approx(q[0])
        assert u == approx(q[1]/q[0])
        assert e == approx(q[2]/q[0] - 0.5 * u * u)

def test_roue_to_q_():
    for _ in range(1000):
        q = np.array([0,0,0.0])
        ro, u, e = np.random.uniform(0.1, 10000), np.random.uniform(-1000, 1000), np.random.uniform(0.1, 10000)
        roue_to_q_(ro, u, e, q)
        assert ro == approx(q[0])
        assert u == approx(q[1]/q[0])
        assert e == approx(q[2]/q[0] - 0.5 * u * u)

def test_roe_to_p():
    for _ in range(1000):
        ro, e, gamma, b = np.random.uniform(0.1, 10000), np.random.uniform(0.1, 10000), np.random.uniform(1.01, 1.9), np.random.uniform(0.00001, 0.0001)
        p = roe_to_p(ro, e, gamma, b)

        e2 = rop_to_e(ro, p , gamma, b)
        assert e == approx(e2)
        assert e > 0
        assert p > 0

def test_csound():
    for _ in range(1000):
        ro, e, gamma, b = np.random.uniform(0.1, 10000), np.random.uniform(0.1, 10000), np.random.uniform(1.01, 1.9), np.random.uniform(0.00001, 0.0001)
        p = roe_to_p(ro, e, gamma, b)
        c = rop_to_csound(ro, p, gamma, b)  
        assert c > 0

def test_AUSM_gas_zero_flux_same_cells_f1():
    for i in range(1000):
        ro, e, gamma, b = np.random.uniform(0.1,1000), np.random.uniform(0.1,1000), np.random.uniform(1.01,1.7), np.random.uniform(0.0001,0.01)
        p = roe_to_p(ro, e, gamma, b)
        c = rop_to_csound(ro, p, gamma, b)
        u = np.random.uniform(-100,100)

        f1, f2, f3 = AUSM_gas_python(
            p1=p, ro1=ro, u1=u, e1=e, c1=c,
            p2=p, ro2=ro, u2=u, e2=e, c2=c,
            vbi=u)
        assert f1 == approx(0)


def test_AUSM_gas_zero_flux_same_cells_f2():
    for i in range(1000):
        ro, e, gamma, b = np.random.uniform(0.1,1000), np.random.uniform(0.1,1000), np.random.uniform(1.01,1.7), np.random.uniform(0.0001,0.01)
        p = roe_to_p(ro, e, gamma, b)
        c = rop_to_csound(ro, p, gamma, b)
        u = np.random.uniform(-100,100)

        f1, f2, f3 = AUSM_gas_python(
            p1=p, ro1=ro, u1=u, e1=e, c1=c,
            p2=p, ro2=ro, u2=u, e2=e, c2=c,
            vbi=u)

        assert f2 == approx(0)


def test_AUSM_gas_zero_flux_same_cells_f3():
    for i in range(1000):
        ro, e, gamma, b = np.random.uniform(0.1,1000), np.random.uniform(0.1,1000), np.random.uniform(1.01,1.7), np.random.uniform(0.0001,0.01)
        p = roe_to_p(ro, e, gamma, b)
        c = rop_to_csound(ro, p, gamma, b)
        u = np.random.uniform(-100,100)

        f1, f2, f3 = AUSM_gas_python(
            p1=p, ro1=ro, u1=u, e1=e, c1=c,
            p2=p, ro2=ro, u2=u, e2=e, c2=c,
            vbi=u)

        assert f3 == approx(0)

def test_AUSM_gas_zero_flux_zero_u_zero_vbi_f1():
    for i in range(1000):
        ro, e, gamma, b = np.random.uniform(0.1,1000), np.random.uniform(0.1,1000), np.random.uniform(1.01,1.7), np.random.uniform(0.0001,0.01)
        p = roe_to_p(ro, e, gamma, b)
        c = rop_to_csound(ro, p, gamma, b)
        u = 0

        f1, f2, f3 = AUSM_gas_python(
            p1=p, ro1=ro, u1=u, e1=e, c1=c,
            p2=p, ro2=ro, u2=u, e2=e, c2=c,
            vbi=u)

        assert f1 == approx(0)

def test_AUSM_gas_zero_flux_zero_u_zero_vbi_f2():
    for i in range(1000):
        ro, e, gamma, b = np.random.uniform(0.1,1000), np.random.uniform(0.1,1000), np.random.uniform(1.01,1.7), np.random.uniform(0.0001,0.01)
        p = roe_to_p(ro, e, gamma, b)
        c = rop_to_csound(ro, p, gamma, b)
        u = 0

        f1, f2, f3 = AUSM_gas_python(
            p1=p, ro1=ro, u1=u, e1=e, c1=c,
            p2=p, ro2=ro, u2=u, e2=e, c2=c,
            vbi=u)

        assert f2 == approx(0)

def test_AUSM_gas_zero_flux_zero_u_zero_vbi_f3():
    for i in range(1000):
        ro, e, gamma, b = np.random.uniform(0.1,1000), np.random.uniform(0.1,1000), np.random.uniform(1.01,1.7), np.random.uniform(0.0001,0.01)
        p = roe_to_p(ro, e, gamma, b)
        c = rop_to_csound(ro, p, gamma, b)
        u = 0
        f1, f2, f3 = AUSM_gas_python(
            p1=p, ro1=ro, u1=u, e1=e, c1=c,
            p2=p, ro2=ro, u2=u, e2=e, c2=c,
            vbi=u)
        assert f3 == approx(0)

def test_AUSM_gas_simmetric_flux_f1():
    for i in range(1000):
        ro1, e1, gamma, b = np.random.uniform(0.1,1000), np.random.uniform(0.1,1000), np.random.uniform(1.01,1.7), np.random.uniform(0.0001,0.01)
        p1 = roe_to_p(ro1, e1, gamma, b)
        c1 = rop_to_csound(ro1, p1, gamma, b)
        u1 = np.random.uniform(-100,100)
        
        ro2, e2 = np.random.uniform(0.1,1000), np.random.uniform(0.1,1000)
        p2 = roe_to_p(ro2, e2, gamma, b)
        c2 = rop_to_csound(ro2, p2, gamma, b)
        u2 = np.random.uniform(-100,100)

        vbi = np.random.uniform(-100,100)

        f11, f12, f13 = AUSM_gas_python(
            p1=p1, ro1=ro1, u1=u1, e1=e1, c1=c1,
            p2=p2, ro2=ro2, u2=u2, e2=e2, c2=c2,
            vbi=vbi)

        f21, f22, f23 = AUSM_gas_python(
            p1=p2, ro1=ro2, u1=-u2, e1=e2, c1=c2, 
            p2=p1, ro2=ro1, u2=-u1, e2=e1, c2=c1,
            vbi=-vbi)
        assert f11 == approx(-f21)

def test_AUSM_gas_simmetric_flux_f2():
    for i in range(1000):
        ro1, e1, gamma, b = np.random.uniform(0.1,1000), np.random.uniform(0.1,1000), np.random.uniform(1.01,1.7), np.random.uniform(0.0001,0.01)
        p1 = roe_to_p(ro1, e1, gamma, b)
        c1 = rop_to_csound(ro1, p1, gamma, b)
        u1 = np.random.uniform(-100,100)

        ro2, e2 = np.random.uniform(0.1,1000), np.random.uniform(0.1,1000)
        p2 = roe_to_p(ro2, e2, gamma, b)
        c2 = rop_to_csound(ro2, p2, gamma, b)
        u2 = np.random.uniform(-100,100)

        vbi = np.random.uniform(-100,100)

        f11, f12, f13 = AUSM_gas_python(
            p1=p1, ro1=ro1, u1=u1, e1=e1, c1=c1,
            p2=p2, ro2=ro2, u2=u2, e2=e2, c2=c2,
            vbi=vbi)

        f21, f22, f23 = AUSM_gas_python(
            p1=p2, ro1=ro2, u1=-u2, e1=e2, c1=c2, 
            p2=p1, ro2=ro1, u2=-u1, e2=e1, c2=c1,
            vbi=-vbi)
        assert f12 == approx(-f22)

def test_AUSM_gas_simmetric_flux_f3():
    for i in range(1000):
        ro1, e1, gamma, b = np.random.uniform(0.1,1000), np.random.uniform(0.1,1000), np.random.uniform(1.01,1.7), np.random.uniform(0.0001,0.01)
        p1 = roe_to_p(ro1, e1, gamma, b)
        c1 = rop_to_csound(ro1, p1, gamma, b)
        u1 = np.random.uniform(-100,100)

        ro2, e2 = np.random.uniform(0.1,1000), np.random.uniform(0.1,1000)
        p2 = roe_to_p(ro2, e2, gamma, b)
        c2 = rop_to_csound(ro2, p2, gamma, b)
        u2 = np.random.uniform(-100,100)

        vbi = np.random.uniform(-100,100)

        f11, f12, f13 = AUSM_gas_python(
            p1=p1, ro1=ro1, u1=u1, e1=e1, c1=c1,
            p2=p2, ro2=ro2, u2=u2, e2=e2, c2=c2,
            vbi=vbi)

        f21, f22, f23 = AUSM_gas_python(
            p1=p2, ro1=ro2, u1=-u2, e1=e2, c1=c2, 
            p2=p1, ro2=ro1, u2=-u1, e2=e1, c2=c1,
            vbi=-vbi)
        assert f13 == approx(-f23)

def test_AUSM_gas_simmetric_flux_zero_vbi_f1():
    for i in range(1000):
        ro1, e1, gamma, b = np.random.uniform(0.1,1000), np.random.uniform(0.1,1000), np.random.uniform(1.01,1.7), np.random.uniform(0.0001,0.01)
        p1 = roe_to_p(ro1, e1, gamma, b)
        c1 = rop_to_csound(ro1, p1, gamma, b)
        u1 = np.random.uniform(-100,100)

        ro2, e2 = np.random.uniform(0.1,1000), np.random.uniform(0.1,1000)
        p2 = roe_to_p(ro2, e2, gamma, b)
        c2 = rop_to_csound(ro2, p2, gamma, b)
        u2 = np.random.uniform(-100,100)

        vbi = 0

        f11, f12, f13 = AUSM_gas_python(
            p1=p1, ro1=ro1, u1=u1, e1=e1, c1=c1,
            p2=p2, ro2=ro2, u2=u2, e2=e2, c2=c2,
            vbi=vbi)

        f21, f22, f23 = AUSM_gas_python(
            p1=p2, ro1=ro2, u1=-u2, e1=e2, c1=c2, 
            p2=p1, ro2=ro1, u2=-u1, e2=e1, c2=c1,
            vbi=-vbi)
        assert f11 == approx(-f21)

def test_AUSM_gas_simmetric_flux_zero_vbi_f2():
    for i in range(1000):
        ro1, e1, gamma, b = np.random.uniform(0.1,1000), np.random.uniform(0.1,1000), np.random.uniform(1.01,1.7), np.random.uniform(0.0001,0.01)
        p1 = roe_to_p(ro1, e1, gamma, b)
        c1 = rop_to_csound(ro1, p1, gamma, b)
        u1 = np.random.uniform(-100,100)

        ro2, e2 = np.random.uniform(0.1,1000), np.random.uniform(0.1,1000)
        p2 = roe_to_p(ro2, e2, gamma, b)
        c2 = rop_to_csound(ro2, p2, gamma, b)
        u2 = np.random.uniform(-100,100)

        vbi = 0

        f11, f12, f13 = AUSM_gas_python(
            p1=p1, ro1=ro1, u1=u1, e1=e1, c1=c1,
            p2=p2, ro2=ro2, u2=u2, e2=e2, c2=c2,
            vbi=vbi)

        f21, f22, f23 = AUSM_gas_python(
            p1=p2, ro1=ro2, u1=-u2, e1=e2, c1=c2, 
            p2=p1, ro2=ro1, u2=-u1, e2=e1, c2=c1,
            vbi=-vbi)
        assert f12 == approx(-f22)

def test_AUSM_gas_simmetric_flux_zero_vbi_f3():
    for i in range(1000):
        ro1, e1, gamma, b = np.random.uniform(0.1,1000), np.random.uniform(0.1,1000), np.random.uniform(1.01,1.7), np.random.uniform(0.0001,0.01)
        p1 = roe_to_p(ro1, e1, gamma, b)
        c1 = rop_to_csound(ro1, p1, gamma, b)
        u1 = np.random.uniform(-100,100)

        ro2, e2 = np.random.uniform(0.1,1000), np.random.uniform(0.1,1000)
        p2 = roe_to_p(ro2, e2, gamma, b)
        c2 = rop_to_csound(ro2, p2, gamma, b)
        u2 = np.random.uniform(-100,100)

        vbi = 0

        f11, f12, f13 = AUSM_gas_python(
            p1=p1, ro1=ro1, u1=u1, e1=e1, c1=c1,
            p2=p2, ro2=ro2, u2=u2, e2=e2, c2=c2,
            vbi=vbi)

        f21, f22, f23 = AUSM_gas_python(
            p1=p2, ro1=ro2, u1=-u2, e1=e2, c1=c2, 
            p2=p1, ro2=ro1, u2=-u1, e2=e1, c2=c1,
            vbi=-vbi)
        assert f13 == approx(-f23)



def AUSM_gas_python(p1, ro1, u1, e1, c1, p2, ro2, u2, e2, c2, vbi):
    r1=ro1
    r2=ro2
    H1 = e1 + 0.5*u1*u1 + p1/r1
    H2 = e2 + 0.5*u2*u2 + p2/r2

    cs = 0.5*(c1+c2)
    Mr1 = (u1-vbi)/cs
    Mr2 = (u2-vbi)/cs

    if abs(Mr1) >= 1.0 :
        M4p = 0.5*(Mr1+abs(Mr1))
        P5p = 0.5*(Mr1+abs(Mr1))/Mr1    
    else:
        M4p = 0.25*((Mr1+1.0)*(Mr1+1.0))*(1.0+2.0*0.25*(Mr1-1.0)*(Mr1-1.0))
        P5p = 0.25*((Mr1+1.0)*(Mr1+1.0))*((2.0-Mr1)+3.0*Mr1*0.25*(Mr1-1.0)*(Mr1-1.0))   

    if abs(Mr2) >= 1.0:
        M4m = 0.5*(Mr2-abs(Mr2))
        P5m = 0.5*(Mr2-abs(Mr2))/Mr2
    else:
        M4m = -0.25*((Mr2-1.0)*(Mr2-1.0))*(1.0+2.0*0.25*(Mr2+1.0)*(Mr2+1.0))
        P5m = 0.25*((Mr2-1.0)*(Mr2-1.0))*((2.0+Mr2)-3.0*Mr2*0.25*(Mr2+1.0)*(Mr2+1.0))
    
    # новое статья Сравнение схем с расщеплением потока для численного решения уравнений Эйлера сжимаемого газа
    # стр 6 внизу
    K_u=0.75
    p_u = -K_u*P5p*P5m*(r1+r2)

    Mrf = M4p + M4m
    pf = P5p*p1 + P5m*p2

    flux1 = 0.5*cs*(Mrf*(r1+r2)-abs(Mrf)*(r2-r1))
    flux2 = copysign(0.5*cs*(Mrf*(r1*u1+r2*u2)-abs(Mrf)*(r2*u2-r1*u1)) + pf, Mrf)
    flux3 = 0.5*cs*(Mrf*(r1*H1+r2*H2)-abs(Mrf)*(r2*H2-r1*H1)) + pf*vbi

    return flux1, flux2, flux3


if __name__ == "__main__":
    ro, e, gamma, b = 1, 2, 1.2, 0.001
    p = roe_to_p(ro, e, gamma, b)
    test_AUSM_gas_()