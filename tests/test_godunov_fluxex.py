import os
import sys
from math import *
wd = os.path.abspath(__file__) # os.path.dirname(os.path.dirname(os.getcwd()))
wd = os.path.dirname(os.path.dirname(wd))
sys.path.append(wd+"\\src\\include\\godunov")
sys.path.append(wd+"\\src")
import numpy as np
import pytest
from pytest import approx
from godunov_python import Godunov_fluxes_python, Godunov_fluxes_cython, flux_wall_border,flux_wall_border_cython
from gaslayer import rop_to_csound, rop_to_e, roe_to_p

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
def get_fluxes_foo(**kwargs):
    # return Godunov_gas_(p1, ro1, u1, e1, c1, p2, ro2, u2, e2, c2, vbi)
    return Godunov_fluxes_cython(**kwargs)
    # return Godunov_fluxes_python(**kwargs)

def get_fluxes_border_foo(**kwargs):
    return flux_wall_border_cython(**kwargs)
    # return flux_wall_border(**kwargs)


def get_rnd_values(n=1000):
    result = []
    for _ in range(n):
        roo = np.random.uniform(0.1, 1000) 
        e = np.random.uniform(0.1, 10000)
        gamma = np.random.uniform(1.01, 1.9)
        b = np.random.uniform(0.00001, 0.0001)
        pp = roe_to_p(roo, e, gamma, b)
        c = rop_to_csound(roo, pp, gamma, b) 
        u=np.random.uniform(-100, 100)
        vbi=np.random.uniform(-100, 100)
        def e_foo(p, ro):
            return p*(1/ro-b)/(gamma-1)
        result.append((roo, e, pp, c, u, vbi, gamma, b, e_foo)) 
    return result

@pytest.fixture(scope="function")
def random_values4fluxes(request):
    return get_rnd_values()

@pytest.fixture(scope="function")
def random_values4fluxes2(request):
    return get_rnd_values()

@pytest.fixture(scope="function")
def random_values4fluxes3(request):
    return get_rnd_values()

def test_Godunov_gas_zero_flux_same_cells_f1(random_values4fluxes):
    for ro, e, p, c, u, vbi, gamma, b, e_foo in random_values4fluxes:
        f1, f2, f3 = get_fluxes_foo(
            p1=p, ro1=ro, u1=u, c1=c,
            p2=p, ro2=ro, u2=u, c2=c,
            vbi=u, e_foo=e_foo, gamma=gamma)
        assert f1 == approx(0, abs=1e-9)


def test_Godunov_gas_zero_flux_same_cells_f2(random_values4fluxes):
    for ro, e, p, c, u, vbi, gamma, b, e_foo in random_values4fluxes:
        f1, f2, f3 = get_fluxes_foo(
            p1=p, ro1=ro, u1=u, e1=e, c1=c,
            p2=p, ro2=ro, u2=u, e2=e, c2=c,
            vbi=u, e_foo=e_foo, gamma=gamma)

        assert f2 == approx(p)


def test_Godunov_gas_zero_flux_same_cells_f3(random_values4fluxes):
    for ro, e, p, c, u, vbi, gamma, b, e_foo in random_values4fluxes:
        f1, f2, f3 = get_fluxes_foo(
            p1=p, ro1=ro, u1=u, e1=e, c1=c,
            p2=p, ro2=ro, u2=u, e2=e, c2=c,
            vbi=u, e_foo=e_foo, gamma=gamma)
        assert f3 == approx(p*u)

def test_Godunov_gas_zero_flux_zero_u_zero_vbi_f1(random_values4fluxes):
    for ro, e, p, c, u, vbi, gamma, b, e_foo in random_values4fluxes:
        u = 0
        f1, f2, f3 = get_fluxes_foo(
            p1=p, ro1=ro, u1=u, e1=e, c1=c,
            p2=p, ro2=ro, u2=u, e2=e, c2=c,
            vbi=u, e_foo=e_foo, gamma=gamma)

        assert f1 == approx(0)

def test_Godunov_gas_zero_flux_zero_u_zero_vbi_f2(random_values4fluxes):
    for ro, e, p, c, u, vbi, gamma, b, e_foo in random_values4fluxes:
        u = 0
        f1, f2, f3 = get_fluxes_foo(
            p1=p, ro1=ro, u1=u, e1=e, c1=c,
            p2=p, ro2=ro, u2=u, e2=e, c2=c,
            vbi=u, e_foo=e_foo, gamma=gamma)
        assert f2 == approx(p)

def test_Godunov_gas_zero_flux_zero_u_zero_vbi_f3(random_values4fluxes):
    for ro, e, p, c, u, vbi, gamma, b, e_foo in random_values4fluxes:
        u = 0
        f1, f2, f3 = get_fluxes_foo(
            p1=p, ro1=ro, u1=u, e1=e, c1=c,
            p2=p, ro2=ro, u2=u, e2=e, c2=c,
            vbi=u, e_foo=e_foo, gamma=gamma)
        assert f3 == approx(0)

def test_Godunov_gas_border_wall_vbiZero_f1(random_values4fluxes):
    for ro, e, p, c, u, vbi, gamma, b, e_foo in random_values4fluxes:
        f1, f2, f3 = get_fluxes_border_foo(
            p=p, ro=ro, u=u, c=c,
            left_border=True,
            vbi=0, e_foo=e_foo, gamma=gamma)

        assert f1 == approx(0, abs=1e-5)
        f1, f2, f3 = get_fluxes_border_foo(
            p=p, ro=ro, u=u, c=c,
            left_border=False,
            vbi=0, e_foo=e_foo, gamma=gamma)
        assert f1 == approx(0, abs=1e-5)


def test_Godunov_gas_border_wall_vbiZero_f2(random_values4fluxes):
    for ro, e, p, c, u, vbi, gamma, b, e_foo in random_values4fluxes:
        f11, f12, f13 = get_fluxes_border_foo(
            p=p, ro=ro, u=u, c=c,
            left_border=True,
            vbi=0, e_foo=e_foo, gamma=gamma)

        f21, f22, f23 = get_fluxes_border_foo(
            p=p, ro=ro, u=-u, c=c,
            left_border=False,
            vbi=0, e_foo=e_foo, gamma=gamma)
        assert f22 == approx(f12, abs=1e-5)

def test_Godunov_gas_border_wall_vbiZero_f3(random_values4fluxes):
    for ro, e, p, c, u, vbi, gamma, b, e_foo in random_values4fluxes:
        f1, f2, f3 = get_fluxes_border_foo(
            p=p, ro=ro, u=u, c=c,
            left_border=True,
            vbi=0, e_foo=e_foo, gamma=gamma)

        assert f3 == approx(0, abs=1e-5)
        f1, f2, f3 = get_fluxes_border_foo(
            p=p, ro=ro, u=u, c=c,
            left_border=False,
            vbi=0, e_foo=e_foo, gamma=gamma)

        assert f3 == approx(0, abs=1e-5)


def test_Godunov_gas_border_wall_f1(random_values4fluxes):
    for ro, e, p, c, u, vbi, gamma, b, e_foo in random_values4fluxes:
        f1, f2, f3 = get_fluxes_border_foo(
            p=p, ro=ro, u=u, c=c,
            left_border=True,
            vbi=vbi, e_foo=e_foo, gamma=gamma)

        assert f1 == approx(0, abs=1e-5)
        f1, f2, f3 = get_fluxes_border_foo(
            p=p, ro=ro, u=u, c=c,
            left_border=False,
            vbi=vbi, e_foo=e_foo, gamma=gamma)

        assert f1 == approx(0, abs=1e-5)

def test_Godunov_gas_border_wall_f2(random_values4fluxes):
    for ro, e, p, c, u, vbi, gamma, b, e_foo in random_values4fluxes:
        f11, f12, f13 = get_fluxes_border_foo(
            p=p, ro=ro, u=u, c=c,
            left_border=True,
            vbi=u, e_foo=e_foo, gamma=gamma)

        f21, f22, f23 = get_fluxes_border_foo(
            p=p, ro=ro, u=u, c=c,
            left_border=False,
            vbi=u, e_foo=e_foo, gamma=gamma)
        assert f22 == approx(f12, abs=1e-5)

def test_Godunov_gas_border_wall_f3(random_values4fluxes):
    for ro, e, p, c, u, vbi, gamma, b, e_foo in random_values4fluxes:
        f11, f12, f13 = get_fluxes_border_foo(
            p=p, ro=ro, u=u, c=c,
            left_border=True,
            vbi=vbi, e_foo=e_foo, gamma=gamma)

        f21, f22, f23 = get_fluxes_border_foo(
            p=p, ro=ro, u=-u, c=c,
            left_border=False,
            vbi=-vbi, e_foo=e_foo, gamma=gamma)
        assert f23 == approx(-f13, abs=1e-5)

def test_Godunov_gas_simmetric_flux_f1(random_values4fluxes, random_values4fluxes2):
    for (ro1, e1, p1, c1, u1, vbi, gamma, b, e_foo), (ro2, e2, p2, c2, u2, vbi2, gamma2, b2, e_foo2) in zip(random_values4fluxes, random_values4fluxes2):
        f11, f12, f13 = get_fluxes_foo(
            p1=p1, ro1=ro1, u1=u1, e1=e1, c1=c1,
            p2=p2, ro2=ro2, u2=u2, e2=e2, c2=c2,
            vbi=vbi, e_foo=e_foo, gamma=gamma)
        f21, f22, f23 = get_fluxes_foo(
            p1=p2, ro1=ro2, u1=-u2, e1=e2, c1=c2, 
            p2=p1, ro2=ro1, u2=-u1, e2=e1, c2=c1,
            vbi=-vbi, e_foo=e_foo, gamma=gamma)
        assert f11 == approx(-f21)

def test_Godunov_gas_simmetric_flux_f2(random_values4fluxes, random_values4fluxes2):
    for (ro1, e1, p1, c1, u1, vbi, gamma, b, e_foo), (ro2, e2, p2, c2, u2, vbi2, gamma2, b2, e_foo2) in zip(random_values4fluxes, random_values4fluxes2):
        f11, f12, f13 = get_fluxes_foo(
            p1=p1, ro1=ro1, u1=u1, e1=e1, c1=c1,
            p2=p2, ro2=ro2, u2=u2, e2=e2, c2=c2,
            vbi=vbi, e_foo=e_foo, gamma=gamma)

        f21, f22, f23 = get_fluxes_foo(
            p1=p2, ro1=ro2, u1=-u2, e1=e2, c1=c2, 
            p2=p1, ro2=ro1, u2=-u1, e2=e1, c2=c1,
            vbi=-vbi, e_foo=e_foo, gamma=gamma)
        assert f12 == approx(f22)

def test_Godunov_gas_simmetric_flux_f3(random_values4fluxes, random_values4fluxes2):
    for (ro1, e1, p1, c1, u1, vbi, gamma, b, e_foo), (ro2, e2, p2, c2, u2, vbi2, gamma2, b2, e_foo2) in zip(random_values4fluxes, random_values4fluxes2):
        f11, f12, f13 = get_fluxes_foo(
            p1=p1, ro1=ro1, u1=u1, e1=e1, c1=c1,
            p2=p2, ro2=ro2, u2=u2, e2=e2, c2=c2,
            vbi=vbi, e_foo=e_foo, gamma=gamma)

        f21, f22, f23 = get_fluxes_foo(
            p1=p2, ro1=ro2, u1=-u2, e1=e2, c1=c2, 
            p2=p1, ro2=ro1, u2=-u1, e2=e1, c2=c1,
            vbi=-vbi, e_foo=e_foo, gamma=gamma)
        assert f13 == approx(-f23)

def test_Godunov_gas_simmetric_flux_zero_vbi_f1(random_values4fluxes, random_values4fluxes2):
    for (ro1, e1, p1, c1, u1, vbi, gamma, b, e_foo), \
        (ro2, e2, p2, c2, u2, vbi2, gamma2, b2, e_foo2) in zip(random_values4fluxes, random_values4fluxes2):
        
        vbi = 0
        f11, f12, f13 = get_fluxes_foo(
            p1=p1, ro1=ro1, u1=u1, e1=e1, c1=c1,
            p2=p2, ro2=ro2, u2=u2, e2=e2, c2=c2,
            vbi=vbi, e_foo=e_foo, gamma=gamma)

        f21, f22, f23 = get_fluxes_foo(
            p1=p2, ro1=ro2, u1=-u2, e1=e2, c1=c2, 
            p2=p1, ro2=ro1, u2=-u1, e2=e1, c2=c1,
            vbi=-vbi, e_foo=e_foo, gamma=gamma)
        assert f11 == approx(-f21)

def test_Godunov_gas_simmetric_flux_zero_vbi_f2(random_values4fluxes, random_values4fluxes2):
    for (ro1, e1, p1, c1, u1, vbi, gamma, b, e_foo), \
        (ro2, e2, p2, c2, u2, vbi2, gamma2, b2, e_foo2) in zip(random_values4fluxes, random_values4fluxes2):
        
        vbi = 0
        f11, f12, f13 = get_fluxes_foo(
            p1=p1, ro1=ro1, u1=u1, e1=e1, c1=c1,
            p2=p2, ro2=ro2, u2=u2, e2=e2, c2=c2,
            vbi=vbi, e_foo=e_foo, gamma=gamma)

        f21, f22, f23 = get_fluxes_foo(
            p1=p2, ro1=ro2, u1=-u2, e1=e2, c1=c2, 
            p2=p1, ro2=ro1, u2=-u1, e2=e1, c2=c1,
            vbi=-vbi, e_foo=e_foo, gamma=gamma)
        assert f12 == approx(f22)

def test_Godunov_gas_simmetric_flux_zero_vbi_f3(random_values4fluxes, random_values4fluxes2):
    for (ro1, e1, p1, c1, u1, vbi, gamma, b, e_foo), \
        (ro2, e2, p2, c2, u2, vbi2, gamma2, b2, e_foo2) in zip(random_values4fluxes, random_values4fluxes2):
        
        vbi = 0
        f11, f12, f13 = get_fluxes_foo(
            p1=p1, ro1=ro1, u1=u1, e1=e1, c1=c1,
            p2=p2, ro2=ro2, u2=u2, e2=e2, c2=c2,
            vbi=vbi,gamma=gamma, e_foo=e_foo)
        f21, f22, f23 = get_fluxes_foo(
            p1=p2, ro1=ro2, u1=-u2, e1=e2, c1=c2, 
            p2=p1, ro2=ro1, u2=-u1, e2=e1, c2=c1,
            vbi=-vbi, e_foo=e_foo, gamma=gamma)
        assert f13 == approx(-f23)

def test_Godunov_gas_simmetric_flux_TREE_cells_f1(random_values4fluxes, random_values4fluxes2, random_values4fluxes3):
    for (ro1, e1, p1, c1, u1, vbi, gamma, b, e_foo), \
        (ro2, e2, p2, c2, u2, vbi2, gamma2, b2, e_foo2), \
        (ro3, e3, p3, c3, u3, vbi3, gamma3, b3, e_foo3) in zip(random_values4fluxes, random_values4fluxes2, random_values4fluxes3):
        
        f12_1, f12_2, f12_3 = get_fluxes_foo(
            p1=p1, ro1=ro1, u1=u1, e1=e1, c1=c1,
            p2=p2, ro2=ro2, u2=u2, e2=e2, c2=c2,
            vbi=vbi, e_foo=e_foo, gamma=gamma)

        f23_1, f23_2, f23_3 = get_fluxes_foo(
            p1=p2, ro1=ro2, u1=u2, e1=e2, c1=c2, 
            p2=p3, ro2=ro3, u2=u3, e2=e3, c2=c3,
            vbi=vbi2, e_foo=e_foo, gamma=gamma)

        f32_1, f32_2, f32_3 = get_fluxes_foo(
            p1=p3, ro1=ro3, u1=-u3, e1=e3, c1=c3,
            p2=p2, ro2=ro2, u2=-u2, e2=e2, c2=c2,
            vbi=-vbi2, e_foo=e_foo, gamma=gamma)

        f21_1, f21_2, f21_3 = get_fluxes_foo(
            p1=p2, ro1=ro2, u1=-u2, e1=e2, c1=c2, 
            p2=p1, ro2=ro1, u2=-u1, e2=e1, c2=c1,
            vbi=-vbi, e_foo=e_foo, gamma=gamma)

        assert f23_1 - f12_1 == approx(f21_1-f32_1)

def test_Godunov_gas_simmetric_flux_TREE_cells_f2(random_values4fluxes, random_values4fluxes2, random_values4fluxes3):
    for (ro1, e1, p1, c1, u1, vbi, gamma, b, e_foo), \
        (ro2, e2, p2, c2, u2, vbi2, gamma2, b2, e_foo2), \
        (ro3, e3, p3, c3, u3, vbi3, gamma3, b3, e_foo3) in zip(random_values4fluxes, random_values4fluxes2, random_values4fluxes3):
        
        f12_1, f12_2, f12_3 = get_fluxes_foo(
            p1=p1, ro1=ro1, u1=u1, e1=e1, c1=c1,
            p2=p2, ro2=ro2, u2=u2, e2=e2, c2=c2,
            vbi=vbi, e_foo=e_foo, gamma=gamma)

        f23_1, f23_2, f23_3 = get_fluxes_foo(
            p1=p2, ro1=ro2, u1=u2, e1=e2, c1=c2, 
            p2=p3, ro2=ro3, u2=u3, e2=e3, c2=c3,
            vbi=vbi2, e_foo=e_foo, gamma=gamma)

        f32_1, f32_2, f32_3 = get_fluxes_foo(
            p1=p3, ro1=ro3, u1=-u3, e1=e3, c1=c3,
            p2=p2, ro2=ro2, u2=-u2, e2=e2, c2=c2,
            vbi=-vbi2, e_foo=e_foo, gamma=gamma)

        f21_1, f21_2, f21_3 = get_fluxes_foo(
            p1=p2, ro1=ro2, u1=-u2, e1=e2, c1=c2, 
            p2=p1, ro2=ro1, u2=-u1, e2=e1, c2=c1,
            vbi=-vbi, e_foo=e_foo, gamma=gamma)

        assert f23_2 - f12_2 == approx(-f21_2+f32_2)

def test_Godunov_gas_simmetric_flux_TREE_cells_f3(random_values4fluxes, random_values4fluxes2, random_values4fluxes3):
    for (ro1, e1, p1, c1, u1, vbi, gamma, b, e_foo), \
        (ro2, e2, p2, c2, u2, vbi2, gamma2, b2, e_foo2), \
        (ro3, e3, p3, c3, u3, vbi3, gamma3, b3, e_foo3) in zip(random_values4fluxes, random_values4fluxes2, random_values4fluxes3):
        
        f12_1, f12_2, f12_3 = get_fluxes_foo(
            p1=p1, ro1=ro1, u1=u1, e1=e1, c1=c1,
            p2=p2, ro2=ro2, u2=u2, e2=e2, c2=c2,
            vbi=vbi, e_foo=e_foo, gamma=gamma)

        f23_1, f23_2, f23_3 = get_fluxes_foo(
            p1=p2, ro1=ro2, u1=u2, e1=e2, c1=c2, 
            p2=p3, ro2=ro3, u2=u3, e2=e3, c2=c3,
            vbi=vbi2, e_foo=e_foo, gamma=gamma)

        f32_1, f32_2, f32_3 = get_fluxes_foo(
            p1=p3, ro1=ro3, u1=-u3, e1=e3, c1=c3,
            p2=p2, ro2=ro2, u2=-u2, e2=e2, c2=c2,
            vbi=-vbi2, e_foo=e_foo, gamma=gamma)

        f21_1, f21_2, f21_3 = get_fluxes_foo(
            p1=p2, ro1=ro2, u1=-u2, e1=e2, c1=c2, 
            p2=p1, ro2=ro1, u2=-u1, e2=e1, c2=c1,
            vbi=-vbi, e_foo=e_foo, gamma=gamma)

        assert f23_3 - f12_3 == approx(f21_3-f32_3)





    


if __name__ == "__main__":
    ro, e, gamma, b = 1, 2, 1.2, 0.001
    p = roe_to_p(ro, e, gamma, b)
    
