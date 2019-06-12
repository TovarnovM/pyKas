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


@pytest.fixture(scope="function")
def tube(request):
    return Tube([1,2,3], [0.1, 0.2, 0.15])

@pytest.fixture(scope="function")
def flux_calculator(request):
    return GasFluxCalculator(1,1,1)

@pytest.fixture(scope="function")
def gasEOS(request):
    return GasEOS(1.4, 0.01)

@pytest.fixture(scope="function")
def gasLayer(request, tube, gasEOS, flux_calculator, grid_strecher):
    return GasLayer(n_cells=7, tube=tube, gasEOS=gasEOS, 
        flux_calculator=flux_calculator, grid_strecher=grid_strecher, n_qs=3)

@pytest.fixture(scope="function")
def gasLayer2(request, gasLayer):
    return gasLayer

@pytest.fixture(scope="function")
def grid_strecher(request):
    return GridStrecher(1)

def test_create_GasEOS(gasEOS):
    assert gasEOS is not None

def test_create_Tube(tube):
    assert tube is not None

def test_create_GridStrecher(grid_strecher):
    assert grid_strecher is not None

def test_create_gasLayer(gasLayer, gasLayer2):
    assert gasLayer is not None
    assert gasLayer2 is not None

def test_grid_strecher_type3():
    gs3 = GridStrecher(strech_type=3, index_anchor=-1)
    gs1 = GridStrecher(strech_type=1, index_anchor=-1)
    
    xs3 = np.linspace(-10, 10)
    xs1 = np.linspace(-10, 10)
    vs3 = np.zeros_like(xs3)
    vs1 = np.zeros_like(xs1)

    for i in range(1000):
        vl, vr = np.random.uniform(-100, 100, size=2)
        gs3.fill_Vs_borders_anchor(vl, vr, xs3, vs3)
        gs1.fill_Vs_borders_proportional(vl, vr, xs1, vs1)
        assert vs3 == approx(vs1)

def test_grid_strecher_type3_2():
    gs3 = GridStrecher(strech_type=3, index_anchor=50)
    gs1 = GridStrecher(strech_type=1, index_anchor=-1)
    
    xs3 = np.linspace(-10, 10, 100)
    xs1 = np.linspace(-10, 10, 100)
    vs3 = np.zeros_like(xs3)
    vs1 = np.zeros_like(xs1)

    for i in range(1000):
        gs3.index_anchor = np.random.randint(1,99)
        vl, vr = np.random.uniform(-100, 100, size=2)
        gs3.fill_Vs_borders_anchor(vl, vr, xs3, vs3)
        gs1.fill_Vs_borders_proportional(vl, vr, xs1, vs1)
        assert vs3 != approx(vs1)
        assert vs3[gs3.index_anchor] == 0
        assert vs3[0] == approx(vs1[0])
        assert vs3[-1] == approx(vs1[-1])

def test_grid_strecher_type3_3():
    gs3 = GridStrecher(strech_type=3, index_anchor=50)
    gs1 = GridStrecher(strech_type=1, index_anchor=-1)
    n=5
    xs3 = np.linspace(-10, 10, n)
    xs1 = np.linspace(-10, 10, n)
    vs3 = np.zeros_like(xs3)
    vs1 = np.zeros_like(xs1)

    for i in range(1000):
        
        ai = np.random.randint(1,n-1)
        gs3.index_anchor = ai
        vl, vr = np.random.uniform(-100, 100, size=2)
        gs3.fill_Vs_borders_anchor(vl, vr, xs3, vs3)
        gs1.fill_Vs_borders_proportional(0, vr, xs1[ai:], vs1[ai:])
        gs1.fill_Vs_borders_proportional(vl, 0, xs1[:ai+1], vs1[:ai+1])
        ai = 1
        assert vs3 == approx(vs1)
        

def test_GridStrecher_init_regular(grid_strecher):
    for _ in range(33):
        v1,v2 = np.random.uniform(-9999,9999,2)
        n = np.random.randint(2,99)

        dummy = np.empty(n)
        grid_strecher.init_regular(v1, v2, dummy)
        assert dummy == approx(np.linspace(v1, v2, n))
	