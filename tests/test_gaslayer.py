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
from gaslayer import GasLayer, GasEOS


@pytest.fixture(scope="function")
def tube(request):
    return Tube([1,2,3], [0.1, 0.2, 0.15])

@pytest.fixture(scope="function")
def gasEOS(request):
    return GasEOS(1.4, 0.01)

@pytest.fixture(scope="function")
def gasLayer(request, tube, gasEOS):
    return GasLayer(n_cells=7, tube=tube, gasEOS=gasEOS)

def test_create_GasEOS(gasEOS):
    assert gasEOS is not None

def test_create_Tube(tube):
    assert tube is not None

def test_create_gasLayer(gasLayer):
    assert gasLayer is not None