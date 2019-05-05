import sys
import os
wd = os.path.abspath(__file__) # os.path.dirname(os.path.dirname(os.getcwd()))
wd = os.path.dirname(os.path.dirname(wd))
sys.path.append(wd+"\\src\\")
from tube import Tube, InterpXY
from pistonlayer import ElPistEOS
from math import pi, sqrt
from pytest import approx
import numpy as np
import random

def _main1():
    eos = ElPistEOS()
    print(eos)
    

def test_pe_synch():
    eos = ElPistEOS(zeroP=False)
    for i in range(10000):
        ro = np.random.uniform(0.75*eos.ro_0, 2.5*eos.ro_0)
        emax = eos.c_0*eos.c_0/(eos.gamma-1)
        e = np.random.uniform(0, emax)

        p = eos.get_p(ro, e)
        e1 = eos.get_e(ro, p)
        assert e == approx(e1, abs=100)

def test_p():
    eos = ElPistEOS()
    for i in range(10000):
        ro = np.random.uniform(0.75*eos.ro_0, 2.5*eos.ro_0)
        emax = eos.c_0*eos.c_0/(eos.gamma-1)*1000
        e = np.random.uniform(0, emax)

        p = eos.get_p(ro, e)
        assert p >= 0 

if __name__ == "__main__":
    _main1()