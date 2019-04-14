import json
import os
import sys
from math import *
wd = os.path.abspath(__file__) # os.path.dirname(os.path.dirname(os.getcwd()))
wd = os.path.dirname(os.path.dirname(wd))
sys.path.append(wd+"\\src\\include\\godunov")
sys.path.append(wd+"\\src")
import numpy as np

from gaslayer import Powder, PowderOvLayer, GridStrecher, GasFluxCalculator
from tube import Tube

d = 23 * 1e-3
barrel = Tube([0,1], [d, d])


with open('gpowders.json') as f:
    all_powders = json.load(f)

print(all_powders.keys())