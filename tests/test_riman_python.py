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
from riman_python import get_distrs_to_time
from gaslayer import rop_to_csound, rop_to_e, roe_to_p
from copy import deepcopy


def get_values(n=1000):
    v = [
         ( 1,       0,      1,      0.125,  0,      0.1,    0.15  ),
         ( 0.445,   0.698,  3.528,  0.5,    0,      0.571,  0.15  ),       
         ( 1,       0,      1,      0.02,   3.55,   1,      0.15  ),
         ( 3.857,   0.920,  10.333, 1,      3.55,   1,      0.09  ),
         ( 1,       0,      0.5,    0.5,    0,      0.5,    0.42  ),
         ( 1,       0.5,    0.5,    0.5,    0.5,    0.5,    0.43  ),
         ( 1,       -1,     1,      0.9275, -1.0781,0.9,    0.18  ),
         ( 1,       0,      1000,   1,      0,      0.01,   0.012  ),
         ( 10,      2000,   500,    20,     0,      500,    0.012  ),
         ( 1,       -2,     0.4,    1,      2,      0.4,    0.15  )
        ]
    for ro_1, u_1, p_1, ro_2, u_2, p_2, t in v:
        yield {
            'p_1' : p_1,          # давление слева
            'ro_1' : ro_1,         # плотность слева   
            'u_1' : u_1,          # скорость слева      
            'p_2' : p_2,        # давление справа     
            'ro_2': ro_2,      # плотность справа    
            'u_2' : u_2,          # скорость справа        
            'p_0' : 0,          # параметр в ур-ии состояния        
            'gamma' : 1.4,      # параметр в ур-ии состояния          
            'c_0' : 0,          # параметр в ур-ии состояния     
            'eps_F':1e-6,       # точность определения решения           
            'n_iter_max':100,   # максимальное количество итераций при определении решения              
            'x0' : 0.0,         # положение разрыва в момент t=0        
            'ts': [t],        # времена, в которых нужно построить графики распределения параметров         
            'n': 100          # кол-во точек, в которых ищутся параметры волны разрежения         
        } 

@pytest.fixture(scope="function")
def random_values4fluxes(request):
    return list(get_values())

def test_simm_get_distrs_to_time(random_values4fluxes):
    for init_cond in random_values4fluxes:
        t = init_cond['ts'][0]
        res1 = get_distrs_to_time(t, **init_cond)
        simm_ic = deepcopy(init_cond)
        simm_ic['p_1'] = init_cond['p_2']
        simm_ic['ro_1'] = init_cond['ro_2']
        simm_ic['u_1'] = -init_cond['u_2']
        simm_ic['p_2'] = init_cond['p_1']
        simm_ic['ro_2'] = init_cond['ro_1']
        simm_ic['u_2'] = -init_cond['u_1']
        res2 = get_distrs_to_time(t, **simm_ic)
        assert res1['xs'] == approx(-np.flip(res2['xs']))
        assert res1['ros'] == approx(np.flip(res2['ros']))
        assert res1['ps'] == approx(np.flip(res2['ps']))
        assert res1['us'] == approx(-np.flip(res2['us']))
        assert res1['es'] == approx(np.flip(res2['es']))
        assert res1['ms'] == approx(-np.flip(res2['ms']))